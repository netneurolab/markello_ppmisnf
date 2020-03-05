# -*- coding: utf-8 -*-
"""
Contains codebits and functions supporting various analyses
"""

import os.path as op

from mapalign import align
import numpy as np
import pandas as pd
from scipy import stats as sstats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import fdrcorrection

from ppmi_snf import directories, utils
import pypmi
from pypmi._info import VISITS


def get_zrand_mask(zrand, percentile=95):
    """
    Calculates stable regions of `zrand` based on `percentile`

    Parameters
    ----------
    zrand : (C, K, M) array_like
        Local similarity of clustering solutions in SNF parameter space
    percentile : [0, 100] float, optional
        Percentile of `zrand` at which to threshold. Default: 95

    Returns
    -------
    mask : (K, M) numpy.ndarray
        Boolean mask to subset stable regions of parameter space
    """

    # by default, use the 95%ile of ALL the z-rand scores (across all different
    # cluster numbers) to threshold the similarity matrices and extract
    # "stable" regions of hyperparameter space
    threshold = np.percentile(zrand, 95)
    masked = [utils.cluster_img_2d(z, threshold)[1] != 0 for z in zrand]

    # we'll use these stable regions as a mask for `fusion` / `labels` matrices
    mask = np.sum(masked, axis=0) > 0

    return mask


def get_embedding_variance(fusion, verbose=True):
    """
    Calculated variance explained by diffusion embedding of `fusion` networks

    Parameters
    ----------
    fusion : (A, N, N) array_like
        Fused matrices generated from SNF gridsearch

    Returns
    -------
    embedding : (N, C) numpy.ndarray
        Average scores for `N` patients across `C` embedding dimensions
    realigned : (A, N, C) numpy.ndarray
        Aligned scores of `N` patients across `C` embedding dimensions for `A`
        different embeddings
    """

    # we calculate embedding with consistent # of dimensions so we can average
    # across them (if left to be "random" we get different #s)
    embeddings, results = zip(*[
        utils.dme(net, n_components=10, return_result=True) for net in fusion
    ])

    # align embeddings w/generalized Procrustes and average across embeddings
    realigned, xfms = align.iterative_alignment([e for e in embeddings])
    embedding = np.mean(realigned, axis=0)

    # normalize lambdas based on sum of lambdas from all components and check
    # variance explained by first five components
    if verbose:
        lambdas = [res['lambdas'] for res in results]
        varexp = [np.sum(var[:5] / var.sum()) for var in lambdas]
        mvar, sdvar = np.mean(varexp) * 100, np.std(varexp, ddof=1) * 100
        print('Variance explained by 5 components: '
              f'{mvar:.2f}% +/- {sdvar:.2f}%')

    return embedding, realigned


def load_longitudinal_behavior(participants):
    """
    Loads longitudinal behavioral data

    Parameters
    ----------
    participants : list
        List of participant IDs for whom behavioral data should be loaded

    Returns
    -------
    behavior : pandas.DataFrame
        Longitudinal behavioral data for `participants`
    """

    # ensure this is a list (DataFrame.query doesn't work with arrays)
    participants = list(participants)

    # load ALL behavioral data (not just first visit!)
    # we'll use this to run some longitudinal models to see how the
    # clusters and embeddings we generate might relate to changes in
    # symptomatology over time`
    behavior = pypmi.load_behavior(directories.ppmi, measures='all')
    first = (behavior.drop_duplicates(['participant', 'visit'], 'first')
                     .reset_index(drop=True))
    last = (behavior.drop_duplicates(['participant', 'visit'], 'last')
                    .reset_index(drop=True))
    behavior = first.combine_first(last)
    behavior = behavior.query(f'participant in {participants}') \
                       .assign(visit=behavior['visit'].astype(VISITS))

    # generate a continuous "time" variable (visit_date - base_visit_date)
    # this is what regressions will be performed against
    base = behavior.groupby('participant').min()['date'].rename('base')
    behavior = pd.merge(behavior, base, on='participant')
    time = (behavior['date'] - behavior['base']) / np.timedelta64(1, 'Y')
    behavior = behavior.assign(time=time).dropna(subset=['visit'])

    return behavior


def run_univariate_anova(data, demographics, verbose=True, run_tukey=True):
    """
    Runs mass-univariate one-way ANOVAs comparing `data` across clusters

    Parameters
    ----------
    data : list of pandas.DataFrame
        Data provided as input to SNF
    demographics : pandas.DataFrame
        Demographic information with at least columns 'cluster'
    verbose : bool, optional

    Returns
    -------
    stats : pandas.DataFrame
        With columns 'datatype', 'feature', 'fval', and 'qval' for features
        differentiable across clusters after FDR-correction
    """

    # get datatypes, features, and concatenated data
    datatypes = np.repeat([
        'cortical thickness',
        'subcortical volume',
        'dopamine binding',
        'CSF assays',
        'clinical measures'
    ], [df.columns.size for df in data])
    features = np.hstack([df.columns for df in data])
    data = sstats.zscore(np.column_stack(data))

    # run mass-univariate anova b/w clusters to examine differentiable features
    clusters = np.asarray(demographics['cluster'])
    fstats, pvals = sstats.f_oneway(
        *(data[clusters == cl] for cl in np.unique(clusters))
    )

    # print most discriminable feature per datatype and run post-hoc Tukey test
    if verbose:
        print('Most discriminating feature for each datatype:')
        print('-' * 46)
        for dt in np.unique(datatypes):
            mask = datatypes == dt
            idx = fstats[mask].argmax()
            feat = features[mask][idx]
            print('{:<20} {}'.format(dt + ':', feat))

            if run_tukey:
                tukey = pairwise_tukeyhsd(data.T[mask][idx], clusters, 0.01)
                print()
                print(tukey.summary())
                print()

    # FDR-correct p-values and check # of differentiable features
    reject, qvals = fdrcorrection(pvals)
    if verbose:
        print(f'\n{reject.sum()} features differentiable across groups')

    # create output table with differentiable features
    idx = fstats[reject].argsort()[::-1]
    stats = pd.DataFrame(dict(datatype=datatypes[reject][idx],
                              feature=features[reject][idx],
                              fval=fstats[reject][idx],
                              qval=qvals[reject][idx]))

    return stats


def run_pdatrophy_anova(demographics, verbose=True, run_tukey=True):
    """
    Runs one-way ANOVA examining diff b/w `clusters` for PD-ICA atrophy scores

    Parameters
    ----------
    demographics : pandas.DataFrame
        With at least columns 'participant' and 'cluster'
    verbose : bool, optional

    Returns
    -------
    pdatrophy : pandas.DataFrame
        With columns 'participant', 'session', 'and 'atrophy'
    """

    # read in pre-calculated PD-ICA atrophy scores
    pdatrophy = pd.read_csv(op.join(directories.parcels, 'pdica_atrophy.csv'))

    # get only subjects who we have cluster labels for and grab first session
    subs = list(np.asarray(demographics['participant']))
    pdatrophy = pdatrophy.query(f'session == 1 & participant in {subs}')

    # z-score PD-ICA atrophy scores
    pdatrophy = pdatrophy.assign(atrophy=sstats.zscore(pdatrophy['atrophy']))

    # run one-way ANOVA to test for cluster differences in PD-ICA atrophy score
    clusters = np.asarray(demographics['cluster'])
    f, p = sstats.f_oneway(
        *(pdatrophy[clusters == cl]['atrophy'] for cl in np.unique(clusters))
    )
    if verbose:
        print(f'Atrophy z-score (cluster diffs): F = {f:.2f}, p = {p:.3f}')

    # post-hoc Tukey test for subgroup differences in PD-ICA atrophy scores
    if run_tukey:
        tukey = pairwise_tukeyhsd(np.asarray(pdatrophy['atrophy']), clusters)
        if verbose:
            print()
            print(tukey.summary())

    return pdatrophy
