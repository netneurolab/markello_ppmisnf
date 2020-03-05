# -*- coding: utf-8 -*-
"""
This script generates results and figures used to support supplementary
analyses demonstrating the robustness of the primary results
"""

import itertools
import os.path as op
import warnings

from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.linalg import orthogonal_procrustes
from scipy.ndimage import center_of_mass
from sklearn import metrics

from ppmi_snf import defaults, directories, structures, utils
import snf.metrics
from netneurotools import cluster
from pypmi.cluster import cluster_fereshtehnejad2017

from analysis import (get_embedding_variance, get_zrand_mask,
                      load_longitudinal_behavior,
                      run_pdatrophy_anova, run_univariate_anova)

# set seaborn, numpy, matplotlib, warnings options
sns.set(style='white', context='notebook', font_scale=1.5)
np.set_printoptions(suppress=True)
warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message='Calling np.sum')
plt.rcParams['svg.fonttype'] = 'none'
fontd = dict(size=36)

# miscellaneous settings
SAVE_FIGS = True  # whether to save figures
SEED = 1234  # random seed for reproducibility
np.random.seed(SEED)  # set random seed


def get_cluster_demographics(hdf):
    """
    Generates table with summary statistics for clusters

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file containing SNF gridsearch outputs

    Returns
    -------
    table : pandas.DataFrame
        Cluster demographics
    """

    # load demographics + cluster labels from hdf files
    demographics = hdf.load('/raw/pd_demographics')
    consensus = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/consensus')
    demographics = demographics.assign(cluster=pd.Categorical(consensus))

    # just a quick check...
    required_cols = [
        'date_diagnosis', 'date_birth', 'date_enroll', 'cluster', 'diagnosis',
        'gender', 'race', 'family_history', 'education'
    ]
    missing = np.setdiff1d(required_cols, demographics.columns)
    if len(missing) > 0:
        raise ValueError('Provided demographics dataframe missing required '
                         'columns: {}'.format(missing))

    # don't modify input
    demo = demographics.copy()

    # calculate some subject-level measures from existing data
    age_onset = np.asarray((demo['date_diagnosis'] - demo['date_birth'])
                           / np.timedelta64(1, 'Y'))
    symp_dur = np.asarray((demo['date_enroll'] - demo['date_diagnosis'])
                          / np.timedelta64(1, 'M'))

    # we need the full behavioral data frame to calculate UPDRS total scores
    behavior = load_longitudinal_behavior(demo.index)
    updrs = behavior.query('visit == "BL"') \
                    .set_index('participant') \
                    .loc[demo.index]
    updrs = updrs[['updrs_i', 'updrs_ii', 'updrs_iii', 'updrs_iv']].sum(axis=1)

    # generate some percent data for clusters
    num_clust = np.asarray(demo.groupby('cluster')['diagnosis'].count())
    num_male = np.asarray(demo.groupby('cluster')['gender']
                              .apply(lambda x: x[x == 'm'].count()))
    perc_male = (num_male / num_clust) * 100
    num_white = np.asarray(demo.groupby('cluster')['race']
                               .apply(lambda x: x[x == 'white'].count()))
    perc_white = (num_white / num_clust) * 100
    perc_famhist = 100 * np.asarray(
        demo.groupby('cluster')['family_history'].mean()
    )

    # get means + SDs for non-percent variables
    demo = demo.assign(age_onset=age_onset,
                       symptom_duration=symp_dur,
                       updrs=updrs)
    cols = ['age_onset', 'education', 'symptom_duration', 'updrs']
    means = demo.groupby('cluster').mean()[cols]
    stds = demo.groupby('cluster').std()[cols]

    table = pd.DataFrame(dict(
        age_onset_mean=means['age_onset'],
        age_onset_std=stds['age_onset'],
        perc_male=perc_male,
        perc_white=perc_white,
        education_mean=means['education'],
        education_std=stds['education'],
        symp_dur_mean=means['symptom_duration'],
        symp_dur_std=stds['symptom_duration'],
        perc_famhist=perc_famhist,
        updrs_mean=means['updrs'],
        updrs_std=stds['updrs']
    ))

    print(table)

    return table


def compare_fereshtehnejad2017(hdf):
    """
    Assesses how clustering criteria from [1]_ performs

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file containing SNF gridsearch outputs

    Returns
    -------
    feresh_labels : (N,) numpy.ndarray
        Cluster labels for `N` patients according to criteria from [1]_

    References
    ----------
    .. [1] Fereshtehnejad, S. M., Zeighami, Y., Dagher, A., & Postuma, R. B.
       (2017). Clinical criteria for subtyping Parkinson’s disease: biomarkers
       and longitudinal progression. Brain, 140(7), 1959-1976.
    """

    # generate clusters from "raw" behavioral data
    raw_beh = hdf.load('/raw/pd_behavioral_measures')
    feresh_labels = cluster_fereshtehnejad2017(raw_beh)

    # examine group distribution (incl. percentages from 2017 clusters)
    groups = ['Mild motor-predominant', 'Intermediate', 'Diffuse malignant']
    uniq_labs, feresh_counts = np.unique(feresh_labels, return_counts=True)
    feresh_orig = [223, 146, 52]  # original group distributions for comparison
    for grp, count, orig in zip(groups, feresh_counts, feresh_orig):
        print(f'{grp + ":" :<23} {count} '
              f'({100 * count / sum(feresh_counts):.2f}% vs'
              f' {100 * orig / sum(feresh_orig):.2f}%)')

    # NMI of feresh_labels with original consensus clustering solution
    consensus = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/consensus')
    feresh_nmi = metrics.v_measure_score(feresh_labels, consensus)
    print(f'\nNMI with consensus clustering solution: {feresh_nmi:.3f}')

    # load in data for running ANOVAs
    keys = [
        'cortical_thickness',
        'subcortical_volume',
        'dat_scans',
        'csf_assays',
        'behavioral_measures'
    ]
    data = [hdf.load(f'/processed/pd_{key}') for key in keys]
    demographics = pd.DataFrame(dict(participant=raw_beh.index,
                                     cluster=feresh_labels))

    # check whether groups differ with respect to PD atrophy
    run_pdatrophy_anova(demographics, verbose=True, run_tukey=False)

    # run a mass-univariate anova on the Fereshtehnejad 2017 clusters
    stats = run_univariate_anova(data, demographics,
                                 verbose=False, run_tukey=False)
    print(f'{len(stats)} features differentiable across groups (q < 0.05)')

    return feresh_labels


def compare_alt_distance(hdf):
    """
    Compares SNF results for alternative distance metrics

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file containing SNF gridsearch outputs

    Returns
    -------
    clusterings : list of numpy.ndarray
        Cluster labels for different distance metrics (squeclidean, cityblock,
        cosine)
    embeddings : list of numpy.ndarray
        Embeddings for different distance metrics (sqeuclidean, cityblock,
        cosine)
    """

    # load SNF gridsearch results for different distance metrics
    euclidean = hdf['/snf/processed/all/sqeuclidean/gridsearch']
    cityblock = hdf['/snf/processed/all/cityblock/gridsearch']
    cosinesim = hdf['/snf/processed/all/cosine/gridsearch']

    # grab consensus clustering solutions from results and compute NMI
    clusterings = [
        res['consensus'] for res in (euclidean, cityblock, cosinesim)
    ]
    nmi = snf.metrics.nmi(clusterings)[np.triu_indices(3, k=1)]
    print(f'NMI: {nmi}')

    # we want to compare correlations between aligned embeddings, so align
    # the embeddings from the other distance metrics to the average sqeuclidean
    # embedding
    embeddings = [euclidean['embedding']]
    for n, res in enumerate([cityblock, cosinesim]):
        fusion = res['fusion'][get_zrand_mask(res['zrand'])]
        emb = [utils.dme(network, n_components=10) for network in fusion]
        aligned = [e @ orthogonal_procrustes(e, embeddings[0])[0] for e in emb]
        embeddings.append(np.mean(aligned, axis=0))

    # correlate embeddings for different distance metrics
    metrics = itertools.combinations(['sqeuclidean', 'cityblock', 'cosine'], 2)
    for (embed1, embed2) in itertools.combinations(embeddings, 2):
        embed_corrs = utils.efficient_corr(embed1, embed2)
        print(r'Embedding corrs {}: {:.2f} $\pm$ {:.2f} [{:.2f}--{:.2f}]\n'
              .format(next(metrics),
                      np.mean(embed_corrs),
                      np.std(embed_corrs, ddof=1),
                      np.min(embed_corrs),
                      np.max(embed_corrs)))

    return clusterings, embeddings


def compare_hyperparameter_variation(hdf):
    """
    Compares SNF results from selected regions of hyperparameter space

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file containing SNF gridsearch outputs

    Returns
    -------
    cluster_nmi : numpy.ndarray
        NMI between clustering solutions from different regions of
        hyperparameter space
    embed_corrs : numpy.ndarray
        Correlations between embedding dimensions from different regions of
        hyperparameter space
    """

    fusion = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/fusion')
    labels = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/labels')
    zrand = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/zrand')
    zrand_mask = get_zrand_mask(zrand)
    fusion, labels = fusion[zrand_mask], labels[zrand_mask]
    _, realigned = get_embedding_variance(fusion, verbose=False)

    K = np.arange(5, 105)
    mu = np.logspace(np.log10(0.3), np.log10(10), 100)
    K, mu = np.meshgrid(K, mu)
    K[~zrand_mask] = 0
    mu[~zrand_mask] = 0

    # get centers of a few "stable" clusters
    mulabels = utils.cluster_img_2d(mu, threshold=0)[1]
    cm = np.around(center_of_mass(mu, mulabels, [3, 7, 9, 10])).T.astype(int)
    hv_idxs = []
    for f in tuple(zip(*cm)):
        where = np.all(np.column_stack(np.where(zrand_mask)) == f, axis=1)
        hv_idxs.append(np.where(where)[0][0])

    # now get NMI + correlations for all pairs of selected hyperparameters
    cluster_nmi, embed_corrs = [], []
    for e1, e2 in itertools.combinations(hv_idxs, 2):
        nmi = [snf.metrics.nmi(labels[[e1, e2], f])[0, 1] for f in range(3)]
        corrs = utils.efficient_corr(realigned[e1], realigned[e2])[:5]
        cluster_nmi.append(nmi)
        embed_corrs.append(corrs)

    # print summary stats for this nonsense
    cluster_nmi, embed_corrs = np.array(cluster_nmi), np.array(embed_corrs)
    print(r'Cluster two   = {:.2f} $\pm$ {:.2f} [{:.2f}--{:.2f}]'
          .format(cluster_nmi[:, 0].mean(), cluster_nmi[:, 0].std(ddof=1),
                  cluster_nmi[:, 0].min(), cluster_nmi[:, 0].max()))
    print(r'Cluster three = {:.2f} $\pm$ {:.2f} [{:.2f}--{:.2f}]'
          .format(cluster_nmi[:, 1].mean(), cluster_nmi[:, 1].std(ddof=1),
                  cluster_nmi[:, 1].min(), cluster_nmi[:, 1].max()))
    print(r'Cluster four  = {:.2f} $\pm$ {:.2f} [{:.2f}--{:.2f}]'
          .format(cluster_nmi[:, 2].mean(), cluster_nmi[:, 2].std(ddof=1),
                  cluster_nmi[:, 2].min(), cluster_nmi[:, 2].max()))
    print(r'Embeddings    = {:.2f} $\pm$ {:.2f} [{:.2f}--{:.2f}]'
          .format(embed_corrs.mean(), embed_corrs.std(ddof=1),
                  embed_corrs.min(), embed_corrs.max()))

    return cluster_nmi, embed_corrs


def compare_hcpd_snf(hdf):
    """
    Assesses SNF clustering / embedding results when both HC + PD patients used

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file containing SNF gridsearch outputs

    Returns
    -------
    fig : matplotlib.figures.Figure
        Plotted figure
    """

    # load results from the combined PD/HC SNF gridsearch
    res = hdf['/snf/processed/pdhc/sqeuclidean/gridsearch']
    dimension1, dimension2 = res['embedding'].T[:2]
    consensus = res['consensus']

    # ground truth labels
    true = np.hstack((np.zeros(len(hdf['/raw/pd_demographics'])),
                      np.ones(len(hdf['/raw/hc_demographics']))))

    # plot scatterplot of ground truth vs HC/PD clusterings olution
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    for ax, clusters, title in zip(axes, [true, consensus],
                                   ['Ground truth', 'SNF clustering']):
        ax.scatter(dimension1, dimension2, c=clusters, s=60, linewidth=0.5,
                   cmap=ListedColormap(defaults.three_cluster_cmap),
                   edgecolor=defaults.edgegray)
        ax.set(title=title, xlabel='dimension 1', ylabel='dimension 2',
               yticklabels=[], xticklabels=[])
        sns.despine(ax=ax)
    ax.set_ylabel('')

    # calculate F1 scores
    # for F1 the actual number representing the cluster assignment matters so
    # we need to re-number the `consensus` cluster labels to match the ground
    # truth as best as possible
    consensus = cluster.match_cluster_labels(consensus, true)
    f1_scores = metrics.f1_score(true, consensus % 2, average=None)
    print('F1 scores: PD = {:.2f}, HC = {:.2f}'.format(*f1_scores))

    return fig


def main():
    # load HDF file
    fname = op.join(directories.snf, f'scale500_deterministic.h5')
    hdf = structures.Frog(fname)

    print('=' * 80)
    print('Calculating supplementary cluster demographic information\n')
    get_cluster_demographics(hdf)

    print('\n' + '=' * 80)
    print('Comparing Fereshtehnejad et al., 2017 clustering results\n')
    compare_fereshtehnejad2017(hdf)

    print('\n' + '=' * 80)
    print('Comparing alternative distance metrics\n')
    compare_alt_distance(hdf)

    print('\n' + '=' * 80)
    print('Comparing SNF results across regions of hyperparameter space\n')
    compare_hyperparameter_variation(hdf)

    print('\n' + '=' * 80)
    print('Comparing SNF results of HC & PD clustering with ground truth\n')
    fig = compare_hcpd_snf(hdf)
    if SAVE_FIGS:
        fname = op.join(directories.figs, 'pdhc_clustering')
        utils.savefig(fname, fig)


if __name__ == "__main__":
    main()
