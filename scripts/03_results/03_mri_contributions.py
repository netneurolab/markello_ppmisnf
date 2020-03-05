# -*- coding: utf-8 -*-
"""
This script generates results and figures used to assess the impact of
neuroimaging data on PD patient clustering solutions generated with SNF
"""

import os.path as op
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import metrics
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import fdrcorrection

from netneurotools import cluster
from ppmi_snf import directories, structures, utils

from analysis import run_pdatrophy_anova, run_univariate_anova

# set seaborn, numpy, matplotlib, warnings options
sns.set(style='white', context='notebook', font_scale=1.5)
warnings.filterwarnings('ignore', category=PendingDeprecationWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning,
                        message='Calling np.sum')
plt.rcParams['svg.fonttype'] = 'none'
fontd = dict(size=36)

# miscellaneous settings
SAVE_FIGS = True  # whether to save figures
SEED = 1234  # random seed for reproducibility
np.random.seed(SEED)  # set random seed


def compare_nomri_clusters(hdf):
    """
    Compares ALL versus NO-MRI consensus cluster labels

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file containing SNF gridsearch outputs

    Returns
    -------
    all_consensus : (N,) numpy.ndarray
        Consensus cluster labels generated from SNF with ALL (i.e., five)
        data modalities
    nomri_consensus : (N,) numpy.ndarray
        Consensus cluster labels generated from SNF with no MRI (i.e., three)
        data modalities
    """

    # load consensus cluster assignments and average embedding for ALL data
    path = '/snf/processed/{}/sqeuclidean/gridsearch/{}'
    all_consensus = hdf.load(path.format('all', 'consensus'))
    all_embedding = hdf.load(path.format('all', 'embedding'))

    # now, load consensus + average embedding for no-MRI data
    nomri_consensus = hdf.load(path.format('nomri', 'consensus'))
    nomri_embedding = hdf.load(path.format('nomri', 'embedding'))

    # match cluster labels for ALL and no-MRI consensus assignments
    nomri_consensus = cluster.match_cluster_labels(nomri_consensus,
                                                   all_consensus)
    nmi = metrics.v_measure_score(all_consensus, nomri_consensus)
    print(f'NMI between clustering solutions: {nmi:.3f}')

    # check correlation between ALL and no-MRI average embeddings
    embedding_corrs = utils.efficient_corr(all_embedding, nomri_embedding)
    print('Correlations b/w embeddings: {:.2f} +/- {:.2f} [{:.2f} - {:.2f}]'
          .format(np.mean(embedding_corrs),
                  np.std(embedding_corrs, ddof=1),
                  np.min(embedding_corrs),
                  np.max(embedding_corrs)))

    return all_consensus, nomri_consensus


def run_twoway_anova(data, all_consensus, nomri_consensus):
    """
    Runs two-way ANOVA to assess clustering solution difference in `data`

    Parameters
    ----------
    data : list of pandas.DataFrame
        Data provided as input to SNF
    all_consensus : (N,) numpy.ndarray
        Consensus cluster labels generated from SNF with ALL (i.e., five)
        data modalities
    nomri_consensus : (N,) numpy.ndarray
        Consensus cluster labels generated from SNF with no MRI (i.e., three)
        data modalities

    Returns
    -------
    features : list of str
        List of features for which the two-way ANOVA is significant (after FDR
        correction)
    """
    # generate dataframe to test differences between clustering solutions
    features = np.hstack([d.columns for d in data])
    data = pd.concat(data, axis=1)
    data = data.assign(mri=pd.Categorical(all_consensus),
                       nomri=pd.Categorical(nomri_consensus))
    features = np.array([f.replace('-', '_') for f in features])
    data.columns = np.hstack([features, 'mri', 'nomri'])

    # run two-way anova b/w clusters and clustering solutions
    omnibus, interaction = [], []
    for feat in features:
        model = smf.ols(f'{feat} ~ C(nomri) * C(mri)', data).fit()
        omnibus.append(model.f_pvalue)
        interaction.append(sm.stats.anova_lm(model)['PR(>F)'][2])

    # FDR correct both omnibus + interaction term and then determine for which
    # features, if any, BOTH terms are q < 0.05
    omnibus = fdrcorrection(omnibus)[1]
    interaction = fdrcorrection(interaction)[1]
    sig_features = features[np.logical_and(omnibus < 0.05, interaction < 0.05)]

    print('{} features significantly different between clustering solutions'
          .format(len(sig_features)))

    return sig_features


def main():
    keys = [
        'cortical_thickness',
        'subcortical_volume',
        'dat_scans',
        'csf_assays',
        'behavioral_measures'
    ]

    # load processed data
    fname = op.join(directories.snf, f'scale500_deterministic.h5')
    hdf = structures.Frog(fname)
    data = [hdf.load(f'/processed/pd_{key}') for key in keys]

    # compare clustering + embedding results for ALL vs NO-MRI data
    print('=' * 80)
    print('Comparing SNF outputs with and without MRI features\n')
    all_consensus, nomri_consensus = compare_nomri_clusters(hdf)

    # generate demographic dataframe with NO-MRI cluster labels
    demographics = hdf.load('/raw/pd_demographics').reset_index()
    demographics = demographics.assign(cluster=pd.Categorical(nomri_consensus))

    # run one-way and two-way ANOVA to assess cluster discriminability
    print('\n' + '=' * 80)
    print('Testing no-MRI cluster differences in PD-ICA atrophy score\n')
    run_pdatrophy_anova(demographics)

    print('\n' + '=' * 80)
    print('Running mass-univariate ANOVA for no-MRI cluster differences\n')
    run_univariate_anova(data, demographics, run_tukey=False)

    print('\n' + '=' * 80)
    print('Running two-way ANOVA comparing SNF w/ and w/o MRI features\n')
    run_twoway_anova(data, all_consensus, nomri_consensus)


if __name__ == "__main__":
    main()
