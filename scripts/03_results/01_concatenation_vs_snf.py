# -*- coding: utf-8 -*-
"""
This script generates results and figures used to compare SNF with data
concatenation
"""

import os.path as op
import warnings

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sstats
import seaborn as sns
from sklearn.cluster import spectral_clustering
from sklearn import metrics

from netneurotools import modularity, plotting
from ppmi_snf import structures, directories, utils
import snf.metrics

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


def gen_mod(affinities, labels):
    """
    Generates modularity estimates for `affinities` using `labels`

    Parameters
    ----------
    affinities : list of (N, N) array_like
        Affinity matrices
    labels : (N,) array_like
        Cluster labels for matrices in `affinities`

    Yields
    ------
    modularity : float
    """

    for aff in affinities:
        yield modularity.get_modularity(aff, labels).sum()


def _print_summary(data, metric):
    """
    Prints summary statistics for `data`

    Parameters
    ----------
    data : array_like
        Input data array
    metric : str
        Metric represented by `data`
    """

    print(u'Cortical thickness {}: {:.2f} \u00B1 {:.2f} [{:.2f}--{:.2f}]'
          .format(metric, data[:, 0].mean(), data[:, 0].std(ddof=1),
                  data[:, 0].min(), data[:, 0].max()))
    print('Other modalities {}:   {:.2f} \u00B1 {:.2f} [{:.2f}--{:.2f}]'
          .format(metric, data[:, 1:].mean(), data[:, 1:].std(ddof=1),
                  data[:, 1:].min(), data[:, 1:].max()))
    print('Overall {}:            {:.2f} \u00B1 {:.2f} [{:.2f}--{:.2f}]'
          .format(metric, data.mean(), data.std(ddof=1),
                  data.min(), data.max()))


def get_nmi_mod(method, print_summary=True):
    """
    Gets normalized mutual information and modularity for `method`

    Parameters
    ----------
    method : {'snf', 'rbf'}
        Method to use for calculating metrics
    print_summary : bool, optional
        Whether to print summary statistics (mean, SD, ranges) for generated
        metrics

    Returns
    -------
    nmi : numpy.ndarray
        Normalized mutual information
    mod : numpy.ndarray
        Modularity estimates
    """

    methods = ['snf', 'rbf']
    if method not in methods:
        raise ValueError(f'Provided `method` {method} invalid.')

    scales = [f'scale{f}' for f in ['033', '060', '125', '250', '500']]
    keys = [
        'cortical_thickness',
        'subcortical_volume',
        'dat_scans',
        'csf_assays',
        'behavioral_measures',
        'all'
    ]

    # iterate over all CT dimensionalities and generate NMI / mod estimates
    nmi, mod = [], []
    for scale in scales:
        # get data for provided scale
        fname = op.join(directories.snf, f'{scale}_deterministic.h5')
        hdf = structures.Frog(fname)
        pd_data = [hdf.load(f'/processed/pd_{key}') for key in keys[:-1]]

        # generate affinity matrix and cluster labels
        # if we're using SNF we can just pre-load the matrices + labels
        if method == 'snf':
            path = '/snf/processed/{}/sqeuclidean/gridsearch/{}'
            affinities = [
                hdf.load(path.format(key, 'fusion_avg')) for key in keys
            ]
            labels = [
                hdf.load(path.format(key, 'consensus')) for key in keys
            ]
        # otherwise, we have to generate the affinities using cosine similarity
        # and then use spectral clustering to generate the labels
        elif method == 'rbf':
            affinities = [
                metrics.pairwise.cosine_similarity(
                    sstats.zscore(f)
                ) + 1 for f in pd_data
            ] + [
                metrics.pairwise.cosine_similarity(
                    sstats.zscore(np.column_stack(pd_data))
                ) + 1
            ]
            labels = [
                spectral_clustering(aff, n_clusters=3, random_state=1234)
                for aff in affinities
            ]

        # get NMI + modularity estimates
        nmi.append(snf.metrics.nmi(labels)[-1, :-1])
        mod.append(list(gen_mod(affinities[:-1], labels[-1])))

    nmi, mod = np.asarray(nmi), np.asarray(mod)

    if print_summary:
        _print_summary(nmi, 'NMI')
        print()
        _print_summary(mod, 'modularity')
        print()

    return nmi, mod


def gen_figure(snf_nmi, rbf_nmi, snf_mod, rbf_mod):
    """
    Generates figure comparing SNF/data concat. NMI + modularity estimates

    Parameters
    ----------
    snf_nmi : array_like
    rbf_nmi : array_like
    snf_mod : array_like
    rbf_mod : array_like

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    xticklabels = ['CT', 'SV', 'DB', 'CSF', 'CLIN']
    yticklabels = ['68', '114', '219', '448', '1000']

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=-0.05, hspace=0.15)

    # make four circleplots
    ax1 = plotting.circleplot(snf_nmi, vmin=-0.01, vmax=1.01, ax=axes[0][0],
                              xticklabels=[], yticklabels=yticklabels)
    ax2 = plotting.circleplot(rbf_nmi, vmin=-0.01, vmax=1.01, ax=axes[0][1],
                              xticklabels=[], yticklabels=yticklabels,
                              cbar_kws={'ticks': [0.00, 1.00]})
    ax3 = plotting.circleplot(snf_mod, vmin=-0.01, vmax=0.36, ax=axes[1][0],
                              xticklabels=xticklabels, yticklabels=yticklabels)
    ax4 = plotting.circleplot(rbf_mod, vmin=-0.01, vmax=0.36, ax=axes[1][1],
                              xticklabels=xticklabels, yticklabels=yticklabels,
                              cbar_kws={'ticks': [0, 0.35]})

    for ax in axes.flatten():
        for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)

    # set axis titles
    ax1.set_title('Similarity network fusion', pad=20)
    ax2.set_title('Data concatenation', pad=20)

    # set axis labels
    ax1.set_ylabel('Dimensionality of\ncortical thickness data',
                   labelpad=15, x=0, y=0)
    ax3.set_xlabel('Data type', labelpad=15, x=1.1)

    # turn off colorbars on lefthand plots
    ax1.collections[0].colorbar.ax.set_visible(False)
    ax3.collections[0].colorbar.ax.set_visible(False)

    # correct colorbar appearance for righthand plots
    ax2.collections[0].colorbar.ax.tick_params(size=0, labelsize=14)
    ax4.collections[0].colorbar.ax.tick_params(size=0, labelsize=14)
    ax2.collections[0].colorbar.ax.set_ylabel('Normalized mutual\ninformation',
                                              rotation=270, labelpad=30)
    ax4.collections[0].colorbar.ax.set_ylabel('Modularity',
                                              rotation=270, labelpad=15)

    # plot small gray lines to better differentiate plots
    plt.plot([0.4725, 0.4725], [0.55, 0.85], color='gray', lw=0.5,
             transform=fig.transFigure, clip_on=False)
    plt.plot([0.4725, 0.4725], [0.15, 0.45], color='gray', lw=0.5,
             transform=fig.transFigure, clip_on=False)
    plt.plot([0.155, 0.415], [0.5, 0.5], color='gray', lw=0.5,
             transform=fig.transFigure, clip_on=False)
    plt.plot([0.525, 0.795], [0.5, 0.5], color='gray', lw=0.5,
             transform=fig.transFigure, clip_on=False)

    return fig


def main():
    # get NMI and modularity for both SNF + data concatenation methods
    print('=' * 80)
    print('Calculating NMI + modularity of SNF outputs:\n')
    snf_nmi, snf_mod = get_nmi_mod('snf')

    print('=' * 80)
    print('Calculating NMI + modularity of data concatenation outputs\n')
    rbf_nmi, rbf_mod = get_nmi_mod('rbf')

    # use these results to generate what will serve as the basis for figure 2
    fig = gen_figure(snf_nmi, rbf_nmi, snf_mod, rbf_mod)
    if SAVE_FIGS:
        fname = op.join(directories.figs, 'data_concatenation')
        utils.savefig(fname, fig)


if __name__ == "__main__":
    main()
