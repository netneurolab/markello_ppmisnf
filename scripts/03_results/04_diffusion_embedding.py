# -*- coding: utf-8 -*-
"""
This script generates results and figures used to investigate a diffusion map
embedding of the PD patient similarity network
"""

import os.path as op
import warnings

from matplotlib import gridspec
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from matplotlib.text import Text
from nilearn.datasets import fetch_neurovault_ids
from nilearn.plotting import plot_glass_brain
import numpy as np
import seaborn as sns
import scipy.stats as sstats
from sklearn.model_selection import KFold
from statsmodels.stats.multitest import fdrcorrection

from ppmi_snf import defaults, directories, structures, utils
from netneurotools.utils import add_constant

from analysis import (get_embedding_variance, get_zrand_mask,
                      run_pdatrophy_anova)

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

# plot settings
scatter_kws = {
    'color': defaults.gray,
    'edgecolor': defaults.edgegray,
    's': 50,
    'linewidths': 0.5,
    'alpha': 1.0
}
line_kws = {
    'color': 'black',
    'linewidth': 1.5
}
axis_kws = {
    'xticklabels': (),
    'yticks': (-3, 0, 3),
    'ylim': (-3.5, 3.5)
}


def run_prediction_models(hdf, feats=None, verbose=True):
    """
    Runs model using diffusion embedding scores to predict behavioral measures

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file containing SNF gridsearch outputs
    feats : list of str, optional
        List of behavioral features to use as prediction targets

    Returns
    -------
    Y_corrs : (K, F) numpy.ndarray
        Correlation between predicted and actual behavioral values for `F`
        features across `K` folds
    Y_mses : (N, F) numpy.ndarray
        Mean-squared error of predicted and actual behavioral values for `F`
        features across `N` subjects
    """

    if feats is None:
        feats = ['pigd', 'tremor']

    holdout = hdf.load('/snf/processed/holdout/all/sqeuclidean/gridsearch')
    X_holdout = holdout['embedding'][:, :5]
    behavior = hdf.load('/processed/pd_behavioral_measures')
    Y_holdout = np.asarray(behavior[feats])
    consensus = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/consensus')

    # to store out-of-sample correlations and MSE scores
    n_splits = 5
    Y_corrs = np.zeros((n_splits, 2))
    Y_mses = np.zeros_like(Y_holdout)

    # 5-fold CV
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=SEED)
    for n, (train_index, test_index) in enumerate(kf.split(X_holdout)):
        # split X and Y into train/test
        X_train, X_test = X_holdout[train_index], X_holdout[test_index]
        Y_train, Y_test = Y_holdout[train_index], Y_holdout[test_index]

        # zscore / zmap and add constant to X matrix
        X_test = add_constant(sstats.zmap(X_test, X_train, ddof=1))
        X_train = add_constant(sstats.zscore(X_train, ddof=1))
        Y_test = sstats.zmap(Y_test, Y_train, ddof=1)
        Y_train = sstats.zscore(Y_train, ddof=1)

        # fit model and predict out-of-sample
        betas = np.linalg.lstsq(X_train, Y_train, rcond=None)[0]
        Y_pred = X_test @ betas

        # get correlation and MSE
        Y_corrs[n] = utils.efficient_corr(Y_pred, Y_test)
        Y_mses[test_index] = (Y_test - Y_pred) ** 2
        Y_mse_mean = np.mean(Y_mses[test_index], axis=0)

        if verbose:
            print(f'Fold {n + 1}: r = {Y_corrs[n]:}, mse = {Y_mse_mean:}')

    if verbose:
        print('\nAverage correlations across folds:')
        corrs_mean, corrs_std = Y_corrs.mean(0), Y_corrs.std(0, ddof=1)
        for n, t in enumerate(feats):
            print(r'{:<9}: r = {:.3f} $\pm$ {:.3f}'
                  .format(t, corrs_mean[n], corrs_std[n]))

        print('\nGroups differences in MSE:')
        f_hold, p_hold = sstats.f_oneway(
            *(Y_mses[consensus == cl] for cl in np.unique(consensus))
        )
        for n, t in enumerate(feats):
            print('{:<9}: F = {:.2f}, p = {:.3f}'
                  .format(t, f_hold[n], p_hold[n]))

    return Y_corrs, Y_mses


def _clean_ylabel(feature, max_length=20):
    """
    Replaces underscores with spaces and splits `feature` based on line length

    Parameters
    ----------
    feature : str
        String to be cleaned
    max_length : str
        Maximum length (in characters) of each line. If `feature` is longer
        than this length it will be split with a newline character. Default: 20

    Returns
    -------
    feature : str
        Cleaned input `feature`
    """

    feature = feature.replace('_', ' ')
    if len(feature) > max_length:
        ylabel = feature.split(' ')
        idx = len(ylabel) // 2
        feature = '\n'.join([' '.join(ylabel[:idx]), ' '.join(ylabel[idx:])])
    return feature


def gen_scatterplots(data, embedding):
    """
    Generates scatterplots of features with embedding dimensions

    Parameters
    ----------
    data : list of pandas.DataFrame
        Data provided as input to SNF
    embedding (N, C) array_like
        Average scores for `N` patients across `C` embedding dimensions

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plotted figure
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

    fig, (top, bot) = plt.subplots(2, 5, figsize=(25, 10))
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    corr_idxs = np.zeros(5, dtype=int)

    for n, dimension in enumerate(embedding.T[:5]):
        # correlate all features with the selected embedded dimension and
        # extract index of strongest positive + negative correlation
        mc = utils.efficient_corr(dimension, data)
        t_idx, b_idx = mc.argmax(), mc.argmin()

        # store the index of the ABSOLUTE strongest correlation for this dim
        corr_idxs[n] = t_idx if np.abs(mc.max()) > np.abs(mc.min()) else b_idx

        # make regression plot of feature by embedded dimension
        for ax, idx in zip([top[n], bot[n]], [t_idx, b_idx]):
            sns.regplot(dimension, data[:, idx], ax=ax, ci=None,
                        scatter_kws=scatter_kws, line_kws=line_kws)
            ax.set_title('{}\nr = {:>5.2f}'
                         .format(datatypes[idx], mc[idx]), pad=10)
            ax.set(ylabel=_clean_ylabel(features[idx]),
                   xlabel='', **axis_kws)
            sns.despine(ax=ax)

        # only set the x-label for the bottom plots
        ax.set(xlabel=f'dimension {n + 1}')

    for (ax, text) in zip([top[0], bot[0]], ['a', 'b']):
        text = ax.text(-0.2, 1.15, text, transform=ax.transAxes,
                       fontdict=fontd)

    return fig, corr_idxs


def gen_figure(data, embedding, realigned, consensus, pdatrophy, corr_idxs):
    """
    Generates figure

    Parameters
    ----------
    data : list of pandas.DataFrame
        Data provided as input to SNF
    embedding : (N, C) array_like
        Average scores for `N` patients across `C` embedding dimensions
    realigned : (A, N, C) array_like
        Aligned scores of `N` patients across `C` embedding dimensions for `A`
        different embeddings
    consensus : (N,) array_like
        Consensus clustering assignments for `N` patients
    pdatrophy : (N,) array_like
        PD-ICA atrophy scores for `N` patients
    corr_idxs : (5,) array_like
        Indices of data features with maximum correlation to first five
        dimensions of `embedding`

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plotted figure
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

    # gridspec so that brain map takes up 2x space of scatterplot
    fig = plt.figure(figsize=(26, 9))
    gs0 = gridspec.GridSpec(2, 1, figure=fig)
    gs0.update(wspace=0.2, hspace=0.4)
    gs = gs0[0].subgridspec(1, 4)

    ###############################################################################
    # show how clusters are located in embedded space

    ax1, ax2 = plt.subplot(gs[0]), plt.subplot(gs[1])
    for sub, cmap in zip([60, 6, 2], defaults.three_cluster_cmap):
        d = np.column_stack([f[sub, :2] for f in realigned])
        ax1.scatter(*d, c=[cmap], alpha=0.05, rasterized=True)
        ax1.scatter(*d.mean(axis=1), c='black', edgecolor='white')
    ax2.scatter(*embedding.T[:2], c=consensus, rasterized=True,
                cmap=ListedColormap(defaults.three_cluster_cmap),
                edgecolor=defaults.edgegray, s=60, linewidth=0.5)
    for ax in [ax1, ax2]:
        ax.set(xlim=(-35, 37), ylim=(-27, 32),
               xlabel='dimension 1', ylabel='',
               xticklabels=[], yticklabels=[])
        sns.despine(ax=ax)
    ax1.set_ylabel('dimension 2')

    ax1.set_title('subject trajectories\nacross embeddings', pad=15)
    ax2.set_title('average subject\nembeddings', pad=15)
    # ax1.text(-0.2, 1.15, 'b', transform=ax1.transAxes)
    # ax2.text(-0.2, 1.15, 'c', transform=ax2.transAxes)
    # utils.shift_axis(ax1, lshift=-0.015)
    # utils.shift_axis(ax2, lshift=0.015)

    # add colorbar to figure (steal space from axes)
    cbar = fig.colorbar(ax2.collections[0], ax=[ax1, ax2],
                        drawedges=False, ticks=[],
                        boundaries=np.arange(0.5, 4.5))
    cbar.outline.set(linewidth=0)
    cbar.ax.tick_params(axis='both', which='both', length=0)
    cbar.ax.set_ylabel('cluster label', rotation=270, labelpad=25)

    ###############################################################################
    # make embedding scatterplots for first two dimensions

    ax2, ax3 = plt.subplot(gs[2]), plt.subplot(gs[3])
    vmin, vmax = -2.5, 2.5
    plotdata = [data[:, corr_idxs[0]], data[:, corr_idxs[1]]]
    titles = [datatypes[corr_idxs[0]], datatypes[corr_idxs[1]]]
    subtitles = [features[corr_idxs[0]], features[corr_idxs[1]]]
    for n, (ax, ydata) in enumerate(zip([ax2, ax3], plotdata)):
        # use mean of all features for given datatype; this is a bad
        # heuristic but we're just trying to get a feel for things
        ydata = sstats.zscore(ydata)

        # make scatterplot
        ax.scatter(*embedding.T[:2], c=ydata, vmin=vmin, vmax=vmax,
                   edgecolor=defaults.edgegray, s=60, linewidth=0.5)
        ax.set(xlabel='dimension 1',
               ylabel='dimension 2' if n == 0 else '',
               xlim=[-35, 37], ylim=[-27, 32],
               xticklabels=[], yticklabels=[])
        sns.despine(ax=ax)
        ax.set_title('{}\n{}'.format(titles[n], subtitles[n]), pad=15)

    # add colorbar
    cbar = fig.colorbar(ax.collections[0], ax=[ax2, ax3],
                        drawedges=False, ticks=[vmin, vmax])
    cbar.outline.set(linewidth=0)
    cbar.ax.tick_params(axis='both', which='both', length=0)

    ###############################################################################
    # make PD ICA brain plot

    gs = gs0[1].subgridspec(1, 6)
    ax4 = plt.subplot(gs[:3])
    pdica = fetch_neurovault_ids(image_ids=[12551],
                                 data_dir=directories.rois,
                                 verbose=0)
    plot_glass_brain(pdica['images'][0], threshold=3.0, axes=ax4)
    ax4.set_title('atrophy network', pad=15)
    # update the size of the L/R annotations on the brains
    for ax in gs0.figure.get_children()[-3:]:
        [f.set_fontsize(14) for f in ax.get_children() if isinstance(f, Text)]

    ###############################################################################
    # make PD ICA scatter plot with second embedding dimension
    ax5 = plt.subplot(gs[5])
    sns.regplot(embedding[:, 1], pdatrophy, ax=ax5, ci=None,
                scatter_kws=scatter_kws, line_kws=line_kws)
    # correlation of PD-ICA atrophy score with embedded dimensions
    rvals, pvals = zip(*[sstats.pearsonr(e, pdatrophy) for e in embedding.T])
    pvals = fdrcorrection(pvals)[1]
    for n, p in enumerate(pvals):
        if p < 0.05:
            print(f'Dimension {n + 1}: r = {rvals[n]:>5.2f}, p = {p:.3f}')
    ax5.set(ylabel='atrophy z-score', xlabel='dimension 2', **axis_kws)
    ax5.set_title(f'r = {rvals[1]:>5.2f}', pad=10)
    utils.shift_axis(ax5, lshift=-0.2, rshift=-0.16)
    sns.despine(ax=ax5)

    return fig


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

    # also load the gridsearch results back in to memory.
    # here, fusion is shape (K, M, N, N),
    #        zrand is shape (C, K, M)
    # where `K` is the nearest-neighbors parameter of SNF
    #       `M` is the scaling (mu) parameter of SNF, and
    #       `N` is PD patients
    fusion = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/fusion')
    zrand = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/zrand')
    consensus = hdf.load('/snf/processed/all/sqeuclidean/gridsearch/consensus')

    print('=' * 80)
    print('Calculating variance explained by diffusion map embedding\n')
    mask = get_zrand_mask(zrand)
    embedding, realigned = get_embedding_variance(fusion[mask])

    print('\n' + '=' * 80)
    print('Calculating prediction model performance\n')
    run_prediction_models(hdf, feats=['pigd', 'tremor'])

    print('\n' + '=' * 80)
    print('Calculating diffusion map embedding dimension correlations\n')
    fig, corr_idxs = gen_scatterplots(data, embedding)
    if SAVE_FIGS:
        fname = op.join(directories.figs, 'diffusion_correlations')
        utils.savefig(fname, fig)

    # load demographics information
    demographics = hdf.load('/raw/pd_demographics').reset_index()
    demographics = demographics.assign(cluster=consensus)
    pdatrophy = run_pdatrophy_anova(demographics, verbose=False,
                                    run_tukey=False)['atrophy']
    fig = gen_figure(data, embedding, realigned, consensus, pdatrophy,
                     corr_idxs)
    if SAVE_FIGS:
        fname = op.join(directories.figs, 'diffusion_embedding')
        utils.savefig(fname, fig)


if __name__ == "__main__":
    main()
