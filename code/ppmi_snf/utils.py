# -*- coding: utf-8 -*-
"""
Teeny tiny grab-bag functions
"""

import numpy as np
from scipy import ndimage


def savefig(fname, fig, dpi=300, no_png=False, no_svg=False, figsize=None,
            **kwargs):
    """
    Saves `fig` to `fname` as png and svg

    Parameters
    ----------
    fname : str
        Path to desired output file
    fig : matplotlib.figure
        Figure to be saved
    dpi : int, optional
        Resolution at which to save figure. Default: 300
    no_png : boo, optional
        Whether to not save PNG. Default: False
    no_svg : bool, optional
        Whether to not save SVG. Default: False
    figsize : tuple, optional
        Length-2 tuple specifying desired figure (width, height) in inches;
        will resize `fig` before saving. If not specified, default figure size
        will be used. Default: None
    """

    if figsize is not None:
        fig.set_size_inches(figsize)

    fname = fname.rsplit('.', 1)[0]
    save_kwargs = {**dict(bbox_inches='tight'), **kwargs}
    if not no_png:
        fig.savefig(f'{fname}.png', dpi=dpi, **save_kwargs)
    if not no_svg:
        fig.savefig(f'{fname}.svg', dpi=dpi, **save_kwargs)


def colorbar(palette, figsize=(10, 1)):
    """
    Plot the values in a color palette as a horizontal array.

    Lightly modified from `seaborn` implementation so that you can plot dense
    colorbars

    Parameters
    ----------
    pal : sequence of matplotlib colors
        Colors, i.e. as returned by seaborn.color_palette()
    """

    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    import seaborn as sns

    n = len(palette)
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.imshow(np.arange(n).reshape(1, n), cmap=ListedColormap(list(palette)),
              interpolation='nearest', aspect='auto')
    ax.set(xticks=np.arange(n) - .5, yticks=[-.5, .5],
           xticklabels=[], yticklabels=[])
    sns.despine(ax=ax, left=True, bottom=True)

    return fig


def shift_axis(ax, lshift=None, rshift=None, tshift=None, bshift=None):
    """
    Shifts position of `ax` by provided amount

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis to be shifted
    {l,r,t,b}shift : float, optional
        Specifies amount by which the {left,right,top,bottom} of `ax` should be
        moved. If not specified then specified side of `ax` remains the same.
        Default: None

    Returns
    -------
    ax : matplotlib.axes.Axes
        Shifted axis
    """

    from matplotlib.transforms import Bbox

    pos = ax.get_position().get_points()

    # no shift specified
    if all(f is None for f in [lshift, rshift, tshift, bshift]):
        return ax

    # always keep left/right shift equal if not otherwise specified
    if lshift is not None and rshift is None:
        rshift = lshift
    elif rshift is not None and lshift is None:
        lshift = rshift
    # no left/right shift specified
    elif lshift is None and rshift is None:
        lshift = rshift = 0

    # always keep top/bottom shift equal if not otherwise specified
    if tshift is not None and bshift is None:
        bshift = tshift
    elif bshift is not None and tshift is None:
        tshift = bshift
    # no top/bottom shift specified
    elif tshift is None and bshift is None:
        tshift = bshift = 0

    shift = np.array([[lshift, bshift], [rshift, tshift]])
    ax.set_position(Bbox(pos + shift))

    return ax


def rainplot(x, y, data, viol_kws=None, strip_kws=None, ax=None, palette=None):
    """
    Quick and dirty raincloud plot, a la [1]_

    Parameters
    ----------
    {x,y} : str
        Columns in `data` to be plotted on the {x,y}-axis
    data : pandas.DataFrame
    viol_kws : dict, optional
    strip_kws : dict, optional
    ax : matplotlib.axes.Axes, optional
    palette : palette name, list, or dict, optional
        Colors to use for the different levels of the hue variable. Should be
        something that can be interpreted by color_palette(), or a dictionary
        mapping hue levels to matplotlib colors.

    Returns
    -------
    ax : matplotlib.axes.Axes

    References
    ----------
    .. [1] Allen, M., Poggiali, D., Whitaker, K., Marshall, T. R., & Kievit, R.
       (2018). Raincloud plots: a multi-platform tool for robust data
       visualization. PeerJ Preprints, 6, e27137v1.
    """

    import matplotlib.pyplot as plt
    import seaborn as sns

    if ax is None:
        ax = plt.gca()

    # make violin plot
    viol_opts = dict(saturation=1, linewidth=1.25, cut=1, palette=palette,
                     edgecolor=(0.24313725, 0.24313725, 0.24313725, 1.0))
    if viol_kws is not None:
        if 'inner' in viol_kws and viol_kws.get('inner') != 'quartile':
            raise ValueError('Provided value for \'inner\' in `viol_kws` is '
                             'invalid. Must be one of: [\'quartile\']')
        viol_opts.update(**viol_kws)
    ax = sns.violinplot(x=x, y=y, data=data, inner='quartile',
                        ax=ax, **viol_opts)

    # remove bottom half of violin
    num_viols = len(ax.collections)
    for b in ax.collections:
        path = b.get_paths()[0]
        vert = path.vertices[:, 1]
        path.vertices[:, 1] = np.clip(vert, -np.inf, vert.mean())

    # remove bottom half of the quartile lines thing
    for f in ax.lines:
        xpos, ypos = f.get_data()
        ypos[1] = ypos.mean()
        f.set_data(xpos, ypos)

    # make strip plot
    strip_opts = dict(dodge=True, palette=palette, linewidth=0)
    if strip_kws is not None:
        strip_opts.update(**strip_kws)
    ax = sns.stripplot(x=x, y=y, data=data, ax=ax, **strip_opts)
    # offset strips so they don't overlap with half violins (move down by ~0.5)
    for b in ax.collections[num_viols:]:
        offset = b.get_offsets()
        offset[:, 1] += 0.25
        b.set_offsets(offset)

    return ax


def dme(network, threshold=90, n_components=10, return_result=False, **kwargs):
    """
    Threshold, cosine similarity, and diffusion map embed `network`

    Parameters
    ----------
    network : (N, N) array_like
        Symmetric network on which to perform diffusion map embedding
    threshold : [0, 100] float, optional
        Threshold used to "sparsify" `network` prior to embedding. Default: 90
    n_components : int, optional
        Number of components to retain from embedding of `network`. Default: 10
    return_result : bool, optional
        Whether to return result dictionary including eigenvalues, original
        eigenvectors, etc. from embedding. Default: False
    kwargs : key-value pairs, optional
        Passed directly to :func:`mapalign.embed.compute_diffusion_map`

    Returns
    -------
    embedding : (N, C) numpy.ndarray
        Embedding of `N` samples in `C`-dimensional spaces
    res : dict
        Only if `return_result=True`
    """

    from mapalign import embed
    from sklearn import metrics
    from sklearn.utils.extmath import _deterministic_vector_sign_flip

    # threshold
    network = network.copy()
    threshold = np.percentile(network, threshold, axis=1, keepdims=True)
    network[network < threshold] = 0

    # cosine similarity
    network = metrics.pairwise.cosine_similarity(network)

    # embed (and ensure consistent output with regard to sign flipping)
    emb, res = embed.compute_diffusion_map(network, n_components=n_components,
                                           return_result=True, **kwargs)
    emb = _deterministic_vector_sign_flip(emb.T).T

    if return_result:
        return emb, res

    return emb


def efficient_corr(x, y):
    """
    Computes correlation of matching columns in `x` and `y`

    Parameters
    ----------
    x, y : (N, M) array_like
        Input data arrays

    Returns
    -------
    corr : (M,) numpy.ndarray
        Correlations of columns in `x` and `y`
    """

    from scipy.stats import zscore

    # we need 2D arrays
    x, y = np.vstack(x), np.vstack(y)

    # check shapes
    if x.shape != y.shape:
        if x.shape[-1] != 1 and y.shape[-1] != 1:
            raise ValueError('Provided inputs x and y must either have '
                             'matching shapes or one must be a column '
                             'vector.\nProvided data:\n\tx: {}\n\ty: {}'
                             .format(x.shape, y.shape))

    corr = np.sum(zscore(x, ddof=1) * zscore(y, ddof=1), axis=0) / (len(x) - 1)

    # fix rounding errors
    corr = np.clip(corr, -1, 1)

    return corr


def cluster_img_2d(data, threshold=None, cluster=20):
    """
    Thresholds and clusters `data`

    Parameters
    ----------
    data : (N, M) array_like
        Must be 2D!
    threshold : float, optional
        All values in `data` below this are set to zero. If not specified will
        use the 95th percentile of the values in `data`. Default: None
    cluster : int, optional
        Number of adjacent values in `data` that must be above `threshold` to
        be retained

    Returns
    -------
    clusterized : (N, M) numpy.ndarray
        Provided `data` with values less than `threshold` set to 0 and clusters
        smaller than `cluster` removed
    labels : (N, M) numpy.ndarray
        Clusters detected in `clusterized`
    """

    # threshold at 95%ile by default
    data = data.copy()
    if threshold is None:
        threshold = np.percentile(data, 95)
    data[np.abs(data) < threshold] = 0

    # label image and remove groups < cluster size
    labs, nlabels = ndimage.label(data)
    for n in range(1, nlabels + 1):
        if labs[labs == n].size < cluster:
            labs[labs == n] = 0

    # mask data to remove small clusters and then re-label groups
    data = data * (labs > 0)
    labels, nlabels = ndimage.label(data)

    # return masked data and labels
    return data, labels
