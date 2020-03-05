# -*- coding: utf-8 -*-
"""
This script runs an SNF gridsearch on different combinations of input data; it
assumes that `01_prepare_snf_data.py` has been run previously to generate HDF5
files that contain this input data. A gridsearch will be run on each file
separately.

Given a set of input patient x feature arrays from the data file, this script
will generate outputs for each of 10,000 unique combinations of SNF
hyperparameters. It will:

    1. Run similarity network fusion with the provided hyperparameters to
       generate a fused patient similarity network, and
    2. Cluster the fused network to generate two-, three-, and four-cluster
       patient assignments.

Once all the fused networks and cluster assignments are generated for each of
the 10,000 hyperparameters combinations, this script will assess the stability
of clustering assignments for local regions of the hyperparameter space using a
z-Rand function. All these outputs (fused networks, clustering assignments, and
clustering stability) will be saved back to disk to reduce computation time
for later analyses.
"""

import glob
import itertools
import os.path as op

from joblib import delayed, Parallel
from mapalign import align
import numpy as np
import pandas as pd
from sklearn.cluster import spectral_clustering
import tqdm

from netneurotools import cluster, stats
from ppmi_snf import directories, structures
from ppmi_snf.utils import cluster_img_2d, dme
from snf import compute
from snf.cv import zrand_convolve

N_PROC = -1  # amount of parallelization for grid-search (-1 = all CPUs)
KEYS = [
    'cortical_thickness',
    'subcortical_volume',
    'dat_scans',
    'csf_assays',
    'behavioral_measures'
]


def load_and_residualize_data(hdf, groups=None):
    """
    Loads raw data and returns residualized outputs

    Parameters
    ----------
    hdf : structures.Frog
        HDF5 file with saved data
    groups : list of str, optional
        Which groups to load. Must be in ['pd', 'hc']. If not specified all
        groups are loaded. Default: None

    Returns
    -------
    data : list of pandas.DataFrame
        Residualized data
    """

    if isinstance(hdf, str):
        hdf = structures.Frog(hdf)

    if groups is None:
        groups = ['pd', 'hc']

    raw_data = [
        pd.concat([hdf[f'/raw/{group}_{key}'] for group in groups])
        for key in KEYS
    ]
    regressors = [
        pd.concat([hdf[f'/regressors/{group}_{key}'] for group in groups])
        for key in KEYS
    ]

    proc_data = []
    for data, reg in zip(raw_data, regressors):
        resid = pd.DataFrame(stats.residualize(reg, data, normalize=False),
                             index=data.index, columns=data.columns)
        proc_data.append(resid)

    return proc_data


def fuse_and_label(data, K, mu, n_clusters, metric):
    """
    Generates fusion + cluster assignments for given hyperparameters

    Small helper function to be used for parallelization of gridsearch

    Parameters
    ----------
    data : list of numpy.ndarray
    K : int
    mu : float
    n_clusters : list of int
    metric : str

    Returns
    -------
    fusion : numpy.ndarray
    labels : list of numpy.ndarray
    """

    aff = compute.make_affinity(*data, K=K, mu=mu, metric=metric,
                                normalize=True)

    if isinstance(aff, list) and len(aff) > 1:
        fusion = compute.snf(*aff, K=K)
    else:
        fusion = aff

    labels = [
        spectral_clustering(fusion, ncl, random_state=1234)
        for ncl in n_clusters
    ]

    return fusion, labels


def run_gridsearch(data, hdf, path, metrics=None, saveall=True):
    """
    Runs gridsearch on `data` and saves outputs to `hdf`[`path`]

    Parameters
    ----------
    data : list of array_like
        Data on which to run SNF gridsearch
    hdf : str or structures.Frog
        Filepath to or loaded structures.Frog object to save output data
    path : str
        Will be inserted into "/snf/{path}/{metric}/gridsearch", specifying
        the path in `hdf` to which gridsearch results will be saved
    metrics : list of str, optional
        Which distance metrics SNF should be run with. If not specified will
        use ['sqeuclidean', 'cityblock', 'cosine']. Default: None
    saveall : bool, optional
        Whether to save all outputs of gridsearch (i.e., all fused matrices,
        all clustering assignments, z-rand convolved similarity matrices, AND
        consensus clustering assignments) instead of only consensus clusters,
        average fused matrix, and agreement matrix. Default: True

    Returns
    -------
    hdf :  structures.Frog
        Same as provided input but with new gridsearch results!
    """

    if metrics is None:
        metrics = ['sqeuclidean', 'cityblock', 'cosine']
    elif isinstance(metrics, str):
        metrics = [metrics]

    if isinstance(hdf, str):
        hdf = structures.Frog(hdf)

    n_subj = len(data[0])
    fname = op.basename(hdf.filename)
    print(f'Running grid-search for {fname} with {len(data)} datatypes; '
          f'saving to path "{path}"')

    # set K / mu (hyperparameters) that we'll explore (10,000 combinations)
    K = np.arange(5, 105)
    mu = np.logspace(np.log10(0.3), np.log10(10), 100)
    # only consider two-, three-, and four-cluster solutions in this space
    n_clusters = [2, 3, 4]

    for metric in metrics:
        # check that the gridsearch wasn't already run for this combination.
        # no need to repeat needless computations!
        mpath = f'/snf/{path}/{metric}/gridsearch'
        if mpath in hdf.groups():
            check = ['consensus', 'fusion_avg', 'agreement', 'embedding']
            if saveall:
                check += ['fusion', 'labels', 'zrand']
            if all(op.join(mpath, p) in hdf.keys() for p in check):
                continue

        # generate fused networks + cluster assignments for all the parameters
        print(f'Generating outputs from gridsearch with {metric} distance')
        fuse = delayed(fuse_and_label)
        gridres = Parallel(n_jobs=N_PROC)(
            fuse(data, k, m, n_clusters, metric)
            for k, m in tqdm.tqdm(list(itertools.product(K, mu)))
        )

        # wrangle outputs from gridsearch and reshape
        fusion, labels = [np.stack(f, axis=0) for f in zip(*gridres)]
        fusion = fusion.reshape(len(K), len(mu), n_subj, n_subj)
        labels = labels.reshape(len(K), len(mu), len(n_clusters), n_subj)

        # don't parallelize zrand_convolve across cluster solutions because
        # it's already parallelizing at a lower level
        print('Convolving cluster assignments with z-Rand kernel')
        zrand_avg = [
            zrand_convolve(labels[..., n, :], n_proc=N_PROC)
            for n in range(len(n_clusters))
        ]

        # make a record of all the gridsearch outputs if desired
        if saveall:
            results = dict(fusion=fusion, labels=labels, zrand=zrand_avg)
        else:
            results = dict()

        # we'll use the 95%ile of all the z-rand scores to threshold the
        # similarity matrices and extract "stable" regions as a mask of the
        # hyperparameter space
        zrand_thr = np.percentile(zrand_avg, 95)
        mask = [cluster_img_2d(z, zrand_thr)[1] != 0 for z in zrand_avg]
        zrand_mask = np.sum(mask, axis=0) > 0

        # only keep assignments / fused networks from stable regions
        stable_idx = np.where(zrand_mask)
        labels, fusion = labels[stable_idx], fusion[stable_idx]

        # extract stable community assignments and make consensus
        comms = labels.reshape(-1, labels.shape[-1]).T
        cons, ag = cluster.find_consensus(comms, return_agreement=True,
                                          seed=1234)
        results['consensus'] = cons

        # run diffusion map embedding and generate average, aligned embedding
        embeddings = [dme(network, n_components=10) for network in fusion]
        realigned, xfms = align.iterative_alignment(embeddings, n_iters=1)
        results['embedding'] = np.mean(realigned, axis=0)

        # we'll keep the average fused network and the agreement matrix to
        # use for calculating modularity
        results['fusion_avg'] = np.mean(fusion, axis=0)
        results['agreement'] = ag

        hdf.save(results, mpath, overwrite=True)

    return hdf


def main():
    # grab all the HDF5 files that exist
    for hdf in sorted(glob.glob(op.join(directories.snf, '*.h5'))):

        # prepare HDF file and pre-load data
        hdf = structures.Frog(hdf)
        data = [hdf.load(f'/processed/pd_{key}') for key in KEYS]

        # the only gridsearches we need to run for all the resolutions are the
        # basic one where we save out `fusion_avg` and `consensus`
        # we need ALL teh data modalities, and each data modality independently
        run_gridsearch(data=data, hdf=hdf, path='processed/all', saveall=False,
                       metrics='sqeuclidean')
        for n, key in enumerate(KEYS):
            run_gridsearch(data=[data[n]], hdf=hdf, path=f'processed/{key}',
                           saveall=False, metrics='sqeuclidean')

    # for the highest resolution data we want to run a BUNCH of auxiliary
    # analyses, though
    hdf = op.join(directories.snf, 'scale500_deterministic.h5')
    hdf = structures.Frog(hdf)
    data = [hdf.load(f'/processed/pd_{key}') for key in KEYS]

    # SNF for all non-MRI data
    run_gridsearch(data=data[2:], hdf=hdf, path='processed/nomri',
                   saveall=True, metrics='sqeuclidean')

    for n, key in enumerate(KEYS):
        # SNF for all modalities except one
        run_gridsearch(data=[d for i, d in enumerate(data) if i != n],
                       hdf=hdf, path=f'processed/no_{key}', saveall=False,
                       metrics='sqeuclidean')

    # SNF removing "holdout" behavioral variables (all + non-MRI)
    if 'behavioral_measures' in KEYS:
        idx = KEYS.index('behavioral_measures')
        data[idx] = data[idx].drop(['tremor', 'pigd'], axis=1)
    run_gridsearch(data=data, hdf=hdf, path='processed/holdout/all',
                   saveall=False, metrics='sqeuclidean')
    run_gridsearch(data=data[2:], hdf=hdf, path='processed/holdout/nomri',
                   saveall=False, metrics='sqeuclidean')

    # finally, run SNF for combined HC / PD subjects
    data = load_and_residualize_data(hdf)
    run_gridsearch(data=data, hdf=hdf, path='processed/pdhc',
                   saveall=False, metrics='sqeuclidean')


if __name__ == '__main__':
    main()
