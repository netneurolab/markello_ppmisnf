# -*- coding: utf-8 -*-
"""
Small helper functions for preprocessing messy data
"""

import pandas as pd


def intersect_subjects(*inputs):
    """
    Finds intersection of participants between `inputs`

    Parameters
    ----------
    inputs : pandas.DataFrame
        Index of dataframes will be used as identifier

    Returns
    -------
    outputs : (N,) list of pandas.DataFrame
        Input dataframes with only shared participants retained
    """

    # get intersection of subjects (row indices of input dataframes)
    subjects = list(set.intersection(*[set(f.index) for f in inputs]))

    # subset dataframes with overlapping subjects
    outputs = [f[f.index.isin(subjects)].copy().sort_index() for f in inputs]

    return outputs


def clean_data(*inputs, cutoff=0.8):
    """
    Parameters
    ----------
    inputs : pandas.DataFrame
        Input dataframes (all same length!)
    cutoff : (0, 1) float, optional
        Percent of subjects/variables that must have non-NA values to be
        retained. Default: 0.8

    Returns
    -------
    inputs : list of pd.core.frame.DataFrame
        Cleaned inputs
    """

    if cutoff < 0 or cutoff > 1:
        raise ValueError(f'Supplied cutoff {cutoff} not between [0, 1].')

    inputs = intersect_subjects(*inputs)

    inputs = [f.dropna(thresh=cutoff * f.shape[0], axis=1) for f in inputs]
    inputs = [f.dropna(thresh=cutoff * f.shape[1], axis=0) for f in inputs]

    inputs = intersect_subjects(*inputs)

    return inputs


def impute_data(*inputs, strategy='median'):
    """
    Imputes missing data in ``inputs``

    Parameters
    ----------
    inputs : list of pandas.DataFrame
        List of input data
    strategy : {'mean', 'median'}, optional
        Imputation strategy. Default: 'median'

    Returns
    -------
    ouputs : list of pandas.DataFrame
        List of imputed data
    """

    from sklearn.impute import SimpleImputer

    # https://github.com/scikit-learn/scikit-learn/pull/9212
    impute = SimpleImputer(strategy=strategy)
    outputs = [pd.DataFrame(impute.fit_transform(f), columns=f.columns,
                            index=f.index) for f in inputs]

    return outputs
