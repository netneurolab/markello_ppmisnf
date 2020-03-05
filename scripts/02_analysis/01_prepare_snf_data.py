# -*- coding: utf-8 -*-
"""
This script, as the name implies, prepares all the "raw" data for SNF. While
the neuroimaging data have already been pre-processed, the remaining data
modalities require a bit of data wrangling (much of which is handled by the
`pypmi` library).

This script loads in all the "raw" data, organizes it into a pseudo-tidy format
that will make analyses easier, performs some batch correction / residualizes
some various regressors, and saves everything into an HDF5 file. (One file is
created for each resolution of the parcellated neuroimaging data.)

Outputs of this script will be stored in the `data/derivative/snf/` directory.
"""

import glob
import itertools
import os.path as op
import re

import numpy as np
import pandas as pd

from netneurotools import datasets as nndata, stats as nnstats
import neurocombat
from ppmi_snf import directories, preprocess, structures
import pypmi


def snfprep(data, demographics, cutoff=0.8, imputation='median'):
    """
    Preprocesses `data` in preparation for SNF

    Preprocessing steps include:

      1. Removing variables / subjects with too many missing values,
      2. Removing outlier subjects (based on median-absolute-deviation), and
      3. Imputing any remaining missing values (median imputation).

    Parameters
    ----------
    data : list of pandas.DataFrame
        All dataframes should be in "tidy" format and the index of the
        dataframe should be the participant ID
    demographics : pandas.DataFrame
        Demographics of participants in `data`
    cutoff : (0, 1) float, optional
        Percent of subjects/variables that must have non-NA values to be
        retained. Default: 0.8
    strategy : {'mean', 'median'}, optional
        Imputation strategy. Default: 'median'

    Returns
    -------
    data, demographics: pd.DataFrame
    """

    # only keep subjects specified by `demographics`
    data = [df.loc[demographics.index] for df in data]

    # drop NaNs and ensure consistent subjects across datatypes
    cleaned = preprocess.clean_data(*data, cutoff=cutoff)

    # get rid of outliers based on median absolute deviation
    # if subject is an outlier in ANY datatype, we want to remove them
    outliers = np.zeros(len(cleaned[0]), dtype=bool)
    for arr in cleaned:
        outliers = np.logical_or(outliers, nnstats.get_mad_outliers(arr))

    # drop outliers (or, as it were, keep subjects that aren't outliers)
    no_outliers = [arr[~outliers] for arr in cleaned]

    # impute data only once outliers have been removed
    imputed = preprocess.impute_data(*no_outliers, strategy=imputation)

    # subset demographics with remaining subjects
    demographics = demographics.loc[imputed[0].index]

    return imputed, demographics


def batch_correct(data, demographics):
    """
    Batch correct `data` with info in `demographics`

    Parameters
    ----------
    data : (S, F) array_like
        Input data where `S` is subjects and `F` is features
    demographics : pandas.DataFrame
        Must have columns ['site', 'family_history', 'gender', 'race',
        'handedness', 'education'] which will be used as regressors

    Returns
    -------
    corrected : numpy.ndarray
        Input `data` corrected for site effects
    """

    combat_kwargs = dict(
        batch_col='site',
        discrete_cols=[
            'diagnosis', 'family_history', 'gender', 'race', 'handedness'
        ],
        continuous_cols=[
            'education'
        ],
        verbose=False
    )
    covars = (combat_kwargs['discrete_cols']
              + combat_kwargs['continuous_cols']
              + [combat_kwargs['batch_col']])
    combat_kwargs['covars'] = demographics.loc[:, covars]

    corrected = neurocombat.combat(data, **combat_kwargs)

    return corrected


def get_visit(df, participants, visit='SC', cutoff=0.8):
    """
    Extracts specified `visit` and `participants` from `df`

    Parameters
    ----------
    df : pd.DataFrame
        "Raw" dataframe that needs a bit of love and attention
    participants : list-of-int
        List of participant IDs to retain in `df`
    visit : str, optional
        Visit to retain in `df`. Default: 'SC'
    cutoff : (0, 1) float, optional
        Percent of subjects/variables that must have non-NA values to be
        retained. Default: 0.8

    Returns
    -------
    df : pd.DataFrame
        Provided data frame after receiving some love and attention
    date : pd.DataFrame
        Visit date extracted from `df`, with same index as `df`
    """

    df = df.dropna(subset=['visit', 'date']) \
           .query(f'participant in {participants} & visit == "{visit}"') \
           .drop('visit', axis=1) \
           .set_index('participant')
    df = df.dropna(axis=1, thresh=(cutoff * len(df)))
    date = pd.DataFrame(df.pop('date'))

    return df, date


def _one_hot(data):
    """
    One-hot encode categorical `data`

    Parameters
    ----------
    data : (N, 1) array_like
        Where `N` is number of samples and there is only one feature

    Returns
    -------
    encoded : (N, F) array_like
        Where `N` is number of samples and `F` is the number of unique values
        in `data`
    """

    data = np.asarray(data).reshape(-1, 1)
    encoded = np.column_stack([data == y for y in np.unique(data)]).astype(int)
    return encoded[:, 1:]


def gen_regressors(date, demo):
    """
    Generates regressor (IV) matrix from `date` and `demo`

    Parameters
    ----------
    date : pd.DataFrame
        Data frame with at least columns ['date']. Participant ID should be
        index of data frame. Can optionally contain column ['etiv']
    demo : pd.DataFrame
        Data frame with at least columns ['date_birth', 'gender']. Participant
        ID should be index of data frame.

    Returns
    -------
    reg : pd.DataFrame
        Data frame with regressors and interaction effects to be modeled and
        residualized from data
    """

    if not all(f in demo.columns for f in ['date_birth', 'gender']):
        raise ValueError('Provided demographics dataframe must have at least '
                         'columns [\'date_birth\', \'gender\'].')
    if 'date' not in date.columns:
        raise ValueError('Provided date dataframe must have \'date\' column.')

    # first, calculate age at time of visit for given datatype
    age = date.loc[demo.index, 'date'] - demo['date_birth']
    age /= np.timedelta64(1, 'Y')

    # then, generate regressor array (including age*gender interaction)
    gender = _one_hot(demo['gender'].get_values())
    age = np.asarray(age)[:, np.newaxis]
    regressors = np.column_stack([age, gender, age * gender])

    # if ETIV is in `date` dataframe, add it to the regressor array
    if 'etiv' in date.columns:
        etiv = np.asarray(date.loc[demo.index, 'etiv'])
        regressors = np.column_stack([regressors, etiv])

    reg = pd.DataFrame(regressors, index=demo.index)

    return reg


def _load_parcels(fname, session=1, parcellation=None, return_date=False):
    """
    Loads parcellated imaging data stored at `fname`

    Parameters
    ----------
    fname : str
        Filepath to .npy file with parcellated data
    session : int, optional
        Session to extract. Default: 1
    parcellation : pandas.DataFrame, optional
        DataFrame with parcellation info. Default: None
    return_date : bool, optional
        Whether to also return dataframe containing participant scan date.
        Default: False

    Returns
    -------
    parcels : (N, G) pandas.DataFrame
        Wide-format parcel data where N is participants and G is parcels
    """

    # load demographic information and get appropriate session info
    demo = pd.read_csv(fname.replace('.npy', '.csv'))
    if session is not None:
        if not isinstance(session, int):
            raise ValueError(f'Provided session {session} is not type int')
        demo = demo.query(f'session == "{session}"')

    # load data and subset if parcellation was provided
    data = np.load(fname)[demo.index]
    if parcellation is not None:
        if isinstance(parcellation, str):
            parcellation = pd.read_csv(parcellation)
        ids = parcellation['id'].astype(int) - 1
        columns = parcellation['label']
        data = data[:, ids]
    else:
        columns = None

    parcels = pd.DataFrame(data, index=demo['participant'], columns=columns)

    if return_date:
        dates = pd.DataFrame({'date': pd.to_datetime(demo['date']).values},
                             index=demo['participant'])
        if 'etiv' in demo.columns:
            dates = dates.assign(etiv=demo['etiv'].get_values().astype(float))

        return parcels, dates

    return parcels


def get_parcels(fname, session=range(1, 5), parcellation=None,
                return_date=False, valid_participants=None):
    """
    Return list of longitudinal data in ``fname`` for ``sessions``

    Parameters
    ----------
    fname : str
        Filepath to .npy file with parcellated data
    session : list, optional
        Session to extract. Default: [1, 2, 3, 4]
    parcellation : pd.core.frame.DataFrame, optional
        DataFrame with parcellation info. Default: None
    return_date : bool, optional
        Whether to also return dataframe containing participant scan date.
        Default: False
    valid_participants : list, optional
        Which participants to retain. Default: all participants

    Returns
    -------
    parcels : list of pandas.DataFrame
        Longitudinal parcellated data
    dates : list of pandas.DataFrame
        Acquisition dates of data in `parcels`
    """

    # grab all sessions if we didn't specify one
    if session is None:
        session = [1, 2, 3, 4]
    elif isinstance(session, int):
        session = [session]

    # load the parcellated data and update the `date` dataframe
    parcels, dates = [], []
    for ses in session:
        data, date = _load_parcels(fname,
                                   session=ses,
                                   parcellation=parcellation,
                                   return_date=True)
        if valid_participants is not None:
            date = date[date.index.isin(valid_participants)]
            data = data[data.index.isin(valid_participants)]
        parcels.append(data)
        dates.append(date)

    # don't return a length-one list
    if len(parcels) == 1:
        parcels, dates = parcels[0], dates[0]

    if return_date:
        return parcels, dates

    return parcels


def main():
    # N.B. this will NOT work unless you set the environmental variables
    #      $PPMI_USER and $PPMI_PASSWORD prior to running this script.
    #      these variables must be the username and password you received when
    #      registering for the PPMI. for more information on data access see:
    #      https://www.ppmi-info.org/access-data-specimens/download-data/
    pypmi.fetch_studydata('all', path=directories.ppmi, overwrite=False)

    # load demographic data and keep only individuals with PD and healthy
    # individuals. we'll use the information in this data frame to residualize
    # our data against different variables (e.g., age, gender)
    print('Loading demographics information...')
    demographics = pypmi.load_demographics(directories.ppmi) \
                        .query('diagnosis in ["pd", "hc"]') \
                        .set_index('participant')
    demographics['family_history'] = demographics['family_history'].astype(bool)

    # load all non-MRI data
    print('Loading all non-MRI data (this step may take some time)...')
    datscan = pypmi.load_datscan(directories.ppmi, measures='all')
    biospec = pypmi.load_biospecimen(directories.ppmi, measures='all')
    behavior = pypmi.load_behavior(directories.ppmi, measures='all')

    # sometimes, because of how PPMI data were collected, there are slight
    # variations in the recorded date for the same visit, resulting in scores
    # for a single visit being split across two or more rows in the dataframe
    # (i.e., one row might have MoCA scores for visit "V01" and the other has
    # UPDRS scores for visit "V01")
    # to remedy this we use pandas `DataFrame.combine_first()` method, merging
    # scores from both rows and retaining the earliest date as the "true" date
    # (dates were generally only ~1 month different and if that difference
    # makes a significant impact on our results then I quit)
    print('Wrangling non-MRI data into a usable format...')
    first = behavior.drop_duplicates(['participant', 'visit'], 'first') \
                    .reset_index(drop=True)
    last = behavior.drop_duplicates(['participant', 'visit'], 'last') \
                   .reset_index(drop=True)
    behavior = first.combine_first(last)

    # get first visit scores for non-MRI data
    datscan, dat_date = get_visit(datscan, list(demographics.index), visit='SC')
    biospec, bio_date = get_visit(biospec, list(demographics.index), visit='BL')

    # behavioral data acquisition was split across screening + baseline visits
    # so we need to take the earliest visit for each measure
    # that is, not all measures were collected at screening so we need to use
    # the baseline visit scores for those measures
    # unfortunately which visit various measures were initially collected at
    # DIFFERED for PD and HC individuals, so we need to do this separately for
    # the two groups and then merge them back together... ¯\_(ツ)_/¯
    beh, beh_dates = [], []
    for diagnosis in ['pd', 'hc']:
        participants = demographics.query(f'diagnosis == "{diagnosis}"').index
        beh_sc, beh_date = get_visit(behavior, list(participants), visit='SC')
        beh_bl, _ = get_visit(behavior, list(participants), visit='BL')
        drop = np.intersect1d(beh_sc.columns, beh_bl.columns)
        beh += [pd.merge(beh_sc, beh_bl.drop(drop, axis=1), on='participant')]
        beh_dates += [beh_date]
    behavior = pd.concat(beh, join='inner')
    beh_date = pd.concat(beh_dates, join='inner')

    # iterate through all combinations of cortical + subcortical parcellations
    # note: there's only one subcortical parcellation (we had considered doing
    # more but the number of good subcortical parcellations is...limited)
    cth_data = sorted(glob.glob(op.join(directories.parcels, '*thickness.npy')))
    vol_data = sorted(glob.glob(op.join(directories.parcels, '*volume.npy')))
    for cth, vol in itertools.product(cth_data, vol_data):

        # determine what cortical / subcortical parcellation combo we're using
        # this will determine the name of the output file
        # the specific details include the resolution of cortical parcellation
        # and the datatype of the subcortical parcellation
        (scale, ) = re.search(r'res-(\d+)', cth).groups()
        (dtype, ) = re.search(r'_hemi-both_(\S+)_', vol).groups()
        hdf = structures.Frog(op.join(directories.snf,
                                      f'scale{scale}_{dtype}.h5'))
        print(f'Loading MRI data for {op.basename(hdf.filename)}...')

        # load parcellated cortical thickness data
        ct_parc = nndata.fetch_cammoun2012(data_dir=directories.rois,
                                           verbose=0)['info']
        ct_parc = pd.read_csv(ct_parc).query(f'scale == "scale{scale}" '
                                             '& structure == "cortex"')
        ct_parc['label'] = (ct_parc['label'] + '_'
                            + ct_parc['hemisphere'].apply(str.lower))
        cortthick, cth_date = get_parcels(cth, session=1, return_date=True,
                                          parcellation=ct_parc)

        # load parcellated subcortical volume data
        sv_parc = nndata.fetch_pauli2018(data_dir=directories.rois,
                                         verbose=0)['info']
        sv_parc = pd.read_csv(sv_parc)
        subvolume, vol_date = get_parcels(vol, session=1, return_date=True,
                                          parcellation=sv_parc)

        # perform batch correction on MRI data
        # first, grab the demographics of subjects for whom we have neuro data.
        # then, remove all sites where we only have data from one subject since
        # we cannot generate batch correction parameters in these instances.
        # finally, perform the actual batch correction using `neurocombat`
        cortthick, subvolume, demo = \
            preprocess.intersect_subjects(cortthick, subvolume, demographics)
        sites, counts = np.unique(demo['site'], return_counts=True)
        demo = demo[demo['site'].isin(sites[counts > 1])]
        cortthick, subvolume, demo = \
            preprocess.intersect_subjects(cortthick, subvolume, demo)
        cortthick.iloc[:, :] = batch_correct(cortthick, demo)
        subvolume.iloc[:, :] = batch_correct(subvolume, demo)

        # only keep subjects for whom we have all datatypes
        # we preprocess HC and PD data separately because part of the process
        # involves imputation and we want to impute missing data using values
        # from each diagnostic group, separately
        data = [cortthick, subvolume, datscan, biospec, behavior]
        *data, demo = preprocess.intersect_subjects(*data, demo)
        hc_data, hc_demo = snfprep(data, demo.query('diagnosis == "hc"'))
        pd_data, pd_demo = snfprep(data, demo.query('diagnosis == "pd"'))

        # only keep features for which we have both PD and HC data
        for n, (hc_dtype, pd_dtype) in enumerate(zip(hc_data, pd_data)):
            cols = np.intersect1d(hc_dtype.columns, pd_dtype.columns)
            hc_data[n], pd_data[n] = hc_data[n][cols], pd_data[n][cols]

        # "regress out" age, gender, age x gender interactions (and total
        # estimated intracranial volume, if MRI data) from all data.
        # we also want to save all this data to disk so we can load it easily
        # in the future! do that for all the raw data, regressor matrices, and
        # processed (i.e., residualized) data
        # we do this because we don't want these sorts of things to bias our
        # initial analyses when creating the fused networks
        keys = [
            'cortical_thickness',
            'subcortical_volume',
            'dat_scans',
            'csf_assays',
            'behavioral_measures'
        ]
        dates = [cth_date, vol_date, dat_date, bio_date, beh_date]
        for grp, dataset, demo in zip(['pd', 'hc'],
                                      [pd_data, hc_data],
                                      [pd_demo, hc_demo]):
            hdf.save(demo, f'/raw/{grp}_demographics', overwrite=False)
            for n, (df, key, date) in enumerate(zip(dataset, keys, dates)):
                reg = gen_regressors(date, demo)

                # get comparative regressors / data (this is always healthy
                # inviduals -- we use them to estimate the betas for the
                # residualization process)
                comp_reg, comp_df = gen_regressors(date, hc_demo), hc_data[n]

                resid = nnstats.residualize(reg, df, comp_reg, comp_df,
                                            normalize=False)
                resid = pd.DataFrame(resid, index=df.index, columns=df.columns)

                hdf.save(df, f'/raw/{grp}_{key}', overwrite=False)
                hdf.save(reg, f'/regressors/{grp}_{key}', overwrite=False)
                hdf.save(resid, f'/processed/{grp}_{key}', overwrite=False)


if __name__ == '__main__':
    main()
