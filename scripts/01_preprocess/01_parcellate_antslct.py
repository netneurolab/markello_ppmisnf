# -*- coding: utf-8 -*-
"""
This script takes the outputs of the ANTs longitudinal cortical thickness
pipeline (hence, antslct) and parcellates the data, generating cortical
thickness and subcortical volume measurements.

Obviously, this assumes that the ANTs longitudinal cortical thickness pipeline
has been run on the data. Regardless, the outputs of this script should land
in the `data/derivative/parcellated` directory.
"""

import glob
import os
import os.path as op
import re
import time

import nibabel as nib
from nilearn._utils import check_niimg
import numpy as np
import pandas as pd
from scipy import ndimage

from netneurotools import datasets, utils
from ppmi_snf import directories

# filename template for subject-specific warps to warp masks to subject space
_wrp = 'sub-{sub}_ses-{ses}_run-01_T1w_templatetosubjectwarp.nii.gz'
_seg = 'sub-{sub}_ses-{ses}_run-01_T1w_brainsegmentation.nii.gz'
WRP_FNAME = op.join(directories.ants, 'sub-{sub}', _wrp)
SEG_FNAME = op.join(directories.ants, 'sub-{sub}', _seg)

# filename templates for extracting subject / session info
TSV_FNAME = op.join(directories.bids, 'sub-{sub}', 'ses-{ses}',
                    'sub-{sub}_ses-{ses}_scans.tsv')
T1W_FNAME = op.join('anat', 'sub-{sub}_ses-{ses}_run-01_T1w.nii.gz')

REGEX = re.compile(r'sub-(\d+)_ses-(\d+)')
DFMT = '%Y-%m-%d'

# QC files with at least columns: participant, session, seg_rating, reg_rating
QC_RATER1 = op.join(directories.ants, 'qc_rm.csv')
QC_RATER2 = op.join(directories.ants, 'qc_ct.csv')


def warp_to_subject(infile, warpfile, outfile=None, template=None,
                    interpolation='NearestNeighbor', verbose=False,
                    use_cache=True):
    """
    Applys `warpfile` to `infile`

    Just a wrapper around antsApplyTransforms to be called from within Python;
    literally opens a subprocess to run the command

    Parameters
    ----------
    infile : str
        Image to be warped
    warpfile : str
        ANTs-generated warp file to apply to `infile`
    outfile : str, optional
        Name of output file. If not specified, will be determined from provided
        `warpfile` and `infile`. Default: None
    template : str, optional
        Reference image that specifies shape, dimensions, etc to be generated
        when appling `warpfile` to `infile`. If not specified the
        defaults embedded in `warpfile` will be used. Default: None
    interpolation : str, optional
        Type of interpolation to use during warping. Default: 'NearestNeighbor'
    verbose : bool, optional
        Whether to print status messages as transformation is applied. Default:
        False
    use_cache : bool, optional
        Whether to check for existence of `outfile` and use that, if it exists.
        If False, will create a new `outfile` regardless of existence. Default:
        True

    Returns
    -------
    outfile : str
        Path to warped parcellation
    """

    warpcmd = 'antsApplyTransforms -d 3 -e {imagetype} {opts} -v {verbose} ' \
              '-n {interpolation} -i {input} -o {output} -t {warpfile}'

    # scalar or timeseries, depending on the dimensionality of image (3 vs 4)
    imagetype = [0, 3][nib.load(infile).ndim > 3]

    opts = ''
    if template is not None:
        opts += f'-r {template}'

    if outfile is None:
        outfile = op.join(op.dirname(warpfile), op.basename(infile))

    if not op.isfile(outfile) or not use_cache:
        utils.run(warpcmd.format(imagetype=imagetype,
                                 opts=opts,
                                 verbose=int(verbose),
                                 interpolation=interpolation,
                                 input=infile,
                                 output=outfile,
                                 warpfile=warpfile),
                  quiet=not verbose)

    return outfile


def get_data(sub, mask, dtype, extractfunc, sessions=None, warptosub=False,
             interpolation='NearestNeighbor'):
    """
    Extracts data for `sub`

    Parameters
    ----------
    sub : str
        Subject ID from which to extract data
    mask : str
        Path to file for masking
    dtype : str
        Data type from which to extract data
    extractfunc : func
        Function for extracting data; must accept two parameters (mask, img)
    sessions : list of str, optional
        List of session identifiers from which to extract data. If not
        specified all available sessions will be used. Default: None
    warptosub : bool, optional
        Whether to warp provided `mask` to subject-space images using available
        warps files. Default: False
    interpolation : str, optional
        Type of interpolation to use if `warptosub=True`. Default:
        'NearestNeighbor'

    Returns
    -------
    data : (N, M) numpy.ndarray
        Data extracted from imgs, where `N` is the number of images and `M` is
        the number of parcels in `mask`
    demographics : (N, 4) pandas.DataFrame

    """

    # get images for supplied sessions (or all images)
    subj_dir = op.join(directories.ants, sub)
    if sessions is not None:
        imgs = []
        for ses in sorted(sessions):
            sespath = op.join(subj_dir, f'{sub}_ses-{ses}_*{dtype}.nii.gz')
            imgs.extend(glob.glob(sespath))
    else:
        imgs = sorted(glob.glob(op.join(subj_dir, f'*{dtype}.nii.gz')))

    if len(imgs) == 0:
        return None, None

    # get subject / session info and add date/etiv for all sessions
    sessions = np.row_stack([REGEX.findall(i) for i in imgs])
    demographics = _add_info(sessions)

    # warp mask to subject images, if needed
    if warptosub:
        outfile = 'sub-{0}_ses-{1}_run-01_T1w_' + op.basename(mask)
        warps = [WRP_FNAME.format(sub=sub, ses=ses) for (sub, ses) in sessions]
        masks = []
        for img, warp in zip(imgs, warps):
            out = op.join(subj_dir, outfile.format(*REGEX.findall(warp)[0]))
            masks.append(warp_to_subject(mask, warp, out, template=img,
                                         interpolation=interpolation))
    else:
        masks = np.repeat([mask], len(imgs))

    # fit mask to data and stack across sessions
    data = np.row_stack([extractfunc(check_niimg(img, atleast_4d=True), mask)
                        for mask, img in zip(masks, imgs)])

    # remove masks warped to subject space (to save on disk space)
    if warptosub:
        for mask in masks:
            os.remove(mask)

    return data, demographics


def _add_info(demo):
    """
    Adds scan date / estimated total intracranial volume to `demo`

    Parameters
    ----------
    demo : array_like
        First column should be participant ID, second should be session ID

    Returns
    -------
    info : pandas.DataFrame
        `demo` converted to DataFrame and with added 'DATES' and 'ETIV' columns
    """

    info = pd.DataFrame({'participant': demo[:, 0],
                         'session': demo[:, 1],
                         'date': '',
                         'etiv': ''})
    for n, (sub, ses) in enumerate(demo):
        fmt = dict(sub=sub, ses=ses)

        # get date of scan
        acqs = pd.read_csv(TSV_FNAME.format(**fmt), sep='\t')
        t1w = acqs.query('filename == "{}"'.format(T1W_FNAME.format(**fmt)))
        try:
            acqtime = pd.to_datetime(t1w.acq_time.squeeze()).strftime(DFMT)
        except (pd.errors.OutOfBoundsDatetime, AttributeError):
            raise Exception(f'Error finding scan date for {sub} session {ses}')
        info.loc[n, 'date'] = acqtime

        # get estimated total intracranial volume of subject at scan
        try:
            seg = nib.load(SEG_FNAME.format(**fmt))
            voxel_size = np.prod(seg.header.get_zooms())
            etiv = voxel_size * np.sum(seg.get_data() != 0)
        except FileNotFoundError:
            etiv = np.nan
        info.loc[n, 'etiv'] = etiv

    return info


def get_labels(infile, label_img):
    """
    Averages data in `infile` for each region in `label_img`

    Parameters
    ----------
    infile : str
        File with data to extract
    label_img : str
        File with regions to define extraction

    Returns
    -------
    signal : numpy.ndarray
        Averaged data in `infile` defined by regions in `label_img`
    """

    infile = check_niimg(infile, atleast_4d=True)
    label_img = check_niimg(label_img, ensure_ndim=3)
    data = infile.get_data()

    if infile.shape[:3] != label_img.shape:
        raise ValueError('Provided file {} does not match label image {}.'
                         .format(infile.file_map['image'].filename,
                                 label_img.file_map['image'].filename))

    labels_data = label_img.get_data()
    labels = np.trim_zeros(np.unique(labels_data))

    # average over data within each parcel
    signal = np.asarray(ndimage.measurements.mean(data.squeeze(),
                                                  labels=labels_data,
                                                  index=labels))
    return signal


def get_volume(infile, threshold=0.4):
    """
    Calculates volumes of regions defined by `infile`, in mm^3

    Parameters
    ----------
    infile : str
        File with regions to calculate volume
    threshold : [0, 1] float, optional
        Threshold to apply to `infile` before volume calculation. Only used if
        `infile` is 4D (i.e., has probabilistic parcels). Default: 0.4

    Returns
    -------
    out : numpy.ndarray
        Volumes of regions defined by `infile`
    """

    infile = check_niimg(infile, atleast_4d=True)
    voxsize = np.prod(infile.header.get_zooms()[:3])
    data = infile.get_data()

    if infile.shape[-1] > 1:
        # probabilistic atlas so we want to threshold the data and then sum
        # over the first three axes
        # the fourth axis is the different parcels, which we want to keep
        out = np.sum(data >= threshold, axis=(0, 1, 2), dtype=float)
    else:
        # deterministic atlas so we want to sum within each unique parcel
        data = data.squeeze()
        indices = np.trim_zeros(np.unique(data))
        out = np.asarray(ndimage.measurements.sum(data != 0,
                                                  labels=data,
                                                  index=indices))

    # multiply by voxel size to get mm^3
    out *= voxsize

    return out


def get_subjects(qc, return_record=True, min_seg=1, min_reg=1):
    """
    Gets subjects / sessions from `qc` with ratings above `min_seg` & `min_reg`

    Parameters
    ----------
    qc : str
        Filepath to QC information
    return_record : bool, optional
        Whether to return numpy recarray instead of ndarray
    min_seg : {0, 1, 2}, optional
        Minimum segmentation rating for subject/session to be kept. Default: 1
    min_reg : {0, 1, 2}, optional
        Minimum registration rating for subject/session to be kept. Default: 1

    Returns
    -------
    data : numpy.ndarray
        Record array with `participant` and `session` fields for "good" data
    """

    df = pd.read_csv(qc)
    # only keep ratings above minimum, drop NA values, and drop SST sessions
    data = df.query(f'seg_rating >= {min_seg} & reg_rating >= {min_reg}') \
             .dropna(subset=['seg_rating', 'reg_rating']) \
             .query('session != "SST"') \
             .get(['participant', 'session'])
    data = data.to_records(index=False) if return_record else data.get_values()

    return data


def extract_data(mask, dtype, outfile):
    """
    Extracts data for `dtype` using `mask` and saves to `outfile`

    Parameters
    ----------
    mask : str
        Mask to use for parcellating data
    dtype : {'brainsegmentation', 'sscorticalthickness'}
        Datatype from which data should be extracted
    outfile : str
        Path to desired output file. A .npy file will be created with the
        parcellated data and a .csv file will be created with relevant
        demographic information (e.g., participant, session, date)
    """

    # check that specified datatype is one of the two we're expecting
    if dtype not in ['brainsegmentation', 'sscorticalthickness']:
        raise ValueError('Invalid `dtype` {}'.format(dtype))

    # data extraction parameters depend on datatype
    info = dict(
        brainsegmentation=dict(
            # whether we should be masking data in subject space (i.e., do we
            # need to warp the mask to the subject)
            warptosub=True,
            # type of interpolation to apply to `mask` if warping to subject
            interpolation='Linear',
            # how to extract data from img + mask once they're in same space
            # extractfunc _must_ accept two inputs (image, mask) but here we
            # only want to operate on the mask, so we write a stupid lambda
            # function around get_volume() to simply ignore the image input
            extractfunc=lambda image, mask: get_volume(mask)
        ),
        sscorticalthickness=dict(
            warptosub=True,
            interpolation='NearestNeighbor',
            extractfunc=get_labels
        )
    )

    # empty list to store parcellated data
    data = []

    # empty dataframe to store auxiliary info for parcellation data
    demo = pd.DataFrame(columns=['participant', 'session', 'date', 'etiv'])

    # if output files exist, load them in for checkpointing
    if op.exists(outfile + '.npy'):
        data = np.load(outfile + '.npy').tolist()
    if op.exists(outfile + '.csv'):
        demo = pd.read_csv(outfile + '.csv')

    # get record of unique subject/session pairs that have ANTs derivatives
    rec = np.intersect1d(get_subjects(QC_RATER1, return_record=True),
                         get_subjects(QC_RATER2, return_record=True))

    # create output directory, if it doesn't already exist
    os.makedirs(directories.parcels, exist_ok=True)

    # iterate through subjects and mask data
    for sub in np.unique(rec['participant']):
        # grab all sessions for this subject
        sessions = rec[rec['participant'] == sub]['session']

        # brief status update about who we're extracting data from
        print(f'Extracting data for {sub}, {len(sessions)} session(s)')
        start_time = time.time()

        # skip subject if they've already been processed to save time
        proc = demo.query(f'participant == {sub.split("-")[-1]}') \
                   .get('session') \
                   .astype(str)
        ses = np.setdiff1d(sessions, np.asarray(proc))
        if len(ses) == 0:
            continue

        # get data for all (remaining) sessions for subject
        cdata, cdemo = get_data(sub, mask=mask, dtype=dtype, sessions=ses,
                                **info.get(dtype))

        # if we found data, add it to group arrays / dataframes and
        # checkpoint by saving the data to disk
        if cdata is not None and cdemo is not None:
            # add data to group arrays / dataframes
            data.append(cdata)
            demo = demo.append(cdemo, ignore_index=True)

            # save to disk
            np.save(outfile + '.npy', np.row_stack(data))
            demo.to_csv(outfile + '.csv', index=False)
        else:
            print(f'Warning: Cannot find data for {sub} session(s) {ses}')

        # status update
        elapsed = time.time() - start_time
        print(f'Extraction took {elapsed:0.2f} seconds')


def main():
    # cammoun2012 is used to generate cortical thickness measures
    cammoun2012 = datasets.fetch_cammoun2012(data_dir=directories.rois)
    del cammoun2012['info']
    for scale, mask in cammoun2012.items():
        outfile = op.basename(mask).replace('.nii.gz', '_corticalthickness')
        outpath = op.join(directories.parcels, outfile)
        extract_data(mask, 'sscorticalthickness', outpath)

    # pauli2018 is used to generate subcortical volume measures
    pauli2018 = datasets.fetch_pauli2018(data_dir=directories.rois)
    del pauli2018['info'], pauli2018['probabilistic']
    for scale, mask in pauli2018.items():
        outfile = op.basename(mask).replace('.nii.gz', '_subcorticalvolume')
        outpath = op.join(directories.parcels, outfile)
        extract_data(mask, 'brainsegmentation', outpath)


if __name__ == "__main__":
    main()
