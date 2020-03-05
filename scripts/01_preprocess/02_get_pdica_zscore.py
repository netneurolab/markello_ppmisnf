# -*- coding: utf-8 -*-
"""
This script calculates a PD-ICA atrophy score for each subject, originally
described in Zeighami et al., 2015, eLife.

The atrophy score is a weighted average of the jacobian determinant images for
each PPMI subject, using the "PD-ICA" Z-map provided by the aforementioned
authors via NeuroVault as the mask. This weighted average serves as a proxy for
the PD-specific atrophy score in the absence of the mixing matrix generated
from the original ICA decomposition (e.g., Zeighami et al., 2019, NeuroImage:
Clinical).

Weighted averages are saved to the `data/derivative/parcellated` directory for
future use.
"""

import glob
import os.path as op
import re

import nibabel as nib
from nilearn.datasets import fetch_neurovault_ids
from nilearn.image import resample_to_img
import pandas as pd

from ppmi_snf import directories


def apply_weighted_mask(img, mask):
    """
    Returns weighted average of `mask` applied to `img`

    Averaging procedure ignores voxels which are 0 in `mask`

    Parameters
    ----------
    img : niimg_like
        Image to which `mask` should be appplied
    mask : niimg_like
        Weighted (not binarized) mask to apply to `img`

    Returns
    -------
    average : float
        Weighted average of `mask` applied to `img`
    """

    mdata = mask.get_fdata()
    res = resample_to_img(img, mask).get_fdata() * mdata

    return res[mdata != 0].mean()


def main():
    # fetch and load the PD-ICA map from Zeighami et al., 2015, eLife
    pdica = fetch_neurovault_ids(image_ids=[12551], data_dir=directories.rois)
    pdica = nib.load(pdica['images'][0])

    # calculate PD atrophy score for each subject/session
    atrophy = []
    jacobians = op.join(directories.ants, 'sub-????', '*invjacobian.nii.gz')
    for img in sorted(glob.glob(jacobians)):

        # extract the subject and session IDs from each file
        sub = int(re.search(r'/sub-(\d+)/', img).group(1))
        ses = int(re.search(r'ses-(\d+)', img).group(1))

        # some of the jacobians were created incorrectly so we receive an
        # EOFError because of the gzip nonsense. what can you do? (answer:
        # re-run those subjects that failed)
        try:
            print(f'Applying PD-ICA mask to subject {sub}, session {ses}')
            atrophy.append([sub, ses, apply_weighted_mask(img, pdica)])
        except EOFError:
            pass

    # generate little CSV file with information and save to disk
    columns = ['participant', 'session', 'atrophy']
    fname = op.join(directories.parcels, 'pdica_atrophy.csv')
    pd.DataFrame(atrophy, columns=columns).to_csv(fname, index=False)


if __name__ == '__main__':
    main()
