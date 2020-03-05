# BIDSifying the PPMI

The [Brain Imaging Data Structure (BIDS)](https://bids.neuroimaging.io/) provides a standardized format for naming and orgnizing neuroimaging data.
It's *incredibly* helpful, and if you haven't heard about it I highly recommend you go and read up on it now.
(Check out the fantastic [BIDS Starter Kit](https://github.com/bids-standard/bids-starter-kit) resource to see how you can get started with BIDS on your datasets.)

We re-organized the PPMI data into BIDS format so that we could process / wrangle it a bit more easily and automate processing via [BIDS Apps](https://doi.org/10.1371/journal.pcbi.1005209).
While we have provided the outputs of the neuroimaging analyses such that you don't **need** to re-process the data yourself, if you're interested in doing so for whatever reason you should be able to by following the instructions below.

## PPMI to BIDS

Assuming you followed the [previous steps](01_accessing_data.md), you should have downloaded the raw neuroimaging data from the PPMI and unzipped it into the `data/raw/ppmi/dicoms` directory.
(You can put the data wherever you like; you'll just have to update the relevant path below.)

We rely on the amazing [`heudiconv`]() library to convert the raw DICOM images to a BIDS-compatible dataset.
We've tried to make this BIDSification process for the PPMI dataset as easy as possible by wrapping it into the [`pypmi`](https://github.com/rmarkello/pypmi) package with the following command:

```python
import pypmi.bids
pypmi.bids.convert_ppmi(raw_dir='data/raw/ppmi/dicoms', 
                        out_dir='data/raw/ppmi/bids', 
                        ignore_bad=True
                        coerce_study_uids=True)
```

To learn a bit more about what the function does behind-the-scenes you can read on below.
Alternatively, jump to the [next section](./03_antslct_pipeline.md) to begin re-processing the PPMI data.

## The `convert_ppmi()` function

The `pypmi.bids.convert_ppmi` function does some light re-organization of the input `raw_dir` directory and then uses `heudiconv` to convert the DICOMs to BIDS format.

The function assumes the structure of `raw_dir` is something like the following:

```
└── raw_dir/
    ├── SUB-0001/
    |   ├── SCAN-TYPE-1/                      (e.g., MPRAGE_GRAPPA)
    |   |   ├── SCAN-DATETIME-1/
    |   |   |   └── SCAN-SERIES-1/
    |   |   |       └── *dcm
    |   |   └── SCAN-DATETIME-2/
    |   ├── SCAN-TYPE-2/                      (e.g., AX_FLAIR)
    |   |   ├── SCAN-DATETIME-1/
    |   |   |   └── SCAN-SERIES-2/
    |   |   |       └── *dcm
    |   |   └── SCAN-DATETIME-2/
    |   ├── .../
    |   └── SCAN-TYPE-N/
    ├── SUB-0002/
    ├── .../
    └── SUB-NNNN/
```

where `SUB-NNNN` are PPMI subject numbers, `SCAN-DATETIME-N` are timestamps of the format YYYY-MM-DD_HH_MM_SS.0, and `SCAN-SERIES-N` are unique identifiers of the format S######.
(Note: you should not have to do any re-naming or re-organization; this is how the PPMI organizes the data if you download it according to the previous instructions.)

Unfortunately, this "default" sub-directory structure is not conducive to use with `heudiconv`.
Scans are grouped by scan type rather than by session, and there are a number of redundant sub-directories that we don't need.
To resolve this, the function reorganizes the data, moving scans around so that the general hierarchy is `{subject}/{session}/{scan}`, which makes for a much easier time converting the PPMI dataset into BIDS format.

An added complication is that a minority of the scans in the PPMI database are "bad" to some degree.
For most, it is likely that there was some issue with exporting/uploading the DICOM files.
For others, the conversion process we intend to utilize (`heudiconv` and `dcm2niix`) fails to appropriately convert the files due to some idiosyncratic reason (that might be fixed at this point, but was not fixed when we originally converted the data).
These scans need to be removed so that we can run the batch of subjects through `heudiconv` without any abrupt failures.
By default, these scans are moved to a sub-directory of `raw_dir`; setting `ignore_bad` to False will retain these scans and try to convert them (but be warned!).

Once re-organization is done the resulting directory is processed with `heudiconv` and the converted BIDS dataset is stored in `out_dir`.

**NOTE:** this function uses Docker to run `heudiconv`, so make sure you followed the [installation instructions](./00_setting_up.md) before running it!

---

[Click here to continue the walkthrough](./03_antslct_pipeline.md)
