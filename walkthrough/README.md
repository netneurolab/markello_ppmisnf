# walkthrough

The Markdown documents contained in this directory document the various steps that you can follow to (hopefully) re-produce our analyses and results.
You can use the table of contents below to jump around, or [click here](./00_setting_up.md) to just start at the beginning.

## Table of Contents

* [**Step 1**: Set-up and installation](./00_setting_up.md)
  * [Required software](./00_setting_up.md#required-software)
  * [Getting the repository with `git`](./00_setting_up.md#getting-the-repository-with-`git`)
  * [Python dependencies](./00_setting_up.md#python-dependencies)
    * [Using `conda` (recommended)](./00_setting_up.md#using-`conda`-(recommended))
    * [Using `pip`](./00_setting_up.md#using-`pip`)
    * [Internal libraries](./00_setting_up.md#internal-libraries)
  * [Docker and Singularity](./00_setting_up.md#docker-and-singularity)
* [**Step 2**: Downloading data from the PPMI](./01_accessing_data.md)
  * [Getting access to the data](./01_accessing_data.md#getting-access-to-the-data)
  * [Downloading the neuroimaging data](./01_accessing_data.md#downloading-the-neuroimaging-data)
  * [Downloading the rest of the data](./01_accessing_data.md#downloading-the-rest-of-the-data)
* [**Step 3**: BIDSifying the PPMI](./02_converting_to_BIDS.md)
  * [PPMI to BIDS](./02_converting_to_BIDS.md#ppmi-to-bids)
  * [The `convert_ppmi()` function](./02_converting_to_BIDS.md#the-`convert_ppmi()`-function)
* [**Step 4**: Processing the PPMI neuroimaging data](./03_antslct_pipeline.md)
  * [Re-running the ANTs pipeline](./03_antslct_pipeline.md#re-running-the-ants-pipeline)
  * [The modified ANTs BIDS App](./03_antslct_pipeline.md#the-modified-ants-bids-app)
  * [Quality control](./03_antslct_pipeline.md#quality-control)
  * [Parcellating the data](./03_antslct_pipeline.md#parcellating-the-data)
* [**Step 5**: Reproducing the analyses](./04_snf_analyses.md)
  * [Preparing for SNF](./04_snf_analyses.md#preparing-for-snf)
  * [The SNF procedure](./04_snf_analyses.md#the-snf-procedure)
  * [Data concatenation versus SNF](./04_snf_analyses.md#data-concatenation-versus-snf)
  * [PD patient biotypes](./04_snf_analyses.md#pd-patient-biotypes)
  * [MRI contributions](./04_snf_analyses.md#mri-contributions)
  * [Diffusion embedding](./04_snf_analyses.md#diffusion-embedding)
  * [Supplementary analyses](./04_snf_analyses.md#supplementary-analyses)

If you have any questions, simply [open an issue](https://github.com/netneurolab/markello_ppmisnf/issues) on this repository and someone will get back to you as soon as possible!
