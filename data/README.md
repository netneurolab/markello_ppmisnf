# data

This directory contains relevant data required for the reported analyses.

Because the Parkinson's Progression Markers Initiative (PPMI) requires researchers to sign a Data Usage Agreement (DUA), we can't make much of the "raw" data publicly available.
That said, we have, as much as possible, set up our workflows such that as soon as you have access to the PPMI you can automatically run things.

(Check out step two of our [walkthrough](../walkthrough/README.md) for instructions on getting access to the PPMI dataset!)

We include here the derivative neuroimaging data we used in the analyses (so that you don't have to go through the entire procedure of downloading the data and re-processing it yourself):

- [`./raw/ppmi/bids/tesla.csv`](./raw/ppmi/bids/tesla.csv): A CSV file with the scanner strength of the different PPMI participants
- [`./derivative/antslct/qc_*.csv`](./derivative/antslct): Two CSV files containing quality control ratings for the outputs of the ANTs longitudinal cortical thickness pipeline
- [`./derivative/parcellated/*csv`](./derivative/parcellated): NPY and CSV files containing parcellated cortical thickness and subcortical volume measurements, generated via the scripts in [`scripts/01_preprocess`](../scripts/01_preprocess)
