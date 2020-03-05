# code

This directory contains Python libraries used in the reported analyses.
Some of these libraries were developed by the authors explicitly in support of the current project (e.g., [`pypmi`](https://github.com/rmarkello/pypmi) and [`snfpy`](https://github.com/rmarkello/snfpy)), but were separated into standalone projects for easier development and reusability.
Others were borrowed from other groups (e.g., [`neurocombat`](https://github.com/ncullen93/neuroCombat)) or represent more general-use libraries (e.g., [`netneurotools`](https://github.com/netneurolab/netneurotools)).
The remaining project-specific codebits are included as an internal module ([`ppmi_snf`](./ppmi_snf)), containing various project settings and utilities used throughout the analyses.

We briefly describe each of these libraries and their relevance to the current project below.

## Project-specific code

### [`ppmi_snf`](./ppmi_snf)

Contains various settings (e.g., see `ppmi_snf.defaults` and `ppmi_snf.directories`) and functions used throughout the current project.
These codebits aren't necessarily very generalizable, but are re-used at various points throughout preprocessing, analysis, and results generation so they exist here rather than in any one single script.

(N.B. The potentially most reusable bit of code in here is `ppmi_snf.structures.Frog`&mdash;named in honor of my dog, [Frog](https://www.instagram.com/frogandkit/)&mdash; which is similar in spirit, though not robustness, to [`deepdish`](https://github.com/uchicago-cs/deepdish). I just didn't know that package existed when beginning this project so this was the solution I came up with. We live and we learn.)

## Primary analytic libraries

### [`snfpy`](./snfpy)

This library is a Python translation of the codebase implementing [similarity network fusion](https://doi.org/10.1038/nmeth.2810) (SNF), originally released in [MATLAB and R](http://compbio.cs.toronto.edu/SNF/SNF/Software.html) by Bo Wang and colleagues.
It largely recapitulates the functionality of the original SNF packages, but adds a few new features, including some code for conducting the hyperparameter gridsearch reported in the main analyses.
If you're interested in applying SNF to your own datasets then take a look at this package.

Functions from this package are used throughout the analyses and results generation.

### [`pypmi`](./pypmi)

This library is designed to automate aspects of working with the [Parkinson's Progression Markers Initiative](https://www.ppmi-info.org/) (PPMI) dataset.
The PPMI requires researchers sign a Data Usage Agreement (DUA) before accessing the data, which is stored on the USC LONI's Image and Data Archive.
I mostly wrote this because I was tired of having to manually navigate through the LONI user interface every time I wanted to download updated data, and figured no one would balk at having an automatic way to fetch and load datasets.
If you're working with demographic / behavioral / DAT scan / biospecimen data from the PPMI in any capacity I would be happy to have you use this package.

Functions from this package are used throughout the analyses and results generation.

## Supporting libraries

### [`netneurotools`](./netneurotools)

This package contains miscellaneous utilities and functions frequently used in the [Network Neuroscience Lab](netneurolab.github.io/).
Though its documentation is relatively spart, it has a lot of generally useful, re-usable functions that wouldn't quite fit in other mainstream Python modules.
It largely serves as a way for the graduate students and postdocs in the lab to share version-controlled code with one another (i.e., it avoids "hey let me e-mail you this script").

Functions from this package are used throughout preprocessing, analyses, and results generation.

### [`neurocombat`](./neurocombat)

This package implements the ComBat algorithm for correcting batch effects in data.
Published in support of analyses described in Fortin et al., 2017, *NeuroImage*, this fork of the [original repository](https://github.com/ncullen93/neuroCombat) was modified to make it a Python-installable package (i.e., none of the internals of the code were changed).

Functions from this package are only used in [`scripts/02_analysis/01_prepare_snf_data.py`](../scripts/02_analysis/01_prepare_snf_data.py)
