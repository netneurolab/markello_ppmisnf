# Similarity network fusion in Parkinson's disease

## "What's in this repository?"

This repository contains data, code, and results for the manuscript "Integrated morphometric, molecular, and clinical characterization of Parkinson's disease pathology." 
The study examines the application of similarity network fusion to Parkinson's disease using data from the [Parkinson's Progression Markers Initiative](https://www.ppmi-info.org) (PPMI).

We've tried to document the various aspects of this repository with a whole bunch of README files, so feel free to jump around and check things out.

## "Just let me run the things!"

Itching to just run the analyses? 
Get going with the following:

```bash
git clone --recurse-submodules https://github.com/netneurolab/markello_ppmisnf
cd markello_ppmisnf
pip install -r requirements.txt
export PYTHONPATH=$PYTHONPATH:$PWD/code/ppmi_snf:$PWD/code/netneurotools:$PWD/code/neurocombat:$PWD/code/pypmi:$PWD/code/snfpy
make all
```

If you don't want to deal with the hassle of creating a new Python environment, download the Singularity image that we used to run our analyses and run things in there:

```bash
git clone --recurse-submodules https://github.com/netneurolab/markello_ppmisnf
cd markello_ppmisnf
wget -O container/ppmi_snf.simg https://osf.io/h6jwx/download
bash container/run.sh
make all
```

## "I want to take things slow"

If you want a step-by-step through all the methods + analyses, take a look at out our [walkthrough](./walkthrough).

## "I have questions"

[Open an issue](https://github.com/netneurolab/markello_ppmisnf/issues) on this repository and someone will try and get back to you as soon as possible!
