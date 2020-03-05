# Processing the PPMI neuroimaging data

Congratulations, you've BIDS-ified your PPMI dataset!
Now that your PPMI dataset is in BIDS format you can get to work on processing it through any of the [BIDS Apps](https://doi.org/10.1371/journal.pcbi.1005209).
The world is your oyster!

We chose to process the PPMI neuroimaging data with a [slightly-modified](https://github.com/rmarkello/antslct/) version of the [ANTs longitudinal cortical thickness pipeline](https://www.biorxiv.org/content/10.1101/170209v2) [BIDS App](https://github.com/BIDS-Apps/antsCorticalThickness).
Note, this pipeline takes a **LONG** time and a **SIGNIFICANT** amount of computational resources to run, so we do NOT recommend you re-run this yourself (just use the [data](../data) we provide with this repository if you're trying to reproduce our analyses!).

## Re-running the ANTs pipeline

If you're *absolutely* determined to run things from scratch, you can re-run the modified ANTs pipeline with the following command:

```bash
N_CORES=8
singularity pull antslct.simg shub://rmarkello/antslct
for SUBJ in $( ls -d data/raw/ppmi/bids/sub-???? ); do
    SUBJ=$( basename ${SUBJ} | cut -d '-' -f2 )
    singularity run -B data/raw/ppmi/bids:/data                               \
                    -B data/derivative/ppmi/anstlct:/output                   \
                    antslct.simg                                              \
                    -s ${SUBJ} -o /output -c ${N_CORES} -m
done
```

(We strongly recommend you use some sort of HPC to parallelize this process across subjects rather than running it serially, as the code above does.
Otherwise, feel free to pick back up this walkthrough in a few years once your processing is done.)

This command will download a new Singularity image (`antslct.simg`) that operates similarly to a BIDS App and will run the ANTs pipeline on the neuroimaging data.

## The modified ANTs BIDS App

You may be asking: why didn't the authors just use the ANTs BIDS App as-is?

The ungenerous answer would be: the lead author wasn't thinking straight.

The slightly-nicer answer: we didn't *really* modify the ANTs pipeline at all. 
That is, we still run the ANTs longitudinal cortical thickness pipeline exactly as-is; we just deviated from the original *BIDS App implementation* because we wanted to do some extra pre- and post-processing specific to our questions / analyses.

We've tried to describe everything on the associated [GitHub repo](https://github.com/rmarkello/antslct/), but we reproduce relevant aspects of the documentation here:

> In many respects this code is similar to the ANTs BIDS-App, which is designed to automatically run the ANTs cortical thickness pipeline on a formatted BIDS datasets. However, this pipeline has a few notable differences from the BIDS App, including:
>
> 1. Automated use of additional data modalities as inputs to the pipeline, as available,
> 1. Use of the [MNI152 ICBM nonlinear asymmetric 2009c](http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009) atlas as the default group template,
> 1. Use of [labelled template files](https://drive.google.com/drive/folders/0B4SvObeEfaRyZGhlUlJOcmItTVU?usp=sharing) for better image segmentation using ANTs [joint label fusion](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3837555/),
> 1. Generation of additional output files, and
> 1. Generation of [fMRIPrep-style visual reports](https://fmriprep.readthedocs.io/en/stable/outputs.html#visual-reports) to aid in quality control
>
> This code will generate all of the "normal" ANTsLongitudinalCorticalThickness.sh outputs in the specified output directory, but will also create several additional derivatives:
>
> 1. `sub-XXXX_*_jacobian.nii.gz`: log jacobian determinant image of the concatenated non-linear warps (i.e., subject → subject-specific template AND subject-specific template → group template) for a single timepoint, in MNI space, where **positive values indicate relative subject atrophy** and **negative values indicate relative subject expansion** as compared to the MNI template
> 1. `sub-XXXX_*_invjacobian.nii.gz`: inverse log jacobian determinant image of the concatenated non-linear warps (i.e., subject → subject-specific template AND subject-specific template → group template) for a single timepoint, in MNI space, where **negative values indicate relative subject atrophy** and **positive values indicate relative subject expansion** as compared to the MNI template
> 1. `sub-XXXX_*_corticalthickness.nii.gz`: cortical thickness image for a single timepoint, in MNI space
> 1. `sub-XXXX_*_sscorticalthickness.nii.gz`: cortical thickness image for a single timepoint, in native space
> 1. `sub-XXXX_*_subjecttotemplatewarp.nii.gz`: combined warp image for a single timepoint, going from native space to MNI space
> 1. `sub-XXXX_*_templatetosubjectwarp.nii.gz`: combined warp image for a single timepoint, going from MNI space to native space
> 1. `sub-XXXX_*_brainsegmentation.nii.gz`: ANTs-style brain tissue segmentation for a single timepoint
> 1. `sub-XXX.html`: a visual report with images for checking the quality of the segmentation and registration for the generated subject-specific template and all individual timepoints

## Quality control

Two of the authors (RDM + CT) used the HTML reports generated from the slightly-modified ANTs pipeline to quality control the data.
We have provided their [quality control ratings](../data/derivative/antslct) in this repository.
(Note: we were **very conservative** in our ratings, so you're gonna see a lot of rejects.)
These ratings were used to determine which subjects should be carried forward throughout the analyses.

It's possible, due to idiosyncracies in non-linear registration, that re-running the pipeline as described above will change some of the generated neuroimaging derivatives, and subjects who we "failed" or "passed" may need to be "passed" or "failed" when run anew.
Again, we provide the derived neuroimaging data so you don't have to make this call yourself!

## Parcellating the data

Now that you've processed and quality-controlled the data, you need to parcellate it.
You can do that from the main directory with the following command:

```bash
make preprocess
```

This command will run the two files located in the [`scripts/01_preprocess`](../scripts/01_preprocess) directory, creating outputs in the [`data/derivative/parcellated`](../data/derivative/parcellated) directory.
(These outputs should already exist if you cloned this repository from GitHub!)

---

[Click here to continue the walkthrough](./04_snf_analyses.md)
