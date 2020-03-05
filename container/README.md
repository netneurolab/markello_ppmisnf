# container

This directory contains files that were used to ensure analyses were run within a (semi-)reproducible [Singularity](https://sylabs.io/docs/) container.
You can download the Singularity image used for all reported analyses from OSF [here](https://osf.io/h6jwx/).
Placing that image into this directory should allow you to use the `run.sh` without needing to re-build a Singularity container!

- [**`run.sh`**](run.sh): This is a helper script designed to make instantiating the Singularity container a bit easier.
  Running this will drop users into a bash shell within the Singularity container, from which all scripts / files should be accessible and executable.
- [**`Singularity`**](Singularity): The Singularity recipe generated from `run.sh` and used to build the Singularity image in which all analyses were executed.

To reproduce all analyses the following code *should* (famous last words) suffice:

```bash
bash container/run.sh  # will drop you into a Singularity container
make all               # will run everything and re-generate the manuscript PDF
```
