# Set-up and installation

## Required software

Reproducing these analyses require the following software (links go to installation instructions for each dependency):

- [Git](https://git-scm.com/), and
- [Python 3.6+](https://docs.conda.io/en/latest/miniconda.html)

If you intend on running the analyses from scratch (i.e., downloading the raw MRI images and re-processing them), then you will also need:

- [Docker](https://docs.docker.com/install/), and
- [Singularity](https://www.sylabs.io/guides/3.2/user-guide/quick_start.html)

## Getting the repository with `git`

First, you'll need a copy of all the data, code, and whatnot in the repository.
You can make a copy by running the following command:

```bash
git clone --recurse-submodules https://github.com/netneurolab/markello_ppmisnf
```

## Python dependencies

It is recommended that you create a new Python environment to install all the dependencies for the analyses.
If you'd like to simply install all the dependencies in your current Python environment you can do that, but no guarantees that things will work without issue!
(No guarantee things will work without issue even if you use environments and containers, but we're trying!)

### Using `conda` (recommended)

If you are using `conda` you can create a new environment and install all the required Python dependencies with the following command:

```bash
conda env create -f environment.yml
```

Alternatively you can add all the dependencies to your current environment with:

```bash
conda env update -f environment.yml
```

### Using `pip`

If you are using `pip` you can install all the dependencies into your current environment with:

```bash
pip install -r requireqments.txt
```

### Internal libraries

We've packaged some internal libraries with this repository. Once you've created your Python environment (or installed the dependencies as described above), you can install these internal libraries with the following commands:

```bash
pip install code/netneurotools code/neurocombat code/pypmi code/snfpy
export PYTHONPATH=$PYTHONPATH:$PWD/code
```

## Docker and Singularity

You only need these dependencies if you intend on reprocessing the neuroimaging data (the procedure for which is described in following steps).
All other analyses _should_ be deterministic and reproducible by simply using a versioned Python environment.
(_Should_ is doing a lot of work in that sentence.)

That said, we also provide a copy of the Singularity image which we used to run all the primary analyses, [accessible via OSF](https://osf.io/h6jwx/).
If you want to use that image to re-run the analyses you can consider installing Singularity.

---

[Click here to continue the walkthrough](./01_accessing_data.md)
