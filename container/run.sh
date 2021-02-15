#!/usr/bin/env bash
#
# Description:
#
#     This script is used to generate, build, and run the Singularity container
#     that all analyses should take place in. This will ensure appropriate
#     versioning of software and a modicum of reproducibility across machines
#     and collaborators.
#
#     This script was initially written to be used on a Linux box running
#     Ubuntu 18.04. YMMV if you try it on any other system! At the very least,
#     you can extract the relevant codebits for creating the Singularity
#     container and run that (or use the pre-generated Singularity recipe in
#     the same directory as this script).
#
# Usage:
#
#     $ bash container/run.sh
#

curr_dir=$PWD && if [ ${curr_dir##*/} = "container" ]; then cd ..; fi
tag=ppmi_snf

# check for the $PPMI_USER and $PPMI_PASSWORD environmental variables
# if not set, prompt for user to set the dynamically so that we can instantiate
# them in the Singularity container
if [[ -z "${PPMI_USER}" || -z "${PPMI_PASSWORD}" ]]; then
  printf "One (or both) of the environmental variables \$PPMI_USER and "
  printf "\$PPMI_PASSWORD is\nnot set. This will prevent you from being able "
  printf "to programatically access PPMI\ndata. To prevent this from causing "
  printf "any issues it is recommended you download the\nPPMI data and place "
  printf "it directly into the data/raw/ppmi/behavior directory. For\n"
  printf "instructions on downloading PPMI data see https://github.com/"
  printf "rmarkello/pypmi.\nAlternatively, you can set these variables now;"
  printf "leave blank to continue\nunaltered if you followed the download "
  printf "instructions.\n\n"
  read -p "PPMI_USER: " ppmi_user
  read -s -p "PPMI_PASSWORD: " ppmi_password
  printf "\n"
  if [[ ! -z "${ppmi_user}" && -z "${PPMI_USER}" ]]; then
    PPMI_USER=${ppmi_user}
  fi
  if [[ ! -z "${ppmi_password}" && -z "${PPMI_PASSWORD}" ]]; then
    PPMI_PASSWORD=${ppmi_password}
  fi
fi

# use neurodocker (<3) to make a Singularity recipe and build the Singularity
# image. this should only happen if the image doesn't already exist
if [ ! -f container/${tag}.simg ]; then
  if [ ! -f container/license.txt ]; then touch container/license.txt; fi
  singularity --quiet exec docker://kaczmarj/neurodocker:0.4.0                \
    /usr/bin/neurodocker generate singularity                                 \
    --base ubuntu:16.04                                                       \
    --pkg-manager apt                                                         \
    --install                                                                 \
      git less nano libgl1-mesa-glx libglu1-mesa libxi6 libxkbcommon-dev      \
      libxcb-xkb-dev libxslt1-dev libgstreamer-plugins-base0.10-dev           \
    --ants                                                                    \
      version=2.2.0                                                           \
    --freesurfer                                                              \
      version=6.0.1                                                           \
      license_path=container/license.txt                                      \
      exclude_path=subjects/V1_average                                        \
    --copy ./environment.yml /opt/environment.yml                             \
    --miniconda                                                               \
      create_env=${tag}                                                       \
      yaml_file=/opt/environment.yml                                          \
    --add-to-entrypoint "source activate ${tag}"                              \
  > container/Singularity
  sudo singularity build container/${tag}.simg container/Singularity
fi

# export some environmental variables to be used by singularity
# first off is setting the PYTHONPATH in singularity; we want to be able to
# load the external libraries in the `code` directory
export SINGULARITYENV_PYTHONPATH=${PWD}/code
for pkg in ${PWD}/code/*; do
    export SINGULARITYENV_PYTHONPATH=${SINGULARITYENV_PYTHONPATH}:${pkg}
done
# next up is adding ANTS to the default path
export SINGULARITYENV_APPEND_PATH=/opt/ants-2.2.0
# setting the multithreading for blas/mkl so numpy doesn't explode
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
# setting the PPMI user/password variables for easy data fetching
export SINGULARITYENV_PPMI_PASSWORD=${PPMI_PASSWORD}
export SINGULARITYENV_PPMI_USER=${PPMI_USER}
# setting some display things
export SINGULARITYENV_DISPLAY=:1

# check to see if symlink to raw BIDS dataset exists
bids="$( readlink -f "data/raw/ppmi/bids" )"
if [ ! -d "${bids}" ]; then
    mkdir -p "data/raw/ppmi/bids"
fi
# check to see if symlink to BIDS derivatives dataset exists
antslct="$( readlink -f "data/derivative/antslct" )"
if [ ! -d "${antslct}" ]; then
    mkdir -p "data/derivative/antslct"
fi

# shell into container and bind paths, as appropriate
singularity shell --cleanenv --home ${PWD} \
                  -B ${antslct}:${antslct} \
                  -B ${bids}:${bids} \
                  -s "/neurodocker/startup.sh" \
                  container/${tag}.simg

# remove some of the things that are created by the container
for fn in .jupyter .local .ipython .cache .config matlab; do
  if [ -d ${PWD}/${fn} ]; then
    rm -fr ${PWD}/${fn}
  fi
done
if [ -f .bash_history ]; then
  rm -f .bash_history
fi

# unset all our environmental variables
unset SINGULARITYENV_APPEND_PATH
unset SINGULARITYENV_PYTHONPATH
unset SINGULARITYENV_MKL_NUM_THREADS
unset SINGULARITYENV_OPENBLAS_NUM_THREADS
unset SINGULARITYENV_PPMI_PASSWORD
unset SINGULARITYENV_PPMI_USER
unset SINGULARITYENV_DISPLAY
unset SINGULARITYENV_XDG_RUNTIME_DIR
