#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environment
conda create -n CAM \
    -c pytorch -c conda-forge -c bioconda \
    biasaway=3.3.0 \
    biopython=1.78 \
    click=7.1.2 click-option-group=0.5.1 \
    cudatoolkit=11.0.3 pytorch torchaudio torchvision \
    genomepy=0.9.3 \
    jupyterlab=3.0.14 \
    logomaker=0.8 \
    matplotlib=3.2.2 \
    pandas=1.2.3 \
    pybedtools=0.8.2 \
    python=3.8.3 \
    scikit-learn=0.23.1 \
    tqdm=4.47.0 \
    umap-learn

# Pip install
#conda activate deep-motif-discovery-pipeline
#pip install genomepy
#pip install snakemake
