#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environment
conda create -n CAM \
    -c bioconda -c conda-forge -c pytorch \
    biasaway=3.3.0 \
    biopython=1.78 \
    click=7.1.2 click-option-group=0.5.1 \
    cudatoolkit=11.1.1 pytorch torchaudio torchvision \
    # dash=1.20.0 \
    jupyterlab=3.0.14 \
    logomaker=0.8 \
    matplotlib=3.4.1 \
    pandas=1.2.4 \
    # plotly=4.14.3 \
    pybedtools=0.8.2 \
    python=3.8.5 \
    # seaborn=0.11.1 \
    scikit-learn=0.24.1 \
    tqdm=4.60.0

# Pip install
conda activate deep-motif-discovery-pipeline
pip install genomepy
pip install snakemake
