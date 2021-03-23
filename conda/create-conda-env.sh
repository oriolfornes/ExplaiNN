#!/bin/bash

# i.e. enable conda (de)activate
eval "$(conda shell.bash hook)"

# Create conda environment
conda create -n deep-motif-discovery-pipeline \
    -c bioconda -c conda-forge -c pytorch \
    biasaway=3.3.0 biopython=1.78 click=7.1.2 click-option-group=0.5.1 \
    dash=1.19.0 ignite=0.4.4 jupyterlab=3.0.11 logomaker=0.8 matplotlib=3.3.4 \
    meme=4.11.2 pandas=1.2.3 plotly=4.14.3 pybedtools=0.8.2 pynvml=8.0.4 \
    python=3.8.5 pytorch=1.8.0 seaborn=0.11.1 scikit-learn=0.24.1 \
    tensorboard=2.4.1 torchvision=0.9.0 tqdm=4.59.0

# Pip install
conda activate deep-motif-discovery-pipeline
pip install genomepy
pip install snakemake
