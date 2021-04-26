#!/usr/bin/env python

from re import X
from Bio import SeqIO
from Bio import motifs
import click
from click_option_group import optgroup
import copy
import gc
import gzip
from ignite.contrib.handlers.tqdm_logger import ProgressBar
from ignite.contrib.metrics import AveragePrecision, GpuInfo, ROC_AUC
from ignite.engine import (Events, create_supervised_trainer,
    create_supervised_evaluator)
from ignite.metrics import Loss
from io import StringIO
from itertools import cycle, islice
import json
import logging
import numpy as np
import os
from scipy.stats import pearsonr, spearmanr
from sklearn.model_selection import PredefinedSplit
import sys
import torch
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm

# Local imports
from utils.architectures import CAM as Model, NonStrandSpecific
from utils.jaspar import get_figure, reformat_jaspar_motif
from utils.selene import load_model_from_state_dict
from utils.sequence import one_hot_encode, one_hot_decode, \
                           rev_comp_one_hot_encoding

# CUDA
device = "cpu"
if torch.cuda.is_available():
    device = "cuda"

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "model_file",
    type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "training_file",
    type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-b", "--batch-size",
    help="Batch size.",
    type=int,
    default=2**6,
    show_default=True,
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)
@click.option(
    "-t", "--threads",
    help="Number of threads to use.",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "-r", "--rev-complement",
    help="Reverse complement training sequences.",
    is_flag=1
)
@optgroup.group("CAM")
@optgroup.option(
    "--cnn-units",
    help="Number of CNN units.",
    type=int,
    default=5,
    show_default=True,
)
@optgroup.option(
    "--data",
    help="Input data type.",
    type=click.Choice(["binary", "continuous"]),
    required=True,
)
@optgroup.option(
    "--motif-length",
    help="Length for motifs.",
    type=int,
    default=19,
    show_default=True,
)

def main(**params):

    # Create output dirs
    if not os.path.isdir(params["output_dir"]):
        os.makedirs(params["output_dir"])
    # for subdir in ["profiles", "sites", "logos"]:
    for subdir in ["profiles", "logos"]:
        if not os.path.isdir(os.path.join(params["output_dir"], subdir)):
            os.makedirs(os.path.join(params["output_dir"], subdir))

    # Get Xs/ys
    Xs, ys = _get_Xs_ys(params["training_file"])

    # Get DataLoader
    data_loader = _get_data_loader(list(Xs), list(ys), params["batch_size"],
        params["rev_complement"], params["threads"])

    # Load model
    selene_dict = torch.load(params["model_file"], map_location=lambda storage,
        location: storage)
    model = Model(params["cnn_units"], params["motif_length"],
        max(Xs[0].shape))
    if selene_dict["arch"] == "NonStrandSpecific":
        model = NonStrandSpecific(model)
    model = load_model_from_state_dict(selene_dict["state_dict"], model)
    if selene_dict["arch"] == "NonStrandSpecific":
        model = model.model
    model.to(device)

    # Free up memory
    _release_memory(Xs)
    _release_memory(ys)

    # Get profiles and sites
    # profiles, sites = _get_profiles_and_sites(model, data_loader)
    profiles = _get_profiles(model, data_loader)

    # Get weights
    weights = model.final.weight.detach().cpu().numpy().flatten().tolist()

    # Zip
    # zipped_list = zip(profiles, sites, weights)
    zipped_list = zip(profiles, weights)
    for f, z in enumerate(sorted(zipped_list, key=lambda x: x[-1], reverse=True)):
        jaspar_file = os.path.join(params["output_dir"], "profiles",
            "filter=%s.jaspar" % str(f+1))
        z[0].matrix_id = "filter=%s" % str(f+1)
        z[0].name = "CAM;l=%s" % str(params["motif_length"])
        with open(jaspar_file, "w") as o:
            o.write(format(z[0], "jaspar"))
        meme_file = os.path.join(params["output_dir"], "profiles",
            "filter=%s.meme" % str(f+1))
        reformat_jaspar_motif(jaspar_file, "meme", meme_file)
        png_file = os.path.join(params["output_dir"], "logos",
            "filter=%s.png" % str(f+1))
        fig = get_figure(jaspar_file)
        fig.savefig(png_file, bbox_inches="tight", pad_inches=0)
        # sites_file = os.path.join(params["output_dir"], "sites",
        #     "filter=%s.txt" % str(f+1))
        # with open(sites_file, "w") as o:
        #     o.write("\n".join(z[1]))
    with open(os.path.join(params["output_dir"], "weights.txt"), "w") as o:
        for f, w in enumerate(weights):
            o.write("filter=%s\t%s\n" % (str(f+1), str(w)))

# def _get_profiles_and_sites(model, data_loader):
def _get_profiles(model, data_loader):

    # Initialize
    outputs = []
    labels = []
    sequences = []
    profiles = []

    with torch.no_grad():
        for x, label in tqdm(data_loader, total=len(data_loader)):
            for encoded_seq in x:
                sequence = "N" * model._motif_length
                sequence += one_hot_decode(encoded_seq.numpy())
                sequence += "N" * model._motif_length
                sequences.append(one_hot_encode(sequence))
            x = x.to(device)
            out = model(x)
            outputs.extend(out.detach().cpu().numpy())
            labels.extend(label.numpy())

    # Get activations
    activations = __get_activations(model, data_loader)

    # Free memory
    _release_memory(data_loader)

    # Get Position Frequency Matrices (PFMs)
    outputs = np.array(outputs)
    labels = np.array(labels)
    # pfms, sites = __get_pfms_and_sites(sequences, activations,
    pfms = __get_pfms(sequences, activations, model._motif_length)

    # Get motif
    for j in range(len(pfms)):
        handle = StringIO("\n".join(["\t".join(map(str, i)) for i in pfms[j]]))
        profiles.append(motifs.read(handle, "pfm-four-rows"))

    # return(profiles, sites)
    return(profiles)

def __get_activations(model, data_loader):

    # Initialize
    activations = torch.tensor([], dtype=torch.float32)

    with torch.no_grad():
        for x, _ in tqdm(data_loader, total=len(data_loader)):
            x = x.to(device)
            x = x.repeat(1, model._cnn_units, 1)
            activation = model.linears[:3](x)
            activations = torch.cat([activations, activation.cpu()])

    return(activations.numpy())

# def __get_pfms_and_sites(sequences, activations, motif_length=19):
def __get_pfms(sequences, activations, motif_length=19):

    """
    For each filter, build a Position Frequency Matrix (PFM) from all sites
    reaching at least ½ the maximum activation value for that filter across all input sequences
    
     the
    activations and the original sequences, and keep the sites used to
    derive such matrix.

    params :
        actvations (np.array) : (N*N_filters*L) array containing the ourput for each filter and selected sequence of the test set
        sequnces (np.array) : (N*4*200) selected sequences (ACGT)
        y (np.array) : (N*T) original target of the selected sequnces
        output_file_path (str) : path to directory to store the resulting pwm meme file
    """

    # Initialize
    n_filters = activations.shape[1]
    pfms = np.zeros((n_filters, 4, motif_length))
    # sites = [[] for _ in range(n_filters)]

    # Find the threshold value for activations (i.e. 50%)
    activation_thresholds = 0.5*np.amax(activations, axis=(0, 2))

    # For each filter...
    for i in range(n_filters):

        activated_sequences_list = []

        # For each sequence...
        for j in range(len(sequences)):

            # Get indices of sequences that activate the filter
            idx = np.where(activations[j,i,:] > activation_thresholds[i])

            for ix in idx[0]:

                s = sequences[j][:,ix:ix+motif_length]
                # activated_sequences_list.append(s)

                # Build PFM
                pfms[i] += np.sum(s)
                print(pfms[i])

        # # If activated sequences...
        # if activated_sequences_list:

        #     # Convert activated sequences to array
        #     activated_sequences_arr = np.stack(activated_sequences_list)

        #     # Build PFM
        #     pfms[i] = np.sum(activated_sequences_arr, axis=0)

            # # Save sites that activated the filter
            # for s in activated_sequences_list:
            #     sites[i].append(one_hot_decode(s))

    # return(pfms, sites)
    return(pfms)

def _get_Xs_ys(fasta_file):

    # Initialize
    Xs = []
    ys = []

    # Xs / ys
    handle = _get_handle(fasta_file)
    for record in SeqIO.parse(handle, "fasta"):
        _, y_list = record.description.split()
        Xs.append(one_hot_encode(str(record.seq).upper()))
        ys.append([float(y) for y in y_list.split(";")])

    return(np.array(Xs), np.array(ys))

def _get_data_loader(Xs_train, ys_train, batch_size=2**6,
    reverse_complement=False, threads=1):

    # Reverse complement
    if reverse_complement:
        n = len(Xs_train)
        for i in range(n):
            encoded_seq = rev_comp_one_hot_encoding(Xs_train[i])
            Xs_train.append(encoded_seq)
            ys_train.append(ys_train[i])

    # TensorDatasets
    train_set = TensorDataset(torch.Tensor(Xs_train), torch.Tensor(ys_train))

    # DataLoaders
    kwargs = dict(batch_size=batch_size, num_workers=threads)
    train_loader = DataLoader(train_set, **kwargs)

    return(train_loader)

def _get_handle(file_name):
    if file_name.endswith("gz"):
        handle = gzip.open(file_name, "rt")
    else:
        handle = open(file_name, "rt")
    return(handle)

def _release_memory(my_object):
   del my_object
   gc.collect()

if __name__ == "__main__":
    main()