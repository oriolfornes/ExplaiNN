#!/usr/bin/env python

from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import click
import gzip
from io import StringIO
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import torch
from tqdm import tqdm
import sys
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Local imports
from architectures import get_metrics
from jaspar import get_figure, reformat_motif
from sequence import rc_many
from test import _get_acts_labs_outs_preds, __get_fwd_rev, _load_model
from train import _get_Xs_ys_ids, _get_data_loader
from utils import get_file_handle

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "model_file",
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "tsv_files",
    type=click.Path(exists=True, resolve_path=True),
    nargs=-1,
)
@click.option(
    "-b", "--batch-size",
    help="Batch size.",
    type=int,
    default=2**6,
    show_default=True,
)
@click.option(
    "-d", "--debugging",
    help="Debugging mode.",
    is_flag=True,
)
@click.option(
    "-i", "--ignore-negatives",
    help="Derive motifs from positive sequences only.",
    is_flag=True,
)
@click.option(
    "-n", "--name",
    help="Name for `Bio.motifs` class.",
    type=str,
    required=True,
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)
@click.option(
    "-r", "--rev-complement",
    help="Reverse complement sequences.",
    is_flag=True,
)

def main(**args):

    ##############
    # Load Data  #
    ##############

    # Initialize
    Xs = {}
    ys = {}
    ids = {}
    sequences = {}
    data_loaders = {}

    # Get data / sequences
    for i in range(len(args["tsv_files"])):
        tsv_file = args["tsv_files"][i]
        X, y, id = _get_Xs_ys_ids(tsv_file, args["debugging"],
            args["rev_complement"])
        Xs.setdefault(i, X)
        ys.setdefault(i, y)
        ids.setdefault(i, id)
        s = _get_sequences(tsv_file)
        sequences.setdefault(i, np.concatenate((s, rc_many(s))))

    # Get DataLoaders
    for i in Xs:
        data_loader = _get_data_loader(Xs[i], ys[i], args["batch_size"])
        data_loaders.setdefault(i, data_loader)

    # Load model
    model = _load_model(args["model_file"])

    ##############
    # Interpret  #
    ############## 

    # Initialize
    activations = {}
    idxs = {}
    labels = {}
    outputs = {}
    predictions = {}
    filter_importances = {}
    if np.unique(ys[0][:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"

    # Create output dirs
    if not os.path.isdir(args["output_dir"]):
        os.makedirs(args["output_dir"])
    for subdir in ["motifs-new", "logos-new"]:
        if not os.path.isdir(os.path.join(args["output_dir"], subdir)):
            os.makedirs(os.path.join(args["output_dir"], subdir))

    # Get weights
    weights = model.final.weight.detach().cpu().numpy()
    weights_file = os.path.join(args["output_dir"], "weights.tsv")
    if not os.path.exists(weights_file):
        data = []
        for i, weight in enumerate(weights.T):
            data.append([f"filter{i}"] + weight.tolist())
        column_names = ["filter"] + list(range(weights.shape[0]))
        df = pd.DataFrame(data, columns=column_names)
        df.to_csv(weights_file, sep="\t", index=False)

    # # Get activations, labels, outputs, predictions, and filter importances
    # for i in sequences:
    #     data_loader = data_loaders[i]
    #     activation, label, output, prediction = \
    #         _get_acts_labs_outs_preds(model, data_loader)
    #     activations.setdefault(i, activation)
    #     labels.setdefault(i, label)
    #     outputs.setdefault(i, output)
    #     predictions.setdefault(i, prediction)
    #     filter_importances.setdefault(i,
    #         np.array([np.multiply(output, weights[j, :]) \
    #             for j in range(weights.shape[0])]))

    # # Get the indices of well predicted sequences
    # for i in sequences:
    #     idx = _get_idxs_of_well_predicted_sequences(predictions[i],
    #         labels[i], input_data, args["ignore_negatives"],
    #         args["rev_complement"])
    #     idxs.setdefault(i, idx)

    # # For each filter, get the activation threshold (i.e. ≥50%)
    # thresholds = _get_activation_thresholds(activations, idxs,
    #     args["rev_complement"])

    # Get filters & logos
    kernels = model.linears[0].weight
    kwargs = {"bar_format": bar_format, "total": kernels.shape[0]}
    for k in tqdm(range(kernels.shape[0]), **kwargs):

        k=61
        kernel = kernels[k, :, :].detach().cpu().numpy()
        kernel = np.exp(kernel * 100.0).T
        prob = kernel / np.nansum(kernel, axis=1)[:, None]
        df = pd.DataFrame(data=prob, columns=list("ACGT"))
        print(df)
        exit(0)
        # df = pd.DataFrame(data=prob, columns=list("ACGT")).multiply(100)
        # df["Sum"] = df["A"] + df["C"] + df["G"] + df["T"]
        sio = StringIO(df.to_string(header=False, index=False))
        m = motifs.read(sio, "pfm-four-columns")
        sio.close()
        m.matrix_id = f"filter{k}"
        m.name = args["name"]
        fh = get_file_handle(motif_file, "wt")
        fh.write(format(m, "jaspar"))
        fh.close()


        # Get kernel
        motif_file = os.path.join(args["output_dir"], "motifs-new",
            f"filter{k}.jaspar")
        if not os.path.exists(motif_file):
            kernel = kernels[k, :, :].detach().cpu().numpy()
            kernel = np.exp(kernel * 100.0).T
            prob = kernel / np.nansum(kernel, axis=1)[:, None]
            df = pd.DataFrame(data=prob, columns=list("ACGT")).multiply(100)
            # df["Sum"] = df["A"] + df["C"] + df["G"] + df["T"]
            sio = StringIO(df.to_string(header=False, index=False))
            m = motifs.read(sio, "pfm-four-columns")
            sio.close()
            m.matrix_id = f"filter{k}"
            m.name = args["name"]
            fh = get_file_handle(motif_file, "wt")
            fh.write(format(m, "jaspar"))
            fh.close()

        # Get logos
        for reverse_complement in [False, True]:
            if reverse_complement:
                logo_file = os.path.join(args["output_dir"], "logos-new",
                    f"filter{k}.rev.png")
            else:
                logo_file = os.path.join(args["output_dir"], "logos-new",
                    f"filter{k}.fwd.png")
            if not os.path.exists(logo_file):
                fig = get_figure(motif_file, reverse_complement)
                fig.savefig(logo_file, bbox_inches="tight", pad_inches=0)

    # Get motifs in MEME format
    meme_file = os.path.join(args["output_dir"], "motifs-new", "filters.meme")
    if not os.path.exists(meme_file):
        jaspar_motifs = []
        for i in range(model._options["cnn_units"]):
            motif_file = os.path.join(args["output_dir"], "motifs-new",
                f"filter{i}.jaspar")
            with open(motif_file) as fh:
                jaspar_motifs.append(motifs.read(fh, "jaspar"))
        reformat_motif(jaspar_motifs, "meme", meme_file)

    # # Get filter importances
    # tsv_file = os.path.join(args["output_dir"], "filter-importances.tsv")
    # if not os.path.exists(tsv_file):
    #     data = []
    #     for i in tqdm(range(thresholds.size), **kwargs):
    #         for j in outputs:
    #             arr = _filter_filter_importances(filter_importances[j][:, :, i],
    #                 idxs[j], activations[j][:, i, :], thresholds[i])
    #             data.extend([[f"filter{i}"] + k.tolist() for k in arr.T])
    #     column_names = ["filter"] + list(range(filter_importances[0].shape[0]))
    #     df = pd.DataFrame(data, columns=column_names)
    #     df.to_csv(f"{tsv_file}.gz", sep="\t", index=False, compression="gzip")
    #     df = df.groupby(["filter"]).median().sort_values([column_names[-1]],
    #         ascending=False)
    #     df.reset_index(inplace=True)
    #     df.to_csv(tsv_file, sep="\t", index=False)

def _get_sequences(tsv_file):

    # Get sequences
    df = pd.read_csv(tsv_file, sep="\t", header=None)

    return(df[1].values)

def _get_idxs_of_well_predicted_sequences(predictions, labels, input_data,
                                          ignore_negatives=False,
                                          rev_complement=False):
    
    if rev_complement:
        fwd = __get_fwd_rev(predictions, "fwd")
        rev = __get_fwd_rev(predictions, "rev")
        p = np.empty(fwd.shape)
        ys = __get_fwd_rev(labels, "fwd")
        # Average predictions from forward and reverse strands
        for i in range(p.shape[1]):
            p[:, i] = np.mean([fwd[:, i], rev[:, i]], axis=0)
            if input_data == "binary":
                p[:, i] = torch.sigmoid(torch.from_numpy(p[:, i])).numpy()
    else:
        if input_data == "binary":
            p = torch.sigmoid(torch.from_numpy(predictions)).numpy()
        else:
            p = predictions
        ys = labels

    # For binary data such as ChIP-seq, well-predicted sequences are bound
    # regions (i.e. 1) with a score greater than 0.5, and unbound regions
    # (i.e. 0) with a score lower or equal than 0.5.
    if input_data == "binary":
        c = Counter(np.where(ys == (p > .5).astype(int))[0])
        idxs = np.array([k for k, v in c.items() if v == ys.shape[1]])
        if ignore_negatives:
            idxs = np.intersect1d(idxs, np.where(np.sum(ys, axis=1) != 0)[0])

    # For non-binary data such as PBM, well-predicted sequences are defined 
    # as the top `n` percentile probes with the highest signal intensities
    # and a score in the top `n` percentile.
    # else:
    #     l = np.argsort(-labels.flatten())[:int(max(labels.shape)*top)]
    #     o = np.argsort(-outputs.flatten())[:int(max(outputs.shape)*top)]
    #     idxs = np.intersect1d(l, o)

    return(idxs)

def _get_activation_thresholds(activations, idxs, rev_complement=False):

    # Initialize
    x = 0
    y = activations[0].shape[1]
    z = activations[0].shape[2]
    activations_arr = np.empty((x, y, z))

    # Get activations of well predicted sequences
    for i in activations:
        if rev_complement:
            fwd = __get_fwd_rev(activations[i], "fwd")
            rev = __get_fwd_rev(activations[i], "rev")
            activations_arr = np.concatenate((activations_arr,
                fwd[idxs[i], :, :], rev[idxs[i], :, :]))
        else:
            activations_arr = np.concatenate((activations_arr,
                activations[i][idxs[i], :, :]))

    return(0.5 * np.amax(activations_arr, axis=(0, 2)))

def _get_sites(fh, idxs, ids, sequences, activations, threshold,
               kernel_size=19, rev_complement=False):
    """
    For each filter and each sequence, get sites reaching at least ½ of the
    max. activation value for that filter.
    """

    # Initialize
    count = 0

    # For each strand...
    for strand in ["+", "-"]:

        if rev_complement:
            ids_arr = __get_fwd_rev(ids, strand)
            seqs_arr = __get_fwd_rev(sequences, strand)
            activations_arr = __get_fwd_rev(activations, strand)
        else:
            ids_arr = ids
            seqs_arr = sequences
            activations_arr = activations

        # For each sequence...
        for i in idxs:

            # Get start positions of sequence sites activating this filter
            starts = np.where(activations_arr[i, :] > threshold)[0]

            # For each start...
            for j in starts:

                # Get site
                start = j
                end = j + kernel_size
                seq = Seq(seqs_arr[i][start:end])
                seq_id = f"{ids_arr[i]}_{strand}_from={start}_to={end}"
                record = SeqRecord(seq, id=seq_id, name="", description="")
                fh.write(record.format("fasta"))

                # If count reaches 1M, this filter is way too ubiquitous!!!
                count += 1
                if count == 1e6:
                    return

        # Skip
        if not rev_complement:
            return

def _sites_to_motif(sites_file):
    """
    From github.com/biopython/biopython/blob/master/Bio/motifs/__init__.py
    Read the motif from JASPAR .sites file.
    """

    # Initialize
    alphabet = "ACGTN"
    instances = []
    pfm = {}

    with gzip.open(sites_file, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                break
            # line contains the header ">...."
            # now read the actual sequence
            line = next(fh)
            instance = ""
            for c in line.strip().upper():
                if c == c.upper():
                    instance += c
            instance = Seq(instance)
            instances.append(instance)

    instances = motifs.Instances(instances, alphabet)
    motif = motifs.Motif(alphabet=alphabet, instances=instances)

    for nt in alphabet[:-1]:
        pfm.setdefault(nt, motif.counts[nt])

    return(motifs.Motif(counts=pfm))

def _filter_filter_importances(filter_importances, idxs, activations,
                               threshold):

    # Initialize
    idxs = np.intersect1d(idxs, np.where(activations > threshold)[0])

    return(filter_importances[:, idxs])

if __name__ == "__main__":
    main()