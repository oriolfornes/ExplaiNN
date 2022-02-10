#!/usr/bin/env python

from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import click
from functools import partial
import gc
import gzip
import matplotlib.pyplot as plt
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import pickle
import re
import time
import torch
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Local imports
from architectures import get_metrics
from jaspar import get_figure, reformat_motif
from sequence import rc_many
from test import _get_acts_outs_preds, __get_fwd_rev, _load_model
from train import _get_seqs_labels_ids, _get_data_loader
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
    "training_file",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-b", "--batch-size",
    help="Batch size.",
    type=int,
    default=100,
    show_default=True,
)
@click.option(
    "-c", "--cpu-threads",
    help="Number of CPU threads to use.",
    type=int,
    default=1,
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
@click.option(
    "-t", "--time-me",
    help="Return the program's running execution time in seconds.",
    is_flag=True,
)

def main(**args):

    # Start execution
    start_time = time.time()

    ##############
    # Load Data  #
    ##############

    # Initialize
    data_loaders = {}

    # Get data
    seqs, labels, ids = _get_seqs_labels_ids(args["training_file"],
        args["debugging"], args["rev_complement"])

    # Get DataLoader
    data_loader = _get_data_loader(seqs, labels, args["batch_size"])

    # Get sequences
    s = _get_sequences(args["training_file"])
    seqs_str = np.concatenate((s, rc_many(s)))

    # Load model
    exp_model = _load_model(args["model_file"])

    ##############
    # Interpret  #
    ############## 

    # Initialize
    if np.unique(labels[:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"

    # Create output dirs
    if not os.path.isdir(args["output_dir"]):
        os.makedirs(args["output_dir"])
    for subdir in ["sites", "motifs", "logos"]:
        if not os.path.isdir(os.path.join(args["output_dir"], subdir)):
            os.makedirs(os.path.join(args["output_dir"], subdir))

    # Get output layer weights
    weights = exp_model.final.weight.detach().cpu().numpy()
    N_motifs = weights.shape[1]
    tsv_file = os.path.join(args["output_dir"], "output-layer-weights.tsv")
    if not os.path.exists(tsv_file):
        data = []
        for i, weight in enumerate(weights.T):
            data.append([f"filter{i}"] + weight.tolist())
        column_names = ["filter"] + list(range(weights.shape[0]))
        df = pd.DataFrame(data, columns=column_names)
        df.to_csv(tsv_file, sep="\t", index=False)

    # Get filter weights
    filter_weights = exp_model.linears[0].weight
    pickle_file = os.path.join(args["output_dir"], "filter-weights.pickle")
    if not os.path.exists(pickle_file):
        filter2weights = {}
        for i in range(filter_weights.data.shape[0]):
            filter_w = filter_weights.data[i, :, :].detach().cpu().numpy()
            filter2weights.setdefault((args["name"], f"filter{i}"), filter_w.T)
        with open(pickle_file, "wb") as fh:
            pickle.dump(filter2weights, fh, protocol=pickle.HIGHEST_PROTOCOL)

    # Get activations, outputs, and predictions
    acts, outs, preds = _get_acts_outs_preds(exp_model, data_loader)

    # Get the indices of well predicted sequences
    idxs = _get_idxs_of_well_predicted_sequences(preds, labels, input_data,
        args["ignore_negatives"], args["rev_complement"])

    # For each filter, get the activation threshold (i.e. ≥50%)
    thresholds = _get_act_thresholds(acts, idxs, args["rev_complement"])

    # Get filter importances
    kwargs = {"bar_format": bar_format, "total": N_motifs}
    tsv_file = os.path.join(args["output_dir"], "filter-importances.tsv")
    if not os.path.exists(tsv_file):
        data = []
        imps = np.array([np.multiply(outs, weights[i, :]) \
            for i in range(weights.shape[0])])
        for i in tqdm(range(N_motifs), **kwargs):
            filt_imps = _filter_filter_importances(imps[:, :, i], idxs,
                acts[:, i, :], thresholds[i])
            data.extend([[f"filter{i}"] + k.tolist() for k in filt_imps.T])
        column_names = ["filter"] + list(range(imps.shape[0]))
        df = pd.DataFrame(data, columns=column_names)
        df.to_csv(f"{tsv_file}.gz", sep="\t", index=False, compression="gzip")
        df = df.groupby(["filter"]).median().sort_values([column_names[-1]],
            ascending=False)
        df.reset_index(inplace=True)
        df.to_csv(tsv_file, sep="\t", index=False)

    # For each filter...
    for i in tqdm(range(N_motifs), **kwargs):

        # Get sites
        sites_file = os.path.join(args["output_dir"], "sites",
            f"filter{i}.fa.gz")
        if not os.path.exists(sites_file):
            fh = get_file_handle(sites_file, "wt")
            _get_sites(fh, idxs, ids, seqs_str, acts[:, i, :], thresholds[i],
                exp_model._options["kernel_size"], args["rev_complement"])
            fh.close()
 
        # Get motif
        motif_file = os.path.join(args["output_dir"], "motifs",
            f"filter{i}.jaspar")
        if not os.path.exists(motif_file):
            fh = get_file_handle(motif_file, "wt")
            try:
                motif = _sites_to_motif(sites_file)
                motif.matrix_id = f"filter{i}"
                motif.name = args["name"]
                fh.write(format(motif, "jaspar"))
            except:
                # i.e. no sites
                pass
            fh.close()

    # Free memory
    del seqs
    del labels
    del ids
    del data_loader
    del seqs_str
    del exp_model
    del weights
    del filter_weights
    del acts
    del outs
    del preds
    del idxs
    del thresholds
    gc.collect()

    # Get motifs in MEME format
    motif_files = []
    jaspar_motifs = []
    for i in range(N_motifs):
        motif_file = os.path.join(args["output_dir"], "motifs",
            f"filter{i}.jaspar")
        motif_files.append(motif_file)
        try:
            with open(motif_file) as fh:
                jaspar_motifs.append(motifs.read(fh, "jaspar"))
        except:
            # i.e. no motif
            pass
    meme_file = os.path.join(args["output_dir"], "motifs", "filters.meme")
    if not os.path.exists(meme_file):
        reformat_motif(jaspar_motifs, "meme", meme_file)

    # Generate logos
    pool = Pool(args["cpu_threads"])
    p = partial(_generate_logo, output_dir=args["output_dir"])
    for _ in tqdm(pool.imap(p, motif_files), **kwargs):
        pass

    # Finish execution
    if args["time_me"]:
        seconds = format(time.time() - start_time, ".2f")
        txt_file = os.path.join(args["output_dir"],
            f"time-{os.path.basename(__file__)}.txt")
        fh = get_file_handle(txt_file, "wt")
        fh.write(f"{seconds} seconds")
        fh.close()

def _get_sequences(tsv_file):

    # Get sequences
    df = pd.read_csv(tsv_file, sep="\t", header=None)

    return(df[1].values)

def _get_idxs_of_well_predicted_sequences(preds, labels, input_data,
                                          ignore_negatives=False,
                                          rev_complement=False):

    # Initialize
    n = .05

    if rev_complement:
        fwd = __get_fwd_rev(preds, "fwd")
        rev = __get_fwd_rev(preds, "rev")
        p = np.empty(fwd.shape)
        ys = __get_fwd_rev(labels, "fwd")
        # Avg. predictions from forward and reverse strands
        for i in range(p.shape[1]):
            p[:, i] = np.mean([fwd[:, i], rev[:, i]], axis=0)
            if input_data == "binary":
                p[:, i] = torch.sigmoid(torch.from_numpy(p[:, i])).numpy()
    else:
        if input_data == "binary":
            p = torch.sigmoid(torch.from_numpy(preds)).numpy()
        else:
            p = preds
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
    else:
        idxs_ys = np.argsort(-ys.flatten())[:int(max(ys.shape)*n)]
        idxs_p = np.argsort(-p.flatten())[:int(max(p.shape)*n)]
        idxs = np.intersect1d(idxs_ys, idxs_p)

    return(idxs)

def _get_act_thresholds(acts, idxs, rev_complement=False):

    # Get activations of well predicted sequences
    if rev_complement:
        fwd = __get_fwd_rev(acts, "fwd")        
        rev = __get_fwd_rev(acts, "rev")
        acts_idxs = np.concatenate((fwd[idxs, :, :], rev[idxs, :, :]))
    else:
        acts_idxs = acts[idxs, :, :]

    return(0.5 * np.amax(acts_idxs, axis=(0, 2)))

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

def _filter_filter_importances(imps, idxs, acts, threshold):

    # Initialize
    idxs = np.intersect1d(idxs, np.where(acts > threshold)[0])

    return(imps[:, idxs])

def _generate_logo(motif_file, output_dir="./"):

    # Initialize
    m = re.search("^filter(\d+).jaspar$", os.path.split(motif_file)[1])

    for reverse_complement in [False, True]:
        if reverse_complement:
            logo_file = f"filter{m.group(1)}.rev.png"
        else:
            logo_file = f"filter{m.group(1)}.fwd.png"
        logo_file = os.path.join(output_dir, "logos", logo_file)
        if not os.path.exists(logo_file):
            try:
                fig = get_figure(motif_file, reverse_complement)
                fig.savefig(logo_file, bbox_inches="tight", pad_inches=0)
            except:
                # i.e. no motif
                fh = get_file_handle(logo_file, "wt")
                fh.close()

if __name__ == "__main__":
    main()