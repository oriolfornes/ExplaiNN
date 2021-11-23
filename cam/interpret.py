#!/usr/bin/env python

from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import click
import gzip
import numpy as np
import os
import pandas as pd
import torch
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Local imports
from architectures import CAM, get_metrics
from jaspar import get_figure, reformat_motif
from sequence import rc_many
from train import _get_Xs_ys_ids, _get_data_loader

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
@click.argument(
    "validation_file",
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
    "-d", "--debugging",
    help="Debugging mode.",
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
# @click.option(
#     "-r", "--rev-complement",
#     help="Reverse complement training sequences.",
#     is_flag=True,
# )

def main(**args):

    ##############
    # Load Data  #
    ##############

    # Get data
    Xs_train, ys_train, seq_ids_train = _get_Xs_ys_ids(args["training_file"],
        args["debugging"], True)
    Xs_val, ys_val, seq_ids_val = _get_Xs_ys_ids(args["validation_file"],
        args["debugging"], True)

    # Get DataLoaders
    train_loader = _get_data_loader(Xs_train, ys_train, args["batch_size"])
    val_loader = _get_data_loader(Xs_val, ys_val, args["batch_size"])
    data_loaders = dict({"train": train_loader, "validation": val_loader})

    # Get sequences
    seqs = {}
    train_seqs = _get_sequences(args["training_file"])
    seqs.setdefault("train", train_seqs)
    seqs["train"] = np.concatenate(
        (seqs["train"], rc_many(train_seqs)), axis=0
    )
    val_seqs = _get_sequences(args["validation_file"])
    seqs.setdefault("validation", val_seqs)
    seqs["validation"] = np.concatenate(
        (seqs["validation"], rc_many(val_seqs)), axis=0
    )

    ##############
    # Interpret  #
    ############## 

    # Initialize
    activations = {}
    labels = {}
    outputs = {}
    predictions = {}
    if np.unique(ys_train[:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"
    model = _load_model(args["model_file"])

    # Create output dirs
    if not os.path.isdir(args["output_dir"]):
        os.makedirs(args["output_dir"])
    for subdir in ["sites", "motifs", "logos"]:
        if not os.path.isdir(os.path.join(args["output_dir"], subdir)):
            os.makedirs(os.path.join(args["output_dir"], subdir))

    # Get activations, labels, outputs, and predictions
    for train_validation in ["train", "validation"]:
        data_loader = data_loaders[train_validation]
        acts_arr, labs_arr, outs_arr, preds_arr = _get_acts_labs_outs_preds(
            model, data_loader, input_data)
        activations.setdefault(train_validation, acts_arr)
        labels.setdefault(train_validation, labs_arr)
        outputs.setdefault(train_validation, outs_arr)
        predictions.setdefault(train_validation, preds_arr)

    # Get the indices of well predicted sequences
    idxs = _get_well_predicted_sequences(predictions["train"], labels["train"],
        input_data)

    # For each filter, get the activation threshold (i.e. ≥50%)
    thresholds = 0.5 * np.amax(activations["train"][idxs, :, :], axis=(0, 2))

    # For each filter...
    for i in tqdm(range(model._options["cnn_units"]), bar_format=bar_format,
        total=model._options["cnn_units"]):

        # Get sites
        sites_file = os.path.join(args["output_dir"], "sites",
            f"filter{i}.fa.gz")
        if not os.path.exists(sites_file):
            with gzip.open(sites_file, "wt") as handle:
                _get_sites(handle, idxs, seq_ids_train, seqs["train"],
                    activations["train"][:, i, :], thresholds[i],
                    model._options["kernel_size"])
 
        # Get motif
        motif_file = os.path.join(args["output_dir"], "motifs",
            f"filter{i}.jaspar")
        if not os.path.exists(motif_file):
            motif = _sites_to_motif(sites_file)
            motif.matrix_id = f"filter{i}"
            motif.name = args["name"]
            with open(motif_file, "wt") as handle:
                handle.write(format(motif, "jaspar"))

        # Get logos
        for reverse_complement in [False, True]:
            if reverse_complement:
                logo_file = os.path.join(args["output_dir"], "logos",
                    f"filter{i}.rev.png")
            else:
                logo_file = os.path.join(args["output_dir"], "logos",
                    f"filter{i}.fwd.png")
            if not os.path.exists(logo_file):
                fig = get_figure(motif_file, reverse_complement)
                fig.savefig(logo_file, bbox_inches="tight", pad_inches=0)

    # Get weights
    weights = model.final.weight.detach().cpu().numpy()
    weights_file = os.path.join(args["output_dir"], "weights.tsv")
    if not os.path.exists(weights_file):
        data = []
        for i, weight in enumerate(weights.T):
            data.append([f"filter{i}"] + weight.tolist())
        df = pd.DataFrame(data)
        df.to_csv(weights_file, sep="\t", index=False)

    # Get feature effects
    dfs = []
    kwargs = {"total": weights.shape[0], "bar_format": bar_format}
    feat_effects_file = os.path.join(args["output_dir"], "feat-effects.tsv.gz")
    if not os.path.exists(feat_effects_file):
        for i in tqdm(range(weights.shape[0]), **kwargs):
            feat_effects = outputs["train"] * weights[i]
            c = [f"filter{j}" for j in range(model._options["cnn_units"])]
            df = pd.DataFrame(feat_effects[idxs, :].tolist(), columns=c)
            df["class"] = i
            dfs.append(df)
        df = pd.concat(dfs)
        df.to_csv(feat_effects_file, sep="\t", index=False, compression="gzip")

    # Get motifs in MEME format
    meme_file = os.path.join(args["output_dir"], "motifs", "filters.meme")
    if not os.path.exists(meme_file):
        jaspar_motifs = []
        for i in range(model._options["cnn_units"]):
            motif_file = os.path.join(args["output_dir"], "motifs",
                f"filter{i}.jaspar")
            with open(motif_file) as handle:
                jaspar_motifs.append(motifs.read(handle, "jaspar"))
        reformat_motif(jaspar_motifs, "meme", meme_file)

    # Get performances
    metrics = get_metrics(input_data=input_data)
    kwargs = {"total": model._options["cnn_units"], "bar_format": bar_format}
    for m in metrics:
        data = []
        performance_file = os.path.join(args["output_dir"], f"{m}.tsv")
        for i in tqdm(range(model._options["cnn_units"]), **kwargs):
            scores = _get_filter_performances(
                activations["validation"][:, i, :], ys_val, metrics[m])
            data.append([f"filter{i}"] + scores)
        df = pd.DataFrame(data)
        df.to_csv(performance_file, sep="\t", index=False, compression="gzip")

def _get_sequences(tsv_file):

    # Get sequences
    df = pd.read_csv(tsv_file, sep="\t", header=None)

    return(df[1].values)

def _load_model(model_file):

    # Initialize
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Get model
    selene_dict = torch.load(model_file)
    model = CAM(
        selene_dict["options"]["cnn_units"],
        selene_dict["options"]["kernel_size"],
        selene_dict["options"]["sequence_length"],
        selene_dict["options"]["n_features"],
        selene_dict["options"]["weights_file"],
    )
    model.load_state_dict(selene_dict["state_dict"])
    model.to(device)

    return(model)

def _get_acts_labs_outs_preds(model, data_loader, input_data):

    # Counter
    idx = 0
    x = len(data_loader.dataset)
    y = model._options["cnn_units"]
    z = model._options["sequence_length"] - model._options["kernel_size"] + 1
    n_features = model._options["n_features"]
    activations = torch.zeros((x, y, z), dtype=torch.float32)
    labels = torch.zeros((x, n_features), dtype=torch.float32)
    outputs = torch.zeros((x, y), dtype=torch.float32)
    predictions = torch.zeros((x, n_features), dtype=torch.float32)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    with torch.no_grad():
        for seqs, labs in tqdm(iter(data_loader), total=len(data_loader),
                bar_format=bar_format):

            # Prepare input
            seqs = seqs.to(device)
            labs = labs.to(device)
            seqs = seqs.repeat(1, model._options["cnn_units"], 1)

            # Get labels
            labels[idx:idx+seqs.shape[0]] = labs.cpu()

            # Get outputs
            outs = model.linears(seqs)
            outputs[idx:idx+seqs.shape[0], :] = outs.cpu()

            # Get predictions
            preds = model.final(outs)
            if input_data == "binary":
                preds = torch.sigmoid(preds)
            predictions[idx:idx+seqs.shape[0]] = preds.cpu()

            # Get activations
            activations[idx:idx+seqs.shape[0], :, :] = \
                model.linears[:3](seqs).cpu()

            # Index increase
            idx += seqs.shape[0]

    return(activations.numpy(), labels.numpy(), outputs.numpy(),
        predictions.numpy())

def _get_well_predicted_sequences(predictions, labels, input_data, top=.1):

    # Initialize
    n = labels[0].shape[0]

    # For binary data such as ChIP-seq, well-predicted sequences are bound
    # regions (i.e. 1) with a score greater than 0.5, and unbound regions
    # (i.e. 0) with a score lower or equal than 0.5.
    if input_data == "binary":
        c = Counter(np.where(labels == (predictions > .5).astype(int))[0])
        return(np.sort(np.array([k for k, v in c.items() if v == n])))
    # For non-binary data such as PBM, well-predicted sequences are defined 
    # as the top `n` percentile probes with the highest signal intensities
    # and a score in the top `n` percentile.
    else:
        l = np.argsort(-labels.flatten())[:int(max(labels.shape)*top)]
        o = np.argsort(-outputs.flatten())[:int(max(outputs.shape)*top)]
        idxs = np.intersect1d(l, o)

    return(idxs)

def _get_sites(handle, idxs, seq_ids, seqs, activations, threshold,
               kernel_size=19):
    """
    For each filter and each sequence, get sites reaching at least ½ of the
    max. activation value for that filter.
    """

    # Initialize
    count = 0
    strands = ["+"] * (seqs.size // 2) + ["-"] * (seqs.size // 2)

    # For each sequence...
    for i in idxs:

        # Get start positions of sequence sites activating this filter
        starts = np.where(activations[i, :] > threshold)

        # For each start...
        for j in starts[0]:

            # Get site
            start = j
            end = j+kernel_size
            strand = strands[i]
            seq = Seq(seqs[i][start:end])
            seq_id = f"{seq_ids[i][0]}_{strand}_from={start}_to={end}"
            record = SeqRecord(seq, id=seq_id, name="", description="")
            handle.write(record.format("fasta"))

            # If count reaches 1M, this filter is way too ubiquitous!!!
            count += 1
            if count == 1e6:
                return

def _sites_to_motif(sites_file):
    """
    From https://github.com/biopython/biopython/blob/master/Bio/motifs/__init__.py
    Read the motif from JASPAR .sites file.
    """

    # Initialize
    alphabet = "ACGTN"
    instances = []
    pfm = {}

    with gzip.open(sites_file, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                break
            # line contains the header ">...."
            # now read the actual sequence
            line = next(handle)
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

def _get_filter_performances(activations, ys, metric):
   
    # Initialize
    scores = []
    ys = ys[:len(ys)//2]

    # Get max. activations on fwd and rev strand
    fwd = activations[:len(activations)//2]
    rev = activations[len(activations)//2:]
    activations = np.concatenate((fwd, rev), axis=1)
    max_activations = np.max(activations, axis=1)

    # Get performance for each class
    for y in ys.T:
        scores.append(metric(y, max_activations))

    return(scores)

if __name__ == "__main__":
    main()