#!/usr/bin/env python

from Bio import SeqIO
from Bio import motifs
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import click
import gc
import numpy as np
import os
import torch
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Local imports
from architectures import CAM
from jaspar import get_figure, reformat_motif
from sequence import one_hot_encode, rc_one_hot_encoding, rc
from train import _get_data_loaders, __get_handle

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
    "-d", "--debugging",
    help="Debugging mode.",
    is_flag=True,
)
@click.option(
    "-n", "--name",
    help="Name for Biopython's {Motif} class.",
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
    help="Reverse complement training sequences.",
    is_flag=True,
)
@click.option(
    "-t", "--threads",
    help="Number of CPU threads to use.",
    type=int,
    default=1,
    show_default=True,
)

def main(**params):

    ##############
    # Initialize #
    ############## 

    # Initialize
    torch.set_num_threads(params["threads"])

    # Create output dirs
    if not os.path.isdir(params["output_dir"]):
        os.makedirs(params["output_dir"])
    for subdir in ["sites", "motifs", "logos"]:
        if not os.path.isdir(os.path.join(params["output_dir"], subdir)):
            os.makedirs(os.path.join(params["output_dir"], subdir))

    # Get model
    selene_dict = torch.load(params["model_file"])
    model = CAM(
        selene_dict["options"]["cnn_units"],
        selene_dict["options"]["kernel_size"],
        selene_dict["options"]["sequence_length"],
        selene_dict["options"]["n_features"],
        selene_dict["options"]["clamp_weights"],
        selene_dict["options"]["no_padding"],
        selene_dict["options"]["weights_file"],
    )
    model.load_state_dict(selene_dict["state_dict"])
    model.to(device)

    ##############
    # Load Data  #
    ##############

    # Get data
    Xs, ys, seq_ids, seqs = _get_Xs_ys_seq_ids_seqs(params["training_file"],
        params["debugging"], params["rev_complement"])

    # Get DataLoader
    train_loader = _get_data_loaders(list(Xs), list(ys),
        # batch_size=params["batch_size"], threads=params["threads"])
        batch_size=params["batch_size"])
    data_loader = train_loader

    ##############
    # Interpret  #
    ############## 

    # Initialize
    jaspar_motifs = []
    if np.unique(ys[:, 0]).size == 2:
        input_data = "binary"
    else:
        input_data = "linear"
    if selene_dict["options"]["no_padding"]:
        padding = 0
    else:
        padding = selene_dict["options"]["kernel_size"]

    # Fix sequences
    for i in range(len(seqs)):
        seqs[i] = "N" * padding + seqs[i] + "N" * padding
    seqs = np.array(seqs)

    # Get predict
    outputs, labels = _predict(model, data_loader, input_data)

    # Get activations
    x = len(seqs)
    y = model._options["cnn_units"]
    z = len(seqs[0]) - model._options["kernel_size"] + 1
    size = (x, y, z)
    activations = _get_activations(model, data_loader,
        torch.zeros(size, dtype=torch.float32))

    # Get the indices of well predicted sequences
    if input_data == "binary":
        idxs = np.where((labels == 1.) & (outputs >= .5))[0]
    else:
        l_idxs = np.argsort(-labels.flatten())[:int(max(labels.shape) * .1)]
        o_idxs = np.argsort(-outputs.flatten())[:int(max(outputs.shape) * .1)]
        idxs = np.intersect1d(l_idxs, o_idxs)

    # Free memory
    del Xs
    del ys
    del outputs
    del labels
    del data_loader
    gc.collect()

    # For each filter, get the activation threshold (i.e. ≥50%)
    thresholds = 0.5 * np.amax(activations[idxs, :, :], axis=(0, 2))

    # For each filter...
    for i in tqdm(range(selene_dict["options"]["cnn_units"]),
                  total=selene_dict["options"]["cnn_units"],
                  bar_format=bar_format):

        # Get sites
        sites_file = os.path.join(params["output_dir"], "sites",
            f"filter{i}.fa.gz")
        if not os.path.exists(sites_file):
            handle = __get_handle(sites_file, "wt")
            _get_sites(handle, idxs, seq_ids, seqs, activations[:, i, :],
                thresholds[i], model._options["kernel_size"])
            handle.close()

        # Get motif
        motif_file = os.path.join(params["output_dir"], "motifs",
            f"filter{i}.jaspar")
        if not os.path.exists(motif_file):
            handle = __get_handle(sites_file)
            motif = _get_motif(handle)
            handle.close()
            motif.matrix_id = f"filter{i}"
            motif.name = params["name"]
            handle = __get_handle(motif_file, "wt")
            handle.write(format(motif, "jaspar"))
            handle.close()
            jaspar_motifs.append(motif)

        # Get logos
        for reverse_complement in [False, True]:
            if reverse_complement:
                logo_file = os.path.join(params["output_dir"], "logos",
                    f"filter{i}.rev.png")
            else:
                logo_file = os.path.join(params["output_dir"], "logos",
                    f"filter{i}.fwd.png")
            if not os.path.exists(logo_file):
                fig = get_figure(motif_file, reverse_complement)
                fig.savefig(logo_file, bbox_inches="tight", pad_inches=0)

    # Get motifs in MEME format
    meme_file = os.path.join(params["output_dir"], "motifs", "filters.meme")
    if not os.path.exists(meme_file):
        reformat_motif(jaspar_motifs, "meme", meme_file)

    # Get weights
    weights_file = os.path.join(params["output_dir"], "weights.tsv")
    if not os.path.exists(weights_file):
        weights = model.final.weight.detach().cpu().numpy()
        handle = __get_handle(weights_file, "wt")
        for i, weight in enumerate(weights.T):
            s = "\t".join(map(str, weight))
            handle.write(f"filter{i}\t{s}\n")
        handle.close()

def _get_Xs_ys_seq_ids_seqs(fasta_file, debugging=False,
    reverse_complement=False):

    # Initialize
    Xs = []
    ys = []
    seq_ids = []
    seqs = []

    # Xs / ys
    handle = __get_handle(fasta_file)
    for record in SeqIO.parse(handle, "fasta"):
        _, y_list = record.description.split()
        Xs.append(one_hot_encode(str(record.seq).upper()))
        ys.append([float(y) for y in y_list.split(";")])
        seq_ids.append((record.id, "+"))
        seqs.append(str(record.seq))
    handle.close()

    # Reverse complement
    if reverse_complement:
        n = len(Xs)
        for i in range(n):
            Xs.append(rc_one_hot_encoding(Xs[i]))
            ys.append(ys[i])
            seq_ids.append((seq_ids[i][0], "-"))
            seqs.append(rc(seqs[i]))

    # Return 1,000 sequences
    if debugging:
        return(np.array(Xs)[:10000], np.array(ys)[:10000], 
               np.array(seq_ids)[:10000])

    return(np.array(Xs), np.array(ys), np.array(seq_ids), seqs)

def _predict(model, data_loader, input_data):

    # Initialize
    outputs = []
    labels = []

    with torch.no_grad():
        for x, label in tqdm(iter(data_loader), total=len(data_loader),
                             bar_format=bar_format):

            # Get outputs
            x = x.to(device)
            out = model(x)
            if input_data == "binary":
                out = torch.sigmoid(out)
            outputs.extend(out.detach().cpu().numpy())

            # Get labels
            labels.extend(label.numpy())

    return(np.array(outputs), np.array(labels))

def _get_activations(model, data_loader, activations):

    # Counter
    idx = 0

    with torch.no_grad():
        for x, _ in tqdm(iter(data_loader), total=len(data_loader),
                bar_format=bar_format):

            # Get activations
            x = x.to(device)
            x = x.repeat(1, model._options["cnn_units"], 1)
            activations[idx:idx+x.shape[0], :, :] = model.linears[:3](x).cpu()
            idx += x.shape[0]           

    return(activations.numpy())

def _get_sites(handle, idxs, seq_ids, seqs, activations, threshold,
               kernel_size=19):
    """
    For each filter and each sequence, get sites reaching at least ½ of the
    maximum activation value for that filter.
    """

    # For each sequence...
    for i in idxs:

        # Get start positions of sequence sites activating this filter
        starts = np.where(activations[i, :] > threshold)

        # For each start...
        for j in starts[0]:

            # Get site
            record_id = seq_ids[i][0]
            strand = seq_ids[i][1]
            start = j
            end = j+kernel_size
            seq = Seq(seqs[i][start:end])
            seq_id = "%s_%s_from=%s_to=%s" % (record_id, strand, start, end)
            record = SeqRecord(seq, id=seq_id, name="", description="")
            handle.write(record.format("fasta"))

def _get_motif(handle):
    """
    From https://github.com/biopython/biopython/blob/master/Bio/motifs/__init__.py
    Read the motif from JASPAR .sites file.
    """

    # Initialize
    alphabet = "ACGTN"
    instances = []
    pfm = {}

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

if __name__ == "__main__":
    main()