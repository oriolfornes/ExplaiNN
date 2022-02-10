#!/usr/bin/env python

from Bio import SeqIO
from Bio import motifs
import click
from click_option_group import optgroup
from collections import Counter
import copy
import gzip
import numpy as np
import os
from sklearn import metrics
import time
import torch
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Local imports
from damo.helpers import *
from sequence import one_hot_encode


# Globals
device = "cpu"
if torch.cuda.is_available():
    device = "cuda"

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "positive_seqs",
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "negative_seqs",
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "jaspar_profile",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-d", "--debugging",
    help="Debugging mode.",
    is_flag=True,
)
@click.option(
    "-o", "--output-file",
    help="Output file.",
    type=click.Path(resolve_path=True),
    default="./DAMO.txt",
    show_default=True,
)
@click.option(
    "-t", "--threads",
    help="Number of CPU threads to use.",
    type=int,
    default=1,
    show_default=True,
)
@optgroup.group("DAMO")
@optgroup.option(
    "-i", "--iterations",
    help="Number of optimization iterations.",
    type=int,
    default=500,
    show_default=True,
)

def main(**args):

    damo(args["positive_seqs"], args["negative_seqs"], args["jaspar_profile"],
         args["debugging"], args["output_file"], args["threads"],
         args["iterations"])

def damo(positive_seqs_file, negative_seqs_file, jaspar_profile_file,
         debugging=False, output_file=None, threads=1, iterations=500):
    """DAMO (Discriminative Additive Model Optimization)."""

    # Initialize
    LEARNING_RATES = [1.0, 0.55, 0.1]

    start = time.time()

    # Get sequences
    positive_seqs = _get_sequences(positive_seqs_file, debugging)
    negative_seqs = _get_sequences(negative_seqs_file, debugging)
    positive_seqs = [str(s.seq).upper() for s in positive_seqs] + \
        [str(s.reverse_complement().seq).upper() for s in positive_seqs]
    negative_seqs = [str(s.seq).upper() for s in negative_seqs] + \
        [str(s.reverse_complement().seq).upper() for s in negative_seqs]

    # One-hot encode sequences
    positive_1he_seqs = np.array(_one_hot_encode(positive_seqs))
    negative_1he_seqs = np.array(_one_hot_encode(negative_seqs))

    # Sequences to Tensor
    positive_seqs_tensor = torch.Tensor(positive_1he_seqs.tolist()).to(device)
    negative_seqs_tensor = torch.Tensor(negative_1he_seqs.tolist()).to(device)

    # Get PWM
    pwm = _get_PWM(jaspar_profile_file)

    # Get a PyTorch model of the PWM
    model = _get_PWM_model(pwm, positive_1he_seqs.shape[-1], device)

    end = time.time()
    print(end-start)

    start = time.time()

    # Score sequences
    positive_scores, positive_sites = _score_sequences(
        model, positive_seqs_tensor, positive_seqs, pwm.shape[1]
    )
    negative_scores, negative_sites = _score_sequences(
        model, negative_seqs_tensor, negative_seqs, pwm.shape[1],
    )

    # Compute AUC
    original_auc, y_true, y_score = _compute_AUC(positive_scores,
        negative_scores)
    auc = copy.copy(original_auc)

    # For each iteration...
    for _ in tqdm(range(iterations), total=iterations, bar_format=bar_format):

        # Initialize
        pwm_0 = copy.copy(pwm)
        auc_0 = copy.copy(auc)
        y_true_0 = copy.copy(y_true)
        y_score_0 = copy.copy(y_score)
        sites_0 = copy.copy(positive_sites) + copy.copy(negative_sites)

        # For each learning rate...
        for lr in LEARNING_RATES:

            # Update PWM
            pwm_1 = _update_PWM(y_true_0, y_score_0, sites_0, lr, pwm_0)

            # Get a PyTorch model of the updated PWM
            model = _get_PWM_model(pwm_1, positive_1he_seqs.shape[-1], device)

            # Score sequences
            positive_scores, positive_sites = _score_sequences(
                model, positive_seqs_tensor, positive_seqs, pwm_1.shape[1]
            )
            negative_scores, negative_sites = _score_sequences(
                model, negative_seqs_tensor, negative_seqs, pwm_1.shape[1],
            )

            # Compute AUC
            auc, y_true, y_score = _compute_AUC(positive_scores,
                negative_scores)

            # Early stop
            if auc > auc_0:
                pwm = copy.copy(pwm_1)
                break

        # Early stop
        if np.array_equal(pwm_0, pwm):
            break

    end = time.time()
    print(end-start)
    print(f"AUC: {original_auc} >>> {auc}")

    print(pwm)

def __get_handle(file_name, mode="rt"):
    if file_name.endswith(".gz"):
        handle = gzip.open(file_name, mode)
    else:
        handle = open(file_name, mode)

    return(handle)

def _get_sequences(fasta_file, debugging=False):

    # Initialize
    sequences = []

    # Get sequences
    handle = __get_handle(fasta_file)
    for record in SeqIO.parse(handle, "fasta"):
        sequences.append(record)
    handle.close()

    # Return 1,000 sequences
    if debugging:
        return(sequences[:1000])

    return(sequences)

def _get_PWM(jaspar_profile):

    # Initialize
    pwm = []

    # Get PWM
    handle = __get_handle(jaspar_profile)
    m = motifs.read(handle, "jaspar")
    handle.close()
    m.pseudocounts = motifs.jaspar.calculate_pseudocounts(m)
    for nucleotide in "ACGT":
        pwm.append([i for i in m.pwm[nucleotide]])

    return(np.log(pwm))

def _one_hot_encode(sequences):
    return([one_hot_encode(s) for s in sequences])

def _get_PWM_model(pwm, sequence_length, device="cpu"):

    # Get PWM model
    conv1d = torch.nn.Conv1d(4, 1, kernel_size=pwm.shape[1])
    conv1d.bias.data = torch.Tensor([0.])
    conv1d.weight.data = torch.Tensor([pwm])
    for p in conv1d.parameters():
        p.requires_grad = False
    maxpool1d = torch.nn.MaxPool1d(sequence_length - pwm.shape[1] + 1, 1,
        return_indices=True)

    return(torch.nn.Sequential(conv1d, maxpool1d).to(device))

def _score_sequences(model, sequences_tensor, sequences, motif_length):

    # Initialize
    scores = []
    sites = []

    # Get predictions
    with torch.no_grad():
        scores_tmp, sites_tmp = model(sequences_tensor)
    scores_tmp = scores_tmp.flatten().tolist()
    sites_tmp = sites_tmp.flatten().tolist()

    # Zip predictions from the same sequences (i.e. fwd and rev strands)  
    z = list(zip(scores_tmp, sites_tmp, list(range(len(sequences)))))
    z = list(zip(z[:len(z)//2], z[len(z)//2:]))

    # For each prediction...
    for fwd, rev in z:
        score, site = __get_best_site(fwd, rev, sequences, motif_length)
        scores.append(score)
        sites.append(site)

    return(scores, sites)

def __get_best_site(fwd, rev, sequences, motif_length):

    # If fwd is better...
    if fwd[0] > rev[0]:
        score = fwd[0]
        site = sequences[fwd[2]][fwd[1]:fwd[1]+motif_length]
    # If rev is better...
    else:
        score = rev[0]
        site = sequences[rev[2]][rev[1]:rev[1]+motif_length]

    return(score, site)

def _compute_AUC(positive_scores, negative_scores):

    # Initialize
    y_true = [1] * len(positive_scores) + [0] * len(negative_scores)
    y_score = positive_scores + negative_scores

    return(metrics.roc_auc_score(y_true, y_score), y_true, y_score)

def _update_PWM(y_true, y_score, sites, lr, pwm):
    """makeAMove."""

    # Initialize
    positive_mix = []
    negative_mix = []

    # Sort and zip sequences and labels
    idxs = np.argsort(-np.asarray(y_score), kind="mergesort")
    sequences = [sites[ix] for ix in idxs]
    labels = [y_true[ix] for ix in idxs]
    z = list(zip(sequences, labels))

    # Get index of first negative and last positive
    ix_first_negative = labels.index(0)
    ix_last_positive = len(labels) - labels[::-1].index(1)

    # For each label, sequence...
    for s, l in z[ix_first_negative:ix_last_positive]:
        positive_mix.append(s) if l == 1 else negative_mix.append(s)

    # Generate PWMs
    positive_pwm = __generate_PWM(positive_mix, pwm.shape[1])
    negative_pwm = __generate_PWM(negative_mix, pwm.shape[1])
    exit(0)

    return(pwm + lr * (positive_pwm - negative_pwm))

def __generate_PWM(sequences, motif_length):
    """gen_pwm."""

    # Get PFM
    pfm = np.zeros((motif_length, 4))
    for i in range(motif_length):
        c = Counter(s[i:i+1] for s in sequences)
        pfm[i, :] = [c[key] for key in "ACGT"]
    pfm = PFM.Normalize(PFM.Normalize(pfm), 1e-4)

    return(np.log(pfm).T)


def __score_perceptron(old_pwm, old_scores, old_sites, sequences, labels,
    positive_seqs, negative_seqs):
    """Adapted from DAMO."""

    scored_positive, pos_sites = scoreAllSeq(read_seq.pos_splitted, pwm, enc)
    scored_negative, neg_sites = scoreAllSeq(read_seq.neg_splitted, pwm, enc)
    auc, y_true, y_score = cal_auc(scored_positive, scored_negative)

    for learning_rate in LEARNING_RATES:
        new_pwm = makeAMove(y_true, y_score, pos_sites + neg_sites, learning_rate, pwm, level, enc)

        temp_scored_positive, _ = scoreAllSeq(read_seq.pos_splitted, new_pwm, enc)
        temp_scored_negative, _ = scoreAllSeq(read_seq.neg_splitted, new_pwm, enc)
        temp_auc = cal_auc(temp_scored_positive, temp_scored_negative)[0]

        if temp_auc > auc:
            return new_pwm

    return pwm


if __name__ == "__main__":
    main()