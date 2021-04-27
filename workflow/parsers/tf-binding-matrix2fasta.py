#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click
from functools import partial
import gzip
import hashlib
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import pickle
from pybedtools import BedTool
from pybedtools.helpers import cleanup
import random
import shutil
import sparse
import subprocess as sp
from tqdm import tqdm
bar_format = "{percentage:3.0f}%|{bar:20}{r_bar}"

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "tf_binding_matrix_dir", type=click.Path(exists=True, resolve_path=True)
)
@click.argument(
    "genome_file", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-d", "--dummy-dir",
    help="Dummy directory.",
    type=click.Path(resolve_path=True),
    default="/tmp/",
    show_default=True
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True
)
@click.option(
    "-t", "--threads",
    help="Threads to use.",
    type=int,
    default=1,
    show_default=True
)

def main(**params):

    # Create output dirs
    if not os.path.isdir(params["output_dir"]):
        os.makedirs(params["output_dir"])
    for split in ["Train", "Validation"]:
        if not os.path.isdir(os.path.join(params["output_dir"], split)):
            os.makedirs(os.path.join(params["output_dir"], split))

    # Get regions indices
    global regions
    f = os.path.join(params["tf_binding_matrix_dir"], "regions_idx.pickle.gz")
    with gzip.open(f, "rb") as handle:
        regions = np.array(list(pickle.load(handle).keys()))

    # Get TF indices
    global tfs_idx
    f = os.path.join(params["tf_binding_matrix_dir"], "tfs_idx.pickle.gz")
    with gzip.open(f, "rb") as handle:
        tfs_idx = pickle.load(handle)

    # Get 3D matrix
    f = os.path.join(params["tf_binding_matrix_dir"], "matrix3d.npz")
    matrix3d = sparse.load_npz(f)

    # Collapse into a 1817918x163 matrix
    global matrix2d
    matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")

    # ReMap-bound and open regions (4, 2) = 1
    # Open regions (1), UniBind-only bound and open regions (3) = 0
    # Remaining regions (0) = None
    matrix2d[matrix2d == 0] = np.nan
    matrix2d[matrix2d == 1] = 0
    matrix2d[matrix2d == 2] = 1
    matrix2d[matrix2d == 3] = 0
    matrix2d[matrix2d == 4] = 1

    # Get chrom sizes
    global chrom_sizes
    chrom_sizes = __get_chrom_sizes(params["genome_file"])

    # Get FASTA sequences
    kwargs = {"total": len(tfs_idx), "bar_format": bar_format}
    pool = Pool(params["threads"])
    p = partial(__get_FASTA_sequences, genome_file=params["genome_file"],
        dummy_dir=params["dummy_dir"], output_dir=params["output_dir"])
    for _ in tqdm(pool.imap(p, sorted(tfs_idx)), **kwargs):
        pass

def __get_chrom_sizes(genome_file):

    # Initialize
    chrom_sizes = {}

    with open("%s.sizes" % genome_file) as handle:
        for line in handle:
            chrom, size = line.strip("\n").split("\t")
            chrom_sizes.setdefault(chrom, tuple([0, size]))

    return(chrom_sizes)

def __get_FASTA_sequences(tf, genome_file, dummy_dir="/tmp/",
    output_dir="./"):

    # Initialize
    m = hashlib.md5()

    # Get Train/Validation splits
    train_chroms = set(["chr%s" % (c+1) for c in range(22) if (c % 2) == 0])
    val_chroms = set(["chr%s" % (c+1) for c in range(22) if (c % 2) != 0])
    s = "\n".join(["%s\t%s\t%s" % (k, v[0], v[1]) \
        for k, v in chrom_sizes.items() if k in train_chroms])
    t = BedTool(s, from_string=True)
    s = "\n".join(["%s\t%s\t%s" % (k, v[0], v[1]) \
        for k, v in chrom_sizes.items() if k in val_chroms])
    v = BedTool(s, from_string=True)

    # Get TF matrix
    matrix1d = matrix2d[:,tfs_idx[tf]]
    pos_ixs = np.where(matrix1d == 1.)[0]
    neg_ixs = np.where(matrix1d == 0.)[0]

    # Positive sequences
    positives = {}
    s = "\n".join(["\t".join(r) for r in regions[pos_ixs]])
    a = BedTool(s, from_string=True)
    for split in ["Train", "Validation"]:
        positives.setdefault(split, [])
        if split == "Train":
            i = a.intersect(t, wa=True)
        else:
            i = a.intersect(v, wa=True)
        i.sequence(fi=genome_file)
        pos_dummy_file = os.path.join(dummy_dir, "%s.%s.pos.fa" % (tf, split))
        with open(pos_dummy_file, "wt") as handle:
            for s in SeqIO.parse(i.seqfn, "fasta"):
                handle.write(s.format("fasta"))
                record = SeqRecord(s.seq, id=s.id, description="1.")
                positives[split].append(record)
        random.shuffle(positives[split])

    # Negative sequences with matched %GC distribution
    negatives = []
    biasaway_negs = []
    negatives_ix = {}
    s = "\n".join(["\t".join(r) for r in regions[neg_ixs]])
    b = BedTool(s, from_string=True)
    for split in ["Train", "Validation"]:
        if split == "Train":
            i = b.intersect(t)
        else:
            i = b.intersect(v)
        i.sequence(fi=genome_file)
        neg_dummy_file = os.path.join(dummy_dir, "%s.neg.fa" % tf)
        with open(neg_dummy_file, "wt") as handle:
            for s in SeqIO.parse(i.seqfn, "fasta"):
                handle.write(s.format("fasta"))
                record = SeqRecord(s.seq, id=s.id, description="0.")
                negatives.append(record)
                md5 = hashlib.md5(str(s.seq).upper().encode()).hexdigest()
                negatives_ix.setdefault(md5, len(negatives) - 1)
        pos_dummy_file = os.path.join(dummy_dir, "%s.%s.pos.fa" % (tf, split))
        biasaway_file = os.path.join(dummy_dir, "%s.biasaway.fa" % tf)
        cmd = "biasaway g -f %s -r %s -b %s -e 1 > %s" % (pos_dummy_file,
            os.path.join(dummy_dir, tf), neg_dummy_file, biasaway_file)
        _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL)
        for s in SeqIO.parse(biasaway_file, "fasta"):
            md5 = hashlib.md5(str(s.seq).upper().encode()).hexdigest()
            biasaway_negs.append(negatives[negatives_ix[md5]])
        random.shuffle(biasaway_negs)

        # Delete tmp files
        os.remove(pos_dummy_file)
        os.remove(neg_dummy_file)
        os.remove(biasaway_file)
        shutil.rmtree(os.path.join(dummy_dir, tf))

        # Save sequences
        sequences_file = os.path.join(output_dir, split, "%s.fa.gz" % tf)
        pos_seqs = positives[split]
        with gzip.open(sequences_file, "wt") as handle:
            sequences = [j for i in zip(pos_seqs, biasaway_negs) for j in i]
            SeqIO.write(sequences, handle, "fasta")

    # Delete tmp files
    cleanup()

if __name__ == "__main__":
    main()