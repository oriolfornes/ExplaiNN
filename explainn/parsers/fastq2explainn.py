#!/usr/bin/env python

from Bio import SeqIO
import click
import os
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
import subprocess as sp

from utils import get_file_handle

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "fastq_file", type=click.Path(exists=True, resolve_path=True), nargs=-1,
)
@click.option(
    "--clip-left",
    help="Left-clip reads.",
    default=0,
    show_default=True
)
@click.option(
    "--clip-right",
    help="Right-clip reads.",
    default=0,
    show_default=True
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
    "-p", "--prefix",
    help="Output prefix.",
)
@click.option(
    "--paired-end",
    help="Paired-end reads.",
    is_flag=True,
)
@click.option(
    "-s", "--sort",
    help="Sort FASTQ files.",
    is_flag=True,
)
@click.option(
    "-t", "--test",
    help="Create a test set.",
    is_flag=True,
)

def main(**args):

    # Create output dir
    if not os.path.exists(args["output_dir"]):
        os.makedirs(args["output_dir"])

    # Get TSV files for ExplaiNN
    _to_ExplaiNN(args["fastq_file"], args["clip_left"], args["clip_right"],
        args["dummy_dir"], args["output_dir"], args["prefix"],
        args["paired_end"], args["sort"], args["test"])

def _to_ExplaiNN(fastq_files, clip_left=0, clip_right=0, dummy_dir="/tmp/",
    output_dir="./", prefix=None, paired_end=False, sort=False,
    create_test=False):

    # Initialize
    data = []
    if paired_end:
        divisor = 2
    else:
        divisor = 1
    if sort:
        sorted_files = list(sorted(fastq_files))
    else:
        sorted_files = fastq_files

    # Ys
    enc = OneHotEncoder()
    arr = [[i] for i in range(len(sorted_files)//divisor)]
    enc.fit(arr)
    ys = enc.transform(arr).toarray().tolist()

    # Get DataFrame
    for i in range(len(sorted_files)):
        handle = get_file_handle(sorted_files[i], "rt")
        for record in SeqIO.parse(handle, "fastq"):
            y = ys[i//divisor]
            s = str(record.seq).upper()[clip_left:-clip_right]
            data.append([record.id, s] + y)
        handle.close()
    df = pd.DataFrame(data)
    df = df.groupby(1).max().reset_index()
    df = df.reindex(sorted(df.columns), axis=1)

    # Generate negative sequences by dinucleotide shuffling
    if df.shape[1] == 3: # i.e. only one class
        data = []
        dummy_file = os.path.join(dummy_dir, "%s+%s+%s.fa" %
            (os.path.split(__file__)[1], str(os.getpid()), prefix))
        with open(dummy_file, "wt") as handle:
            for z in zip(df.iloc[:, 0]. tolist(), df.iloc[:, 1]. tolist()):
                handle.write(f">{z[0]}\n{z[1]}\n")
        cmd = "biasaway k -f %s -k 2 -e 1 > %s.biasaway" % (dummy_file,
            dummy_file)
        _ = sp.run([cmd], shell=True, cwd=scripts_dir, stderr=sp.DEVNULL)
        for s in SeqIO.parse("%s.biasaway" % dummy_file, "fasta"):
            header = "%s::shuf" % s.description.split(" ")[-1]
            data.append([header, str(s.seq), 0.])
        df = pd.concat((df, pd.DataFrame(data)))
        os.remove(dummy_file)
        os.remove("%s.biasaway" % dummy_file)

    # Get data splits
    train, validation = train_test_split(df, test_size=0.2)
    if create_test:
        validation, test = train_test_split(validation, test_size=0.5)

    # Save sequences
    if prefix is None:
        tsv_file = os.path.join(output_dir, "train.tsv.gz")
    else:
        tsv_file = os.path.join(output_dir, f"{prefix}.train.tsv.gz")
    train.to_csv(tsv_file, sep="\t", header=False, index=False,
        compression="gzip")
    if prefix is None:
        tsv_file = os.path.join(output_dir, "validation.tsv.gz")
    else:
        tsv_file = os.path.join(output_dir, f"{prefix}.validation.tsv.gz")
    validation.to_csv(tsv_file, sep="\t", header=False, index=False,
        compression="gzip")
    if create_test:
        if prefix is None:
            tsv_file = os.path.join(output_dir, "test.tsv.gz")
        else:
            tsv_file = os.path.join(output_dir, f"{prefix}.test.tsv.gz")
        test.to_csv(tsv_file, sep="\t", header=False, index=False,
            compression="gzip")

if __name__ == "__main__":
    main()