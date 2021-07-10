#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import click
import gzip
import os
from sklearn.model_selection import train_test_split

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "sequences_dir", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True
)

def main(**params):

    # Create output dir
    if not os.path.exists(params["output_dir"]):
        os.makedirs(params["output_dir"])

    # Initialize
    train = []
    validation = []
    test = []

    # For species...
    for s in ["human", "mouse"]:

        # Initialize
        sequences = []

        # Sequences
        f = os.path.join(params["sequences_dir"], s, f"{s}.fa")
        for r in SeqIO.parse(open(f), "fasta"):
            if r.seq.count("N") > 0: continue
            if s == "human":
                sequences.append(SeqRecord(r.seq, id=r.id, description="1.;0."))
            else:
                sequences.append(SeqRecord(r.seq, id=r.id, description="0.;1."))

        # Split records 80/20
        t, rest = train_test_split(sequences, test_size=.2, random_state=1)
        train.append(t)

        # Split the 20 split 50/50
        v, t = train_test_split(rest, test_size=.5, random_state=1)
        validation.append(v)
        test.append(t)

    # Save sequences
    f = os.path.join(params["output_dir"], "train.fa.gz")
    with gzip.open(f, "wt") as handle:  
        sequences = [j for i in zip(train[0], train[1][:len(train[0])]) for j in i]
        SeqIO.write(sequences, handle, "fasta")
    f = os.path.join(params["output_dir"], "validation.fa.gz")
    with gzip.open(f, "wt") as handle:
        sequences = [j for i in zip(validation[0], validation[1][:len(validation[0])]) for j in i]
        SeqIO.write(sequences, handle, "fasta")
    f = os.path.join(params["output_dir"], "test.fa.gz")
    with gzip.open(f, "wt") as handle:
        sequences = [j for i in zip(test[0], test[1][:len(train[0])]) for j in i]
        SeqIO.write(sequences, handle, "fasta")


if __name__ == "__main__":
    main()