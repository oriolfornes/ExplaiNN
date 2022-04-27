#!/usr/bin/env python

import click
import copy
import json
import os
import pandas as pd
import random
import re
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder

# Globals
scripts_dir = os.path.dirname(os.path.realpath(__file__))

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "json_file", type=click.Path(exists=True, resolve_path=True)
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)
@click.option(
    "-p", "--prefix",
    help="Output prefix.",
)
@click.option(
    "-s", "--shuffle-lowercase",
    help="Shuffle lowercase letters.",
    is_flag=True
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

    # Get TSV files for CAMs
    _to_CAM(args["json_file"], args["output_dir"], args["prefix"],
        args["shuffle_lowercase"], args["test"])

def _to_CAM(json_file, output_dir="./", prefix=None, shuffle_lc=False,
    create_test=False):

    # Initialize
    regexp = re.compile(r"[a-z]+")

    # Load JSON
    with open(json_file) as handle:
        sequences = json.load(handle)
    sequences.pop(0)

    # Ys
    enc = OneHotEncoder()
    arr = [[i] for i in range(len(sequences[0]) - 1)]
    if len(arr) == 2:
        ys = arr
    else:
        enc.fit(arr)
        ys = enc.transform(arr).toarray().tolist()

    # Get DataFrame
    data = []
    for i in range(len(sequences)):
        for j in range(1, len(sequences[i])):
            seq_id = sequences[i][j][0]
            s = sequences[i][j][1]
            if shuffle_lc:
                # 1) extract blocks of lowercase letters;
                # 2) shuffle the letters; and
                # 3) put the shuffled letters back
                l = list(s)
                for m in re.finditer(regexp, s):
                    sublist = l[m.start():m.end()]
                    random.shuffle(sublist)
                    l[m.start():m.end()] = copy.copy(sublist)
                s = "".join(l)
            data.append([seq_id, s.upper()] + ys[j - 1])
    df = pd.DataFrame(data)
    df = df.groupby(1).max().reset_index()
    df = df.reindex(sorted(df.columns), axis=1)

    # Get data splits
    train, validation = train_test_split(df, test_size=.2)
    if create_test:
        validation, test = train_test_split(validation, test_size=.5)

    # Save sequences
    if prefix is None:
        tsv_file = os.path.join(output_dir, "train.tsv.gz")
    else:
        tsv_file = os.path.join(output_dir, "%s.train.tsv.gz" % prefix)
    train.to_csv(tsv_file, sep="\t", header=False, index=False,
        compression="gzip")
    if prefix is None:
        tsv_file = os.path.join(output_dir, "validation.tsv.gz")
    else:
        tsv_file = os.path.join(output_dir, "%s.validation.tsv.gz" % prefix)
    validation.to_csv(tsv_file, sep="\t", header=False, index=False,
        compression="gzip")
    if create_test:
        if prefix is None:
            tsv_file = os.path.join(output_dir, "test.tsv.gz")
        else:
            tsv_file = os.path.join(output_dir, "%s.test.tsv.gz" % prefix)
        test.to_csv(tsv_file, sep="\t", header=False, index=False,
            compression="gzip")

if __name__ == "__main__":
    main()