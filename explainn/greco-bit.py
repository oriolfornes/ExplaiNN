#!/usr/bin/env python

import click
from click_option_group import optgroup
import os
import re
import subprocess as sp

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "training_dir",
    type=click.Path(exists=True, resolve_path=True),
)
@click.argument(
    "validation_dir",
    type=click.Path(exists=True, resolve_path=True),
)
@click.option(
    "-o", "--output-dir",
    help="Output directory.",
    type=click.Path(resolve_path=True),
    default="./",
    show_default=True,
)

def main(**params):

    # Initialize
    files = {}
    scripts_dir = os.path.dirname(os.path.realpath(__file__))

    # Get files
    for f in os.listdir(params["training_dir"]):
        m = re.search("^(\S+)@\S+.fa.gz$", f)
        files.setdefault(m.group(1), [None, None])
        files[m.group(1)][0] = os.path.join(params["training_dir"], f)
    for f in os.listdir(params["validation_dir"]):
        m = re.search("^(\S+)@\S+.fa.gz$", f)
        files.setdefault(m.group(1), [None, None])
        files[m.group(1)][1] = os.path.join(params["validation_dir"], f)

    # Train
    for tf in sorted(files):
        m = re.search("^(\S+@\S+).fa.gz$", os.path.split(files[tf][0])[1])
        o = os.path.join(params["output_dir"], m.group(1))
        best_model = os.path.join(o, "best_model.pth.tar")
        if not os.path.exists(best_model):
            cmd = f"{scripts_dir}/train.py -i 10 -o {o} -r {files[tf][0]} {files[tf][1]}"
            _ = sp.run([cmd], shell=True)

if __name__ == "__main__":
    main()