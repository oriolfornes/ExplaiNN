#!/usr/bin/env python

import click
import os
import shutil
import torch

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.argument(
    "model_file",
    type=click.Path(exists=True, resolve_path=True)
)

def fixit(model_file):

    # Initialize
    dir_name = os.path.dirname(model_file)
    file_name = os.path.basename(model_file)
    tmp_file = os.path.join(dir_name, ".%s" % file_name)

    # Copy it!
    if not os.path.exists(tmp_file):
        shutil.copy(model_file, tmp_file)

    # Load it!
    selene_dict = torch.load(model_file)

    # Fix it!
    if selene_dict["options"]["apply_sigmoid"]:
        selene_dict["options"].setdefault("input_data", "binary")
    else:
        selene_dict["options"].setdefault("input_data", "linear")
    del selene_dict["options"]["apply_sigmoid"]

    # Save it!
    torch.save(selene_dict, model_file)

if __name__ == "__main__":
    fixit()