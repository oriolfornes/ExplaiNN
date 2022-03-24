#!/usr/bin/env python

import argparse
from Bio import motifs
import logomaker
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# Specify font
scripts_dir = os.path.dirname(os.path.realpath(__file__))
prop = fm.FontProperties(fname=os.path.join(scripts_dir, "fonts", "Arial.ttf"))

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    Parses arguments provided through the command line.
    """

    # Initialize
    parser = argparse.ArgumentParser()

    # Mandatory args
    parser.add_argument("motif_file", metavar="motif.jaspar")
    parser.add_argument("logo_file", metavar="logo.png")

    # Optional args
    parser.add_argument("-r", "--rev-complement", action="store_true",
        help="plot the reverse complement logo")

    return(parser.parse_args())

def get_figure(motif_file, rc=False):

    # From https://biopython.readthedocs.io/en/latest/chapter_motifs.html
    m = motifs.read(open(motif_file), "jaspar")
    pwm = list(m.counts.normalize(pseudocounts=.5).values())

    return(_get_figure(pwm, rc))

def _get_figure(pwm, rc=False):

    # From https://www.bioconductor.org/packages/release/bioc/html/seqLogo.html
    if rc:
        arr = np.array(pwm)
        pwm = np.flip(arr).tolist()
    IC = 2 + np.add.reduce(pwm * np.log2(pwm))
    df = pd.DataFrame({
        "pos": [i + 1 for i in range(len(IC))],
        "A": pwm[0] * IC,
        "C": pwm[1] * IC,
        "G": pwm[2] * IC,
        "T": pwm[3] * IC
    })
    df = df.set_index("pos")

    # From https://logomaker.readthedocs.io/en/latest/examples.html
    fig, ax = plt.subplots(1, 1, figsize=(len(df)/2.0, 2))
    logo = logomaker.Logo(df, ax=ax, show_spines=False)
    logo.style_spines(spines=["left", "bottom"], visible=True)
    logo.ax.set_aspect(1.5)
    logo.ax.xaxis.set_ticks(list(df.index))
    logo.ax.set_xticklabels(labels=list(df.index), fontproperties=prop)
    logo.ax.set_ylabel("Bits", fontproperties=prop)
    logo.ax.set_ylim(0, 2)
    logo.ax.yaxis.set_ticks([0, 1, 2])
    logo.ax.set_yticklabels(labels=[0, 1, 2], fontproperties=prop)

    return(fig)

def main():

    # Parse arguments
    args = parse_args()

    # Get figure
    fig = get_figure(args.motif_file, rc=args.rev_complement)

    # Save
    fig.savefig(args.logo_file, bbox_inches="tight", pad_inches=0)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()