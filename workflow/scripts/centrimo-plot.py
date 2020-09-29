#!/usr/bin/env python

import argparse
from decimal import Decimal
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# Specify font
scripts_dir = os.path.dirname(os.path.realpath(__file__))
prop = fm.FontProperties(fname=os.path.join(scripts_dir, "fonts", "arial.ttf"))

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
    parser.add_argument("centrimo_file", metavar="centrimo.tsv")
    parser.add_argument("counts_file", metavar="site_counts.txt")
    parser.add_argument("pdf_file", metavar="plot.pdf")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Read p-value
    with open(args.centrimo_file, "r") as f:
        for line in f:
            if line.startswith("   1"):
                line = line.split("\t")
                p = Decimal(line[5])
                break

    # Read counts
    df = pd.read_csv(
        args.counts_file,
        # names=["Distance to peak centre", "Number of motif occurrences"],
        names=["x", "y"],
        sep="\t",
        skiprows=1
    )

    # Plot centrality
    fig, ax = plt.subplots()
    ax.plot(list(df["x"]), list(df["y"]/df["y"].max()))
    ax.set_xlabel("Distance to peak centre", fontproperties=prop)
    ax.set_xlim(-500, 500)
    ax.xaxis.set_ticks([-500, -250, 0, 250, 500])
    ax.set_xticklabels(labels=[-500, -250, 0, 250, 500], fontproperties=prop)
    ax.set_ylabel("Relative number of motif occurrences", fontproperties=prop)
    ax.set_ylim(0, 1)
    ax.yaxis.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(labels=[0.0, 0.25, 0.5, 0.75, 1.0], fontproperties=prop)
    ax.text(-450, .9, str(p))

    # Save logo
    fig.savefig(args.pdf_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()