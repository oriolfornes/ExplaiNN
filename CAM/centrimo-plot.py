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
    parser.add_argument("plot_file", metavar="plot.png")

    return(parser.parse_args())

def get_figure(centrimo_file, counts_file):

    # Read p-value
    with open(centrimo_file, "r") as f:
        for line in f:
            if line.startswith("   1"):
                line = line.split("\t")
                p = Decimal(line[5])
                break

    # Read counts
    df = pd.read_csv(
        counts_file,
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
    ax.set_yticklabels(labels=[0.0, 0.25, 0.5, 0.75, 1.0], fontproperties=prop)
    ax.text(-450, .9, str(p))

    return(fig)

def main():

    # Parse arguments
    args = parse_args()

    # Get figure
    fig = get_figure(args.centrimo_file, args.counts_file)

    # Save logo
    fig.savefig(args.plot_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()