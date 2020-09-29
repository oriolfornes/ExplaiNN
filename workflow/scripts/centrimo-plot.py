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

    # Get centrimo plot

    # Save logo
    fig.savefig(args.pdf_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()