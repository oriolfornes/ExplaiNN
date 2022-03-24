#!/usr/bin/env python

import argparse
import os
import re

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
    parser.add_argument("meme_file", metavar="file.meme")
    parser.add_argument(
        "--out-dir", default="./",
        help="output directory (default: ./)"
    )
    parser.add_argument(
        "--prefix", default="motif",
        help="prefix for LPM file (default: \"motif\")"
    )

    return(parser.parse_args())

def meme_to_lpm(meme_to_lpm, out_dir="./", prefix="motif"):

    # Initialize
    m = []
    p = []

    # Create output dirs
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # MOTIF 1-CCASYAGRKGGCRSY STREME-1
    # letter-probability matrix: alength= 4 w= 15 nsites= 52413 S= 2.7e-9718
    # MOTIF AICE.filter0 ENCODE.IRF4-AICE CCTCCTTTGCAGCCGGGATGCAGTTC
    # letter-probability matrix: alength= 4 w= 26 nsites= 54 E= 0
    # >letter-probability matrix filter1 TF-Binding-Matrix.CTCF: alength= 4 w= 19 nsites= 577 E= 0
    with open(meme_to_lpm) as handle:
        for line in handle:
            line = line.strip("\n")
            if line.startswith("MOTIF"):
                m.append(line)
                p.append([])
            elif line.startswith("letter-probability matrix:"):
                m[-1] = f">letter-probability matrix {m[-1]}: {line[27:]}"
            else:
                n = re.search("^\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$", line)
                if n:
                    p[-1].append("%s\t%s\t%s\t%s" % (n.group(1).rjust(10),
                                                     n.group(2).rjust(10),
                                                     n.group(3).rjust(10),
                                                     n.group(4).rjust(10)))

    for i, z in enumerate(zip(m, p)):
        motif_file = os.path.join(out_dir, f"{prefix}{i}.lpm")
        with open(motif_file, "wt") as handle:
            handle.write("%s\n" % z[0])
            handle.write("%s\n" % "\n".join(z[1]))

def main():

    # Parse arguments
    args = parse_args()

    meme_to_lpm(args.meme_file, args.out_dir, args.prefix)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()