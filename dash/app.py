#!/usr/bin/env python

import argparse
from Bio import motifs
import copy
import dash
import dash_html_components as html
import os
import numpy as np
import pandas as pd
import re
import shutil
import sys

# Specify imports
scripts_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, scripts_dir)
from sortByCentrimoPvalue import get_sorted_motifs

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
    parser.add_argument("results_dir", metavar="results-dir")

    return(parser.parse_args())

def get_data(results_dir):

    # Initialize
    data = []

    # For each subdir...
    for subdir in sorted(os.listdir(results_dir)):

        # Initialize
        base_dir = os.path.join(results_dir, subdir)
        motifs_dir = os.path.join(base_dir, "motifs", "jaspar", "pfm")
        logos_dir = os.path.join(base_dir, "motifs", "logos")
        centrimo_dirs = os.path.join(base_dir, "centrimo")

        # Get TF, experiment, peaks
        m = re.search("^(\w+)\W(.+)\W(PEAKS\d+)$", subdir)
        tf = m.group(1)
        experiment = m.group(2)
        peaks = m.group(3)

        # Row:
        # 1) TF
        # 2) ExperimentId
        # 3) PeaksId
        # 4) Sites
        # 5) PWM
        # 6) Consensus
        # 7) Logo
        # 8) Centrality
        # 9) Plot
        row = [tf, experiment, peaks, None, None, None, None, None, None]

        # Read motifs "assets"
        extension = {}
        motif_file = os.path.join(motifs_dir, "%s_matrix_list.tab" % subdir)
        with open(motif_file, "r") as f:
            for line in f:
                if line.startswith(";"):
                    continue
                line = line.strip().split("\t")
                extension.setdefault(line[-1], line[1])

        # i.e. non-resolved TF
        if len(extension) == 0:
            data.append(row)
            continue

        # For each motif, pvalue...
        for motif, pvalue in get_sorted_motifs(centrimo_dirs):

            # Initialize
            motif_files = [
                os.path.join(
                    motifs_dir, "%s_%s.jaspar" % (subdir, extension[motif])
                ),
                os.path.join(
                    "assets", "%s_%s.jaspar" % (subdir, extension[motif])
                )
            ]
            logo_files = [
                os.path.join(
                    logos_dir, "%s_%s.png" % (subdir, extension[motif])
                ),
                os.path.join(
                    "assets", "%s_%s.logo.png" % (subdir, extension[motif])
                )
            ]
            centrimo_files = [
                os.path.join(
                    centrimo_dirs,
                    "%s_%s.501bp" % (subdir, extension[motif][12:]),
                    "plot.png"
                ),
                os.path.join(
                    "assets", "%s_%s.centrimo.png" % (subdir, extension[motif])
                )
            ]

            # https://biopython.readthedocs.io/en/latest/chapter_motifs.html
            m = motifs.read(open(motif_files[0]), "jaspar")
            m.pseudocounts = {"A": 0.5, "C": 0.5, "G": 0.5, "T": 0.5}

            # Row:
            motif_row = copy.copy(row)
            motif_row[3] = int(sum([m.counts[n][0] for n in "ACGT"]))
            motif_row[4] = motif_files[1]
            shutil.copy(motif_files[0], motif_files[1])
            motif_row[5] = str(m.pssm.consensus)
            motif_row[6] = logo_files[1]
            shutil.copy(logo_files[0], logo_files[1])
            motif_row[7] = str(pvalue)
            motif_row[8] = centrimo_files[1]
            shutil.copy(centrimo_files[0], centrimo_files[1])
            data.append(motif_row)

    return(np.array(data))

def get_html_table_with_hyperlinks(df):
    """
    From: https://github.com/plotly/dash-recipes/
    Script: dash-html-table-hyperlinks.py
    """

    # Initialize
    rows = [html.Tr([html.Th(col) for col in df.columns])]
    # style = {"height": "25%", "width": "25%"}
    

    for i in range(len(df)):
        row = []
        for col in df.columns:
            value = df.iloc[i][col]
            if os.path.isfile(str(value)):
                file_name = os.path.basename(value)
                # if col == "Logo":
                #     cell = html.Td(
                #         html.Img(src=get_src(file_name), style=style)
                #     )
                # else:
                cell = html.Td(
                    html.A(href=get_src(file_name), children=col)
                )
            else:
                cell = html.Td(children=value)
            row.append(cell)
        rows.append(html.Tr(row))

    return html.Table(rows)

def get_src(source):

    return(app.get_asset_url(source))

# Parse arguments
args = parse_args()

# Create assets dir
if not os.path.isdir("assets"):
    os.makedirs("assets")

# Get data frame
cols = [
    "TF", "ExperimentId", "PeaksId", "Sites", "PWM", "Consensus", "Logo",
    "Centrality", "Plot"
]
df = pd.DataFrame(get_data(args.results_dir), columns=cols)

# Use plotly example stylesheet
external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

# Initialize app
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Create html table with hyperlinks
app.layout = html.Div(
    children=[get_html_table_with_hyperlinks(df)]
)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    app.run_server(host="127.0.0.1")