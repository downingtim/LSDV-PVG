#!/usr/bin/env python
import seaborn as sns
import argparse
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd
import sys
import os

def load_data(filename):
    with open(filename, 'r') as f:
        data = [int(x.strip()) for x in f.readlines()]
        n = len(data)
        y = data
        y = [y[i] for i in range(len(y))]
        x = [(i+1) for i in range(n)]
        name = os.path.basename(filename)
        if os.path.dirname(filename) != "":
            name = os.path.dirname(filename)+"/"+name
        df = pd.DataFrame({"x":x, "kmers":y, "species": name})
        return df

def main():
    parser = argparse.ArgumentParser(description="Plot an histogram from a single histogram file.")
    parser.add_argument("input_file", help="Input histogram file.")
    parser.add_argument("output_file", help="Output figure file in pdf.")

    args = parser.parse_args()

    sns.set_context("paper")
    sns.set_style("whitegrid")
    _, ax = plt.subplots(figsize=(6, 3))

    filename = args.input_file
    df = load_data(filename)

    sns.barplot(x='x', y='kmers', hue='species', data=df)
    sns.set_style("whitegrid")

    plt.gca().spines["top"].set_visible(True)
    plt.gca().spines["right"].set_visible(True)
    plt.gca().spines["bottom"].set_visible(True)
    plt.gca().spines["left"].set_visible(True)
    plt.gca().spines["left"].set_visible(True)
    line_width = 0.28
    plt.gca().spines["left"].set_linewidth(line_width)
    plt.gca().spines["right"].set_linewidth(line_width)
    plt.gca().spines["bottom"].set_linewidth(line_width)
    plt.gca().spines["top"].set_linewidth(line_width)
    plt.grid(True, linewidth=line_width)

    n = len(df.index)
    if n > 25:
        xticklabels = [i for i in range(1, n, n//20)] +[n]
        xticks = [x -1 for x in xticklabels]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)

    plt.xticks(fontsize=8)
    plt.xlabel('Number of Genomes') 
    plt.ylabel('Number of $k$-mers')

    # Legend
    plt.tight_layout()
    font = FontProperties()
    font.set_style('italic')
    plt.legend(loc="upper left",prop=font)

    root, ext = os.path.splitext(args.output_file)
    if ext.lower() != '.pdf':
        args.output_file = args.output_file + '.pdf'
    plt.savefig(args.output_file, bbox_inches='tight')

if __name__ == "__main__":
    main()
