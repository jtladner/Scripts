#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.patches as patches

import argparse
import pandas as pd
#import seaborn as sns
import numpy as np
import inout as io
import itertools as it

# This script is used for generating scatterplots using matplotlib

colPalette = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument('-e', '--epiInfo',  help='Tab-delimited file containing information about each set of epitopes for which you want to generate a map. There should be one row per map. Expected headers: "ProteinName", "ProteinLength", "EpiCoordsFilepath", "EpiColor"', required=True)

    parser.add_argument('--figWidth', default=15, type=int, help="Figure height.")
    parser.add_argument('--figHeight', default=2, type=int, help="Figure height.")
    parser.add_argument('--mapHeight', default=10, type=int, help="Figure height.")

    opts = parser.parse_args()
    
    # Read in epi group info
    df = pd.read_csv(opts.epiInfo, sep="\t", header=0)
    
    # Step through each epitope group
    for i, row in df.iterrows():
        
        # Initiate figure and axis objects
        fig,ax = plt.subplots(1, 1, figsize=(opts.figWidth,opts.figHeight),facecolor='w')

        thisDF = pd.read_csv(row["EpiCoordsFilepath"], sep="\t", header=0)

        # Create a Rectangle patch representing the full protein
        rect = patches.Rectangle((0, 0), row["ProteinLength"], opts.mapHeight, linewidth=3, edgecolor='k', facecolor='none')
        ax.add_patch(rect)

        # Create a Rectangle patch for each epitope
        for j, thisRow in thisDF.iterrows():
            thisL = thisRow["End"]-thisRow["Start"]+1
            rect = patches.Rectangle((thisRow["Start"], 0), thisL, opts.mapHeight, linewidth=1, edgecolor='k', facecolor=row["EpiColor"], alpha=0.7)
            ax.add_patch(rect)

        for each in [0, row["ProteinLength"]]:
            ax.text(each, -3, str(f"{each:3d}"), rotation=45, ha="right", va="center")

        # Modify plot limits
        ax.set_xlim(-5, row["ProteinLength"]+5)
        ax.set_ylim(-1, opts.mapHeight+1)

        # Hide axes
        ax.axis("off")

        # ylabs = [0.05, 0.01, 0.001]
        # ax.set_yticks([-np.log10(each) for each in ylabs])
        # ax.set_yticklabels(ylabs)

        # ax.set_xlabel("Zscore difference", fontsize=30)
        # ax.set_ylabel("T-test p-value", fontsize=30)

        # ax.axhline(-np.log10(0.05), c="k", ls=":")
        # ax.axvline(0, c="k", ls=":")
        
        outBase = f'{row["ProteinName"]}_epiMap'
        fig.savefig(f'{outBase}.pdf',dpi=300,bbox_inches='tight')
        fig.savefig(f'{outBase}.png',dpi=300,bbox_inches='tight')
#----------------------End of main()

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

