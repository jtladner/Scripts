#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


import argparse
import pandas as pd
import seaborn as sns
import numpy as np
import inout as io
import matrixtools as mt
#import statistics as stat
#import itertools as it
#from collections import defaultdict



# This script looks for peptides that show correlated reactivity across samples

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    # Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
    reqArgs.add_argument('-d', '--data',  help='Data matrix for with scores that will be used as the basis for correaltions.', required=True)
    reqArgs.add_argument('-o', '--out', help='Name for output file.', required=True)
#    reqArgs.add_argument('-m', '--meta',  help='Header in data file for x-axis categories.', required=True)
#    reqArgs.add_argument('-y', '--yHead',  help='Header in data file for y-axis.', required=True)

    parser.add_argument('-m', '--meta', help='Optional metadata file that can be used for specifying a subset of peptides to consider.')

#    parser.add_argument('--logY', default=False, action="store_true", help="Use this option if you want the y-axis to be coversted to log scale. If some of the values are <=0, an integer will be added to all y-axis values prior to conversion")
#    parser.add_argument('--delim', default="\t", help="Delimiter used in the data file.")
#    parser.add_argument('--xRotate', default=False, action="store_true", help="Use this option if you want the x-labels to be rotated vertically.")
#    parser.add_argument('--width', default=5, type=int, help="Figure width.")
#    parser.add_argument('--height', default=4, type=int, help="Figure height.")
#    parser.add_argument('--include', help="Header,Value pairs used to indicate a subset of rows to include.", nargs='*')

    opts = parser.parse_args()

    pyMat = mt.readMatrix(opts.data)
    df = pd.DataFrame(pyMat,columns=sorted(list(pyMat.keys())))
    corrMatrix = df.T.corr()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

