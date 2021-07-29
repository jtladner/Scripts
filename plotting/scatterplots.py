#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


import argparse
#import pandas as pd
#import seaborn as sns
import numpy as np
import inout as io
import itertools as it

# This script is used for generating scatterplots using matplotlib

colPalette = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    reqArgs = parser.add_argument_group('required arguments')
    reqArgs.add_argument('-d', '--data',  help='Data matrix for generating scatterlots.', required=True)

    parser.add_argument('-o', '--outfile', help='Name for out file.')
    parser.add_argument('-x', '--xHead',  help='Header in data file for x-axis.')
    parser.add_argument('-y', '--yHead',  help='Header in data file for y-axis.')

    parser.add_argument('-c', '--color',  help='Header name to use to color points in plot.')
    parser.add_argument('--cMap',  help='Optional way of mapping the "color" variable in the data matrix to another categorical variable. 3 comma-separated variables should be provided: map file, key header, value header.')
    parser.add_argument('--xLog', default=False, type=float, help="Use if you want x-axis to be shown on a log-scale. Argument provided should be a float to add to the y values before calculating the log value.")
    parser.add_argument('--yLog', default=False, type=float, help="Use if you want y-axis to be shown on a log-scale. Argument provided should be a float to add to the y values before calculating the log value.")
    parser.add_argument('--delim', default="\t", help="Delimiter used in the data file.")
    parser.add_argument('--xLab', help="String for x-label. If not provided, --xHead is used.")
    parser.add_argument('--yLab', help="String for y-label. If not provided, --yHead is used.")
#    parser.add_argument('--width', default=5, type=int, help="Figure width.")
#    parser.add_argument('--height', default=4, type=int, help="Figure height.")
    parser.add_argument('--include', help="Header,Value pairs used to indicate a subset of rows to include.", nargs='*')
    parser.add_argument('--exclude', help="Header,Value pairs used to indicate a subset of rows to exclude.", nargs='*')
    parser.add_argument('--xLegend', default=0.1, type=float, help="x-coordinate to use for color legend.")
    parser.add_argument('--yLegend', default='top', help="Indicate whether you want the legend at the 'top' or 'bottom' of the plot.")
    parser.add_argument('--xeqy', default=False, action="store_true", help="Use if you want an x=y line included in the plot.")
    parser.add_argument('--markerSize', default=10, type=int, help="Size of marker used in plot.")
    parser.add_argument('--alpha', default=0.6, type=float, help="Alpha (transparency) value to use in plot.")

    parser.add_argument('-b', '--batch',  help='An alternative way to provide input/output files that allows the generation of multiple plots with a single command. File provided should be tab-delimited, one line per output plot: xHeader, yHeader, outfile, colorHeader. Colorheader column is optional.')
    parser.add_argument('-a', '--allByAll',  help='Optional way to specify xHead and yHead. Should be a list of column headers. A plot will be generated for all pairwise comparisons of the columns in this file. Output names will be generated based on the column name.')

    opts = parser.parse_args()

    
    if opts.cMap:
        mapF, kHead, vHead = opts.cMap.split(",")
        opts.cMap = io.fileDictHeader(mapF, kHead, vHead)
    
    #Read in data file 
    dataD = io.fileDictFull(opts.data, opts.delim)


    if opts.batch:
        with open(opts.batch) as fin:
            for line in fin:

                cols = line.rstrip("\n").split("\t")
                
                opts.xHead = cols[0]
                opts.yHead = cols[1]
                opts.outfile = cols[2]
                if len(cols)==4:
                    opts.color = cols[3]
                
                # Make a subset dict that will be manipulated
                subD = {k:dataD[k] for k in opts.xHead.split(",") + opts.yHead.split(",")}
                
                # Make sure data columns are formatted properly 
                prepData(subD, opts)

                #Generate plot
                scatter(subD, opts)
                
                #Reset Labels
                opts.xLab=None
                opts.yLab=None

    elif opts.allByAll:
        heads = io.fileList(opts.allByAll, header=False)
        for h1, h2 in it.combinations(heads, 2):

                opts.xHead = h1
                opts.yHead = h2
                opts.outfile = "%s_%s.png" % (h1, h2)
                
                # Make a subset dict that will be manipulated
                subD = {k:dataD[k] for k in opts.xHead.split(",") + opts.yHead.split(",")}
                
                # Make sure data columns are formatted properly 
                prepData(subD, opts)

                #Generate plot
                scatter(subD, opts)
                
                #Reset Labels
                opts.xLab=None
                opts.yLab=None

            
    
    else:
        
        if opts.include:
            for each in opts.include:
                header, val = each.split(",")
                dataD = onlyInclude(dataD, header, val)

        if opts.exclude:
            for each in opts.exclude:
                header, val = each.split(",")
                dataD = toExclude(dataD, header, val)

        
        # Make sure data columns are formatted properly 
        prepData(dataD, opts)
    
        #Generate plot
        scatter(dataD, opts)

#----------------------End of main()

def toExclude(dataD, header, val):
    newD = {k:[] for k in dataD}
    for i,v in enumerate(dataD[header]):
        if v!=val:
            for thisK, thisV in dataD.items():
                newD[thisK].append(thisV[i])
    return newD

def onlyInclude(dataD, header, val):
    newD = {k:[] for k in dataD}
    for i,v in enumerate(dataD[header]):
        if v==val:
            for thisK, thisV in dataD.items():
                newD[thisK].append(thisV[i])
    return newD


def prepData(dataD, opts):
    yHeadList = opts.yHead.split(",")
    for each in yHeadList:
        if opts.yLog:
            dataD[each] = [np.log10(float(a)+opts.yLog) if a else None for a in dataD[each]]
        else:
            dataD[each] = [float(a) if a else None for a in dataD[each]]

    xHeadList = opts.xHead.split(",")
    for each in xHeadList:
        if opts.xLog:
            dataD[each] = [np.log10(float(a)+opts.xLog)  if a else None for a in dataD[each]]
        else:
            dataD[each] = [float(a) if a else None for a in dataD[each]]

    #Average values if multiple column headers are provided
    if len(yHeadList) >1: 
        opts.yHead = "|".join([y.split("_")[0] for y in yHeadList])
        dataD[opts.yHead] = [np.mean([dataD[n][i] for n in yHeadList]) for i in range(len(dataD[yHeadList[0]]))]

    if len(xHeadList) >1: 
        opts.xHead = "|".join([x.split("_")[0] for x in xHeadList])
        dataD[opts.xHead] = [np.mean([dataD[n][i] for n in xHeadList]) for i in range(len(dataD[xHeadList[0]]))]

def scatter(dd, opts):
    if opts.xLab==None:
        opts.xLab=opts.xHead
    if opts.yLab==None:
        opts.yLab=opts.yHead

    fig,ax = plt.subplots(1,1,figsize=(5, 4),facecolor='w')    

    if opts.color:
        colD={}
        
        if opts.cMap:
            for i,c in enumerate(sorted(list(set(opts.cMap.values())))):
                colD[c] = colPalette[i]

            sPoints = ax.scatter(dd[opts.xHead], dd[opts.yHead], c=[colD[opts.cMap[a]] for a in dd[opts.color]], alpha=opts.alpha, s=opts.markerSize)

        else:
            for i,c in enumerate(list(set(dd[opts.color]))):
                colD[c] = colPalette[i]
        
            sPoints = ax.scatter(dd[opts.xHead], dd[opts.yHead], c=[colD[a] for a in dd[opts.color]], alpha=opts.alpha, s=opts.markerSize)

        #Make legend
        count=0
        for n, c in colD.items():
            count+=0.08
            if opts.yLegend=="top":
                ax.text(opts.xLegend, 1-count, n, color=c, transform=ax.transAxes, fontsize=15)
            elif opts.yLegend=="bottom":
                ax.text(opts.xLegend, count, n, color=c, transform=ax.transAxes, fontsize=15)

    else:
        ax.scatter(dd[opts.xHead], dd[opts.yHead], alpha=opts.alpha, s=opts.markerSize)
    
    if opts.xeqy:
        thelower = max([min(dd[opts.xHead]), min(dd[opts.yHead])])
        thehigher = min([max(dd[opts.xHead]), max(dd[opts.yHead])])
        ax.plot([thelower, thehigher], [thelower, thehigher], ls=":", c="k")
    
    ax.set_xlabel(opts.xLab, fontsize=20)
    ax.set_ylabel(opts.yLab, fontsize=20)
    ax.tick_params(labelsize=15)
    
    if opts.outfile:
        plt.savefig(opts.outfile,dpi=300,bbox_inches='tight')
    plt.close()

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

