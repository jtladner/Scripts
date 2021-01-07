#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


import optparse
#import pandas as pd
#import seaborn as sns
import numpy as np
import inout as io

# This script is used for generating scatterplots using matplotlib

colPalette = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"]

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-d', '--data',  help='Matrix containing data that will be used to generate scatterplots. [None, REQ]')
    p.add_option('-o', '--outfile', help='Name for out file [None, REQ]')
    p.add_option('-x', '--xHead', help='Header name to use for x axis. [None, OPT]')
    p.add_option('-y', '--yHead', help='Header name to use for y axis. [None, OPT]')
    p.add_option('-c', '--color',  help='Header name to use to color points in plot. [None, OPT]')
    p.add_option('--cMap',  help='Optional way of mapping the "color" variable in the data matrix to another categorical variable. 3 comma-separated variables should be provided: map file, key header, value header [None, OPT]')
    p.add_option('-b', '--batch',  help='An alternative way to provide input/output files that allows the generation of multiple plots with a single command. File provided should be tab-delimited, one line per output plot: xHeader, yHeader, outfile, colorHeader. Colorheader column is optional. [None, OPT]')

    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\t]")
    p.add_option('--xLab', help="String for x-label. If not provided, --xHead is used. [None]")
    p.add_option('--yLab', help="String for y-label. If not provided, --yHead is used.[None]")
    p.add_option('--xLegend', default=0.1, type="float", help="x-coordinate to use for color legend.[0.1]")
    p.add_option('--yLegend', default='top', help="Indicate whether you want the legend at the 'top' or 'bottom' of the plot.[top]")

    p.add_option('--xLog', default=False, action="store_true", help="Use if you want x-axis to be shown on a log-scale. [False]")
    p.add_option('--yLog', default=False, action="store_true",  help="Use if you want y-axis to be shown on a log-scale.[False]")

    opts, args = p.parse_args()
    
    if opts.cMap:
        mapF, kHead, vHead = opts.cMap.split(",")
        opts.cMap = io.fileDictHeader(mapF, kHead, vHead)
    
    if opts.batch:
        with open(opts.batch) as fin:
            for line in fin:

                #Read in data file 
                dataD = io.fileDictFull(opts.data, opts.delim)

                cols = line.rstrip("\n").split("\t")
                opts.xHead = cols[0]
                opts.yHead = cols[1]
                opts.outfile = cols[2]
                if len(cols)==4:
                    opts.color = cols[3]
                    
                # Make sure data columns are formatted properly 
                prepData(dataD, opts)

                #Generate plot
                scatter(dataD, opts)
                
                #Reset Labels
                opts.xLab=None
                opts.yLab=None

    
    else:
        #Read in data file 
        dataD = io.fileDictFull(opts.data, opts.delim)

        # Make sure data columns are formatted properly 
        prepData(dataD, opts)
    
        #Generate plot
        scatter(dataD, opts)

#----------------------End of main()

def prepData(dataD, opts):
    yHeadList = opts.yHead.split(",")
    for each in yHeadList:
        if opts.yLog:
            dataD[each] = [np.log10(float(a)) for a in dataD[each]]
        else:
            dataD[each] = [float(a) for a in dataD[each]]

    xHeadList = opts.xHead.split(",")
    for each in xHeadList:
        if opts.xLog:
            dataD[each] = [np.log10(float(a)) for a in dataD[each]]
        else:
            dataD[each] = [float(a) for a in dataD[each]]

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

            sPoints = ax.scatter(dd[opts.xHead], dd[opts.yHead], c=[colD[opts.cMap[a]] for a in dd[opts.color]], alpha=0.6)

        else:
            for i,c in enumerate(list(set(dd[opts.color]))):
                colD[c] = colPalette[i]
        
            sPoints = ax.scatter(dd[opts.xHead], dd[opts.yHead], c=[colD[a] for a in dd[opts.color]], alpha=0.6)

        #Make legend
        count=0
        for n, c in colD.items():
            count+=0.08
            if opts.yLegend=="top":
                ax.text(opts.xLegend, 1-count, n, color=c, transform=ax.transAxes, fontsize=15)
            elif opts.yLegend=="bottom":
                ax.text(opts.xLegend, count, n, color=c, transform=ax.transAxes, fontsize=15)

    else:
        ax.scatter(dd[opts.xHead], dd[opts.yHead], alpha=0.6)
    
    
    ax.set_xlabel(opts.xLab, fontsize=20)
    ax.set_ylabel(opts.yLab, fontsize=20)
    ax.tick_params(labelsize=15)
    
    if opts.outfile:
        plt.savefig(opts.outfile,dpi=300,bbox_inches='tight')
    plt.close()

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

