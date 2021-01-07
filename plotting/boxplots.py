#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


import optparse
import pandas as pd
import seaborn as sns
import numpy as np
import inout as io
#import statistics as stat
#import itertools as it
#from collections import defaultdict



# This script reads in 1) normalized or raw read counts and 2) bins and calculated zscores

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-d', '--data',  help='Data matrix for generating boxplots. [None, REQ]')
    p.add_option('-o', '--outfile', help='Name for out file [None, REQ]')
    p.add_option('-x', '--xHead',  help='Header in data file for x-axis categories. [None, REQ]')
    p.add_option('-y', '--yHead',  help='Header in data file for y-axis. [None, REQ]')
    p.add_option('--hue',  help='Header in data file for "hue", or subcategories for x-axis. [None]')
    p.add_option('--logY', default=False, action="store_true", help="Use this option if you want the y-axis to be coversted to log scale. If some of the values are <=0, an integer will be added to all y-axis values prior to conversion [None]")


    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\t]")
    p.add_option('--xLab', help="String for x-label. If not provided, --xHead is used. [None]")
    p.add_option('--yLab', help="String for y-label. If not provided, --yHead is used.[None]")
    p.add_option('--xRotate', default=False, action="store_true", help="Use this option if you want the x-labels to be rotated vertically. [None]")
    p.add_option('--width', default=5, type="int", help="Figure width. [5]")
    p.add_option('--height', default=4, type="int", help="Figure height. [4]")

    opts, args = p.parse_args()
    
    #Read in data file and create pandas dataframe
    dataD = io.fileDictFull(opts.data, opts.delim)
    dataD[opts.yHead] = [float(a) for a in dataD[opts.yHead]]
    
    if opts.logY:
        opts.addInt=False
        minY = min(dataD[opts.yHead])
        if minY<=0:
            opts.addInt = 1
            while opts.addInt + minY <=0:
                opts.addInt+=1
            dataD[opts.yHead] = [a+opts.addInt for a in dataD[opts.yHead]]
        dataD[opts.yHead] = [np.log10(a) for a in dataD[opts.yHead]]
    
    df = pd.DataFrame(dataD)
    
    catBoxplot(df, opts, colorHead=opts.hue, xLab=opts.xLab, yLab=opts.yLab, out=opts.outfile)

#----------------------End of main()

def catBoxplot(df, opts, colorHead=None, xLab=None, yLab=None, out=None):
    if xLab==None:
        xLab=opts.xHead
    if yLab==None:
        yLab=opts.yHead
        
    fig,ax = plt.subplots(1,1,figsize=(opts.width, opts.height),facecolor='w')
    
    if colorHead:
        sns.boxplot(x=opts.xHead, y=opts.yHead, hue=colorHead, data=df, ax=ax, fliersize=0)
        sns.stripplot(x=opts.xHead, y=opts.yHead, hue=colorHead, data=df, jitter=True, dodge=True, linewidth=0.5, ax=ax)
    
    else:
        sns.boxplot(x=opts.xHead, y=opts.yHead, data=df, ax=ax, fliersize=0)
        sns.stripplot(x=opts.xHead, y=opts.yHead, data=df, jitter=True, dodge=True, linewidth=0.5, ax=ax)

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    ax.set_xlabel(xLab, fontsize=20)
    if opts.logY:
        if opts.addInt:
            ax.set_ylabel("log10(%s+%d)" % (yLab, opts.addInt), fontsize=20)
        else:
            ax.set_ylabel("log10(%s)" % yLab, fontsize=20)
    else:
        ax.set_ylabel(yLab, fontsize=20)
    ax.tick_params(labelsize=15)
    
    if opts.xRotate:
        plt.xticks(rotation="vertical")
    
    if out:
        plt.savefig(out,dpi=300,bbox_inches='tight')


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

