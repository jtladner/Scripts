#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42


import optparse, glob, os, sys
#import itertools as it
import numpy as np
from scipy.stats import linregress
#import inout as io             #Available at https://github.com/jtladner/Modules

from collections import defaultdict

# This script is used for processing data from competitive FLISA assays with multiple dilutions per sample

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    #Experimental data
    p.add_option('-d', '--data',  help='File containing data from an assay. Should include both a positive and negative control sample. Columns: sample_name, rep#, sample_type, datapoints, associated_blank [None, REQ]')
    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\\t]")

    #Expectation info for flagging potential issues
#    p.add_option('-e', '--enrDir',  help='Directory containing lists of enriched peptides. [None, REQ]')
#    p.add_option('--enrExt', help='Common file ending for files containing enriched sets of peptides. Once removed, the remaining filename should consist only of sample name(s) [None, REQ]')
#    p.add_option('--snDelim', default="~", help='Delimiter used to separate sample names in pEnrich files [~]')

    #Negative contols
#    p.add_option('--negMatrix',  help='Optional way to provide a separate data matrix for negative controls. If not provided, negative controls will be assumed to be from the same matrix as the experiemntal data  [None, OPT]')
#    p.add_option('-c', '--negControls',  help='Comma-delimited names of samples to use as negative controls. If this is not provided, then all of the samples in the negative control matrix will be used [None, REQ]')
#    p.add_option('-x', '--xLabel', default="", help='Label to be used for x-axis. []')

    #Output controls
#    p.add_option('-o', '--outDir', help='Directory name for output files. Will be created. [None]')
#    p.add_option('--plotType', default="png", help='Type of plot to generate. This should be a file extension recognized by matplotlib, e.g., pdf, png, tiff [png]')
#    p.add_option('--plotLog', default=False, action="store_true", help="Use if you want axes to be shown on a log-scale. [False]")
#    p.add_option('--negColor', default="#1b9e77", help='Color to use for Unenriched peptides. [#1b9e77]')
#    p.add_option('--posColor', default="#d95f02", help='Color to use for Enenriched peptides. [#d95f02]')

    opts, args = p.parse_args()

    # Print command used 
    print("Command run: '%s'" % (sys.argv[0]))
    
    #Read in data file 
    dataD, dataLabels, emptyLabel = parseData(opts.data, delim=opts.delim)
    
    abort = 0
    
    if len(dataD["positive"]) != 1:
        print("Expected 1 positive control sample, found %d" % (len(dataD["positive"])))
        abort = 1
    if len(dataD["negative"]) != 1:
        print("Expected 1 negative control sample, found %d" % (len(dataD["negative"])))
        abort = 1
    if len(dataD["experimental"]) == 0:
        print("Expected at least 1 experimental sample, found %d" % (len(dataD["experimental"])))
        abort = 1
    
    if abort:
        print("Aborting run because of the above error(s).")
    else:
            
        
        #Calculate slope, etc. for positive control sample
        posName = list(dataD["positive"].keys())[0]
        posSlope, posPvalue, posRvalue, posDatapoints = bestSlope(dataD["positive"][posName], dataLabels)
        print("Positive", posName, posSlope, posPvalue, posRvalue, posDatapoints)
        
        
        #Calculate slope, etc. for negative control sample
        negName = list(dataD["negative"].keys())[0]
        negSlope, negPvalue, negRvalue, negDatapoints = bestSlope(dataD["negative"][negName], dataLabels)
        print("Negative", negName, negSlope, negPvalue, negRvalue, negDatapoints)

    with open("%s_results.tsv" % (".".join(opts.data.split(".")[:-1])), "w") as fout:
        fout.write("Sample\tCategory\tDatapointsUsed\tSlope\tPvalue\tCorrelationCoefficient\tPosSlopeRatio\tNegSlopeRatio\n")
        slopeD={}
        
        for samp, sampD in dataD["experimental"].items():
            slopeD[samp] = {}
            slope, pvalue, rvalue, datapoints = bestSlope(sampD, dataLabels)
            
#            slopeD[samp]["slope"] = slope
#            slopeD[samp]["rvalue"] = rvalue
#            slopeD[samp]["datapoints"] = datapoints
            
            posRatio = slope/posSlope
            negRatio = slope/negSlope
            
            #Determine category
            if pvalue>0.05 or negRatio < 2:
                cat = "None"
                if negRatio > 2:
                    print("Check data for %s. Linear regression is insignifiacnt (pvalue = %f) but slope is %fx greater than negative control (%f)" % (samp, pvalue, negRatio, slope))
            
            elif posRatio > 1.1:
                cat = "VeryStrong"

            elif posRatio > 0.8:
                cat = "Strong"
            
            elif posRatio > 0.6:
                cat = "Intermediate-Strong"

            elif posRatio > 0.4:
                cat = "Intermediate"

            elif posRatio > 0.2:
                cat = "Intermediate-Weak"
            
            else:
                cat = "Weak"
            
            fout.write("%s\t%s\t%s\t%.2f\t%.6f\t%.3f\t%.3f\t%.3f\n" % (samp, cat, ",".join(datapoints), slope, pvalue, rvalue, slope/posSlope, slope/negSlope))
            

    
#     Prep negative control data for plotting on x-axis, this will be the same across all of the plots
#     if opts.negMatrix:
#         negD = io.fileDictFullRowNames(opts.negMatrix, opts.delim)
#     else:
#         negD = dataD
#     
#     if not opts.negControls:
#         negNames = list(negD.keys())
#     else:
#         negNames = opts.negControls.split(",")
#     
#     peptideNames = list(negD[negNames[0]].keys())
#     
#     x = [np.mean([float(negD[sn][pn]) for sn in negNames]) for pn in peptideNames]
# 
#     
#     Read in lists of enriched peptides
#     enrFiles = glob.glob("%s/*%s" % (opts.enrDir, opts.enrExt))
# 
#     Generate plot for each list of enriched peptides
#     for eF in enrFiles:
#         enrichedD = io.fileEmptyDict(eF, header=False)
#         sNames = eF.split("/")[-1].split(opts.enrExt)[0].split(opts.snDelim)
#         
#         Average values for y-axis
#         y = [np.mean([float(dataD[sn][pn]) for sn in sNames]) for pn in peptideNames]
#         
#         Determine color for each point
#         c = [opts.posColor if p in enrichedD else opts.negColor for p in peptideNames]
#         
#         Generate plot
#         fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w')            
#         
#         if opts.plotLog:
#             ax.scatter(np.log10(x), np.log10(y), c=c, alpha=0.5)
#         else:
#             ax.scatter(x, y, c=c, alpha=0.5)
# 
#         ax.set_ylabel(",".join(sNames), fontsize=15)
#         ax.set_xlabel(opts.xLabel, fontsize=15)
# 
#        if lim:
#            ax.set_xlim(lim[0], lim[1])
#            ax.set_ylim(lim[0], lim[1])
# 
#         fig.savefig("%s/%s.%s" % (opts.outDir, opts.snDelim.join(sNames), opts.plotType), dpi=300, bbox_inches='tight')
# 
#         fig.close()
#----------------------End of main()

def bestSlope(sampD, dataLabels):
    slope = -1000
    for i in range(0, len(dataLabels)-2):
        these = dataLabels[i:i+3]
        
        x=[]
        y=[]
        
        for each in these:
            for point in sampD[each]:
                x.append(int(each))
                y.append(point)
        
        results = linregress(x, y)
        if results.slope > slope:
            slope = results.slope
            pvalue = results.pvalue
            rvalue = results.rvalue
            datapoints = these
        
    return slope, pvalue, rvalue, datapoints

def parseData(file, delim="\t"):

    #To hold data
    dataD={x:{} for x in ["positive", "negative", "experimental"]}

    #Step through file
    with open(file, "r") as fin:
        lc=0
        for line in fin:
            lc+=1
            thisCols = line.rstrip("\n").split(delim)
            if lc==1:
                headMap = {i:x for i,x in enumerate(thisCols)}
                headMapRev = {x:i for i,x in enumerate(thisCols)}
                
                dataLabels = thisCols[3:-1]
                emptyLabel = thisCols[-1]
                
            else:
                sampName = thisCols[headMapRev["sample"]]
                sampType = thisCols[headMapRev["type"]]

                if sampType not in dataD:
                    print("Unrecognized sample type: %s" % (sampType))
                else:
                    if sampName not in dataD[sampType]:
                        dataD[sampType][sampName] = defaultdict(list)
                    
                    for i,x in enumerate(thisCols):
                        if i>2:
                            dataD[sampType][sampName][headMap[i]].append(float(x))

    return dataD, dataLabels, emptyLabel


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

