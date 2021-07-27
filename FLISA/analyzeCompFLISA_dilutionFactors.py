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

errorMessages = {
    "posNS": "The positive control sample resulted in a non-significant slope!",
    "negS": "The negative control sample resulted in a significant slope!",
    "posHigh": "The positive control sample resulted in higher than expected slope!",
    "negLow": "The negative control sample resulted in lower than expected slope!",
    "poorCorr": "The correlation coefficient for this sample is low!",
    "moreDil":  "Consider including more dilutions for this sample!",
    "moreConc":  "Used three highest concentrations, consider changing dilutions!"
}

# This script is used for processing data from competitive FLISA assays with multiple dilutions per sample

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

    p.add_option('-d', '--data',  help='File containing data from an assay. Should include both a positive and negative control sample. Columns: sample_name, rep#, sample_type, datapoints, associated_blank [None, REQ]')
    p.add_option('--delim', default="\t", help="Delimiter used in the data file. [\\t]")
    p.add_option('--maxPosSlope', type="float", default=-1000, help="If the slope for the positive control sample is greater than this, a warning will be issued. [-1000]")
    p.add_option('--minNegSlope', type="float", default=-100, help="If the slope for the negative control sample is lower than this, a warning will be issued. [-100]")
    p.add_option('--minRvalue', type="float", default=0.97, help="If the correlation coefficient for a sample is less than this, a warning will be issued. [0.97]")
    p.add_option('--maxSlopeRatio', type="float", default=1.2, help="If the ratio of the (slope from empty control and most dilute sample)/(best dilution slope) is greater than this, a warning will be issued. [1.2]")

    #Dilution factors
    p.add_option('--expDilFactor', type="int", default=3, help="Dilution factor for experimental samples. [3]")
    p.add_option('--negDilFactor', type="int", default=3, help="Dilution factor for negative sample. [3]")
    p.add_option('--posDilFactor', type="int", default=2, help="Dilution factor for positive sample. [2]")

    opts, args = p.parse_args()
        
    #Define colors for experimental samples
    colL = ["#cc79a7", "#0072b2", "#f0e442"]
    
    # Print command used 
    print("Command run: '%s'" % ("  ".join(sys.argv)))
    
    #Read in data file 
    dataD, dataLabels, emptyLabel = parseData(opts.data, delim=opts.delim)
    
    # Check to make sure there is only one positive control, one negative control and at least one experimental sample
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
        #Create list to hold warnings
        warnings = []
                    
        #Calculate slope, etc. for positive control sample
        posName = list(dataD["positive"].keys())[0]

        posSlope, posPvalue, posRvalue, posYinter, posDatapoints = bestSlope(dataD["positive"][posName], dataLabels, opts.posDilFactor)
        
        # Some quality checks
        if posPvalue > 0.05:
            warnings.append(errorMessages["posNS"])
        if posSlope > opts.maxPosSlope:
            warnings.append(errorMessages["posHigh"])
        if abs(posRvalue) < opts.minRvalue:
            warnings.append("%s: Positive" % (errorMessages["poorCorr"]))
        if posDatapoints == dataLabels[:3]:
            warnings.append("%s: Positive" % (errorMessages["moreConc"]))
        
        
        #Compare difference between empty control and most dilute sample, to see if additional dilutions may be needed
        empPosSlope = simpSlope({float(dataLabels[-1]):dataD["positive"][posName][dataLabels[-1]], float(dataLabels[-1])+1:dataD["positive"][posName]["N"]})
            
        # Check if we may need more dilutions for the positive control
        if empPosSlope/posSlope > opts.maxSlopeRatio:
            warnings.append("%s: Positive" % (errorMessages["moreDil"]))


        #Calculate slope, etc. for negative control sample
        negName = list(dataD["negative"].keys())[0]
        negSlope, negPvalue, negRvalue, negYinter, negDatapoints = bestSlope(dataD["negative"][negName], dataLabels, opts.negDilFactor)

        # Some quality checks
        if negPvalue <= 0.05:
            warnings.append(errorMessages["negS"])
        if negSlope < opts.minNegSlope:
            warnings.append(errorMessages["negLow"])

    #initiate figure
    fig,ax = plt.subplots(1,1,figsize=(5,5),facecolor='w')            


    with open("%s_results.tsv" % (".".join(opts.data.split(".")[:-1])), "w") as fout:
        

        # Write header and control info to output file
        fout.write("Sample\tCategory\tDatapointsUsed\tSlope\ty-intercept\tPvalue\tCorrelationCoefficient\tPosSlopeRatio\tNegSlopeRatio\n")
        fout.write("%s\t%s\t%s\t%.2f\t%.2f\t%.6f\t%.3f\t%.3f\t%.3f\n" % (posName, "Positive",  ",".join(posDatapoints), posSlope, posYinter, posPvalue, posRvalue, posSlope/posSlope, posSlope/negSlope))
        fout.write("%s\t%s\t%s\t%.2f\t%.2f\t%.6f\t%.3f\t%.3f\t%.3f\n" % (negName, "Negative",  ",".join(negDatapoints), negSlope, negYinter, negPvalue, negRvalue, negSlope/posSlope, negSlope/negSlope))

        # Plot control data and set axis labels        
        xIntsPos = np.array([np.log10(1/(opts.posDilFactor**(float(x)-1))) for x in dataLabels])
        xIntsNeg = np.array([np.log10(1/(opts.negDilFactor**(float(x)-1))) for x in dataLabels])
        xIntsExp = np.array([np.log10(1/(opts.expDilFactor**(float(x)-1))) for x in dataLabels])

        plotPoints(ax, {k:dataD["positive"][posName][k] for k in dataLabels}, "#009e73", posDatapoints, "positive", opts.posDilFactor)
        ax.plot(xIntsPos, posYinter + posSlope*xIntsPos, c="#009e73", alpha=0.6, linestyle=":")
        plotPoints(ax, {k:dataD["negative"][negName][k] for k in dataLabels}, "#d55e00", negDatapoints, "negative", opts.negDilFactor)
        ax.plot(xIntsNeg, negYinter + negSlope*xIntsNeg, c="#d55e00", alpha=0.6, linestyle=":")

        
        
        ax.set_ylabel("Signal", fontsize=15)
        ax.set_xlabel("Dilution", fontsize=15)

        slopeD={}
        
        for samp, sampD in dataD["experimental"].items():
            slopeD[samp] = {}
            slope, pvalue, rvalue, yinter, datapoints = bestSlope(sampD, dataLabels, opts.expDilFactor)
            
            #Compare difference between empty control and most dilute sample, to see if additional dilutions may be needed
            empSlope = simpSlope({float(dataLabels[-1]):sampD[dataLabels[-1]], float(dataLabels[-1])+1:sampD["N"]})
                
            # Some quality checks
            if abs(rvalue) < opts.minRvalue:
                warnings.append("%s: %s" % (errorMessages["poorCorr"], samp))
            if empSlope/slope > opts.maxSlopeRatio:
                warnings.append("%s: %s" % (errorMessages["moreDil"], samp))
            if datapoints == dataLabels[:3]:
                warnings.append("%s: %s" % (errorMessages["moreConc"], samp))

#            slopeD[samp]["slope"] = slope
#            slopeD[samp]["rvalue"] = rvalue
#            slopeD[samp]["datapoints"] = datapoints
            
            posRatio = slope/posSlope
            negRatio = slope/negSlope
            
            #Determine category
            if pvalue>0.05 or (0 < negRatio < 2):
                cat = "None"
                if abs(negRatio) > 2:
                    print("\nCheck data for %s. Linear regression is insignifiacnt (pvalue = %f) but slope is %fx greater than negative control (%f)" % (samp, pvalue, negRatio, slope))
            
            elif posRatio > 0.8:
                cat = "VeryStrong"

            elif posRatio > 0.6:
                cat = "Strong"
            
#            elif posRatio > 0.6:
#                cat = "Intermediate-Strong"

            elif posRatio > 0.4:
                cat = "Intermediate"

            elif posRatio > 0.2:
                cat = "Weak"
            
            else:
                cat = "VeryWeak"
            
            fout.write("%s\t%s\t%s\t%.2f\t%.2f\t%.6f\t%.3f\t%.3f\t%.3f\n" % (samp, cat, ",".join(datapoints), slope, yinter, pvalue, rvalue, slope/posSlope, slope/negSlope))
            
            # Plot points
            plotPoints(ax, {k:sampD[k] for k in dataLabels}, colL[len(slopeD)-1], datapoints, samp, opts.expDilFactor)

            # Linestyle is defined based on strength of interaction
            if cat in ["VeryStrong", "Strong"]:
                ls="-"
            elif cat in ["Intermediate-Strong", "Intermediate"]:
                ls="-."
            elif cat in ["Intermediate-Weak", "Weak"]:
                ls="--"
            else:
                ls=":"
            
            # Plot line
            ax.plot(xIntsExp, yinter + slope*xIntsExp, c=colL[len(slopeD)-1], alpha=0.6, linestyle=ls)

            
    #Save figure
    ax.legend()
    fig.savefig("%s_results.pdf" % (".".join(opts.data.split(".")[:-1])),dpi=300,bbox_inches='tight')
    
#        if lim:
#            ax.set_xlim(lim[0], lim[1])
#            ax.set_ylim(lim[0], lim[1])

    # Report encountered warnings
    if warnings:
        print("\n***WARNING***\n%s" % ("\n".join(warnings)))

#----------------------End of main()

#Keys need to be floats or ints
def simpSlope(sampD):
    x=[]
    y=[]
    
    for dil, dL in sampD.items():
        for point in dL:
            x.append(dil)
            y.append(point)
    
    results = linregress(x, y)
        
    return results.slope


def bestSlope(sampD, dataLabels, dilFactor):
    slope = 10000
    for i in range(0, len(dataLabels)-2):
        these = dataLabels[i:i+3]
        
        x=[]
        y=[]
        
        for each in these:
            for point in sampD[each]:
                x.append(np.log10(1/(dilFactor**(float(each)-1))))
                y.append(point)
        
        results = linregress(x, y)
        if results.slope < slope:
            slope = results.slope
            pvalue = results.pvalue
            rvalue = results.rvalue
            yinter = results.intercept
            datapoints = these
        
    return slope, pvalue, rvalue, yinter, datapoints

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

def plotPoints(axObj, aDataD, color, pointsUsed, lab, dilFactor):
    x=[]
    y=[]
    
    
    for k, v in aDataD.items():
        for each in v:
#            x.append(float(k))
            x.append(np.log10(1/(dilFactor**(float(k)-1))))
            y.append(float(each))
    
    axObj.scatter(x, y, c=color, alpha=0.6, zorder=4, label=lab)

    xUsed=[]
    yUsed=[]

    for k in pointsUsed:
        for each in aDataD[k]:
#            xUsed.append(float(k))
            xUsed.append(np.log10(1/(dilFactor**(int(k)-1))))
            yUsed.append(float(each))
    
    axObj.scatter(xUsed, yUsed, facecolors='none', edgecolors="black", alpha=0.6, zorder=5)
    


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

