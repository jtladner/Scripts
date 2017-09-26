#!/usr/bin/env python

# By Jason Ladner

from __future__ import division
import sys, optparse, os, pysam
#For plotting
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as fp
fontP = fp()
#Numpy for stdev
import numpy as np
import scipy
import scipy.stats

#import matplotlib as mpl
#from matplotlib import pylab as plt
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#import numpy as np
#typeface='helvetica neue lt std'
#mpl.rcParams['font.family']=typeface
#mpl.rcParams['font.size']=22

#Changed option to using R2 instead of using both reads
#In v3, I added an option to supply a file with a bunch of bam file locations included. 
#In v4, changed the way I'm normalizing. Now normalizing by just taking the average of all the average coverages calculated for the ORFs

def main():

    #To parse command line
    usage = "usage: %prog [options] bam1 [bam2 ...]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
#    p.add_option('-b', '--bam', help='input bam file. Must be indexed [None]')
    p.add_option('--useR2', default=False, action='store_true', help='Use this flag if you want use R2 instead of R1 [false]')
    p.add_option('--sem', default=False, action='store_true', help='Display confidence intervals as standard error, as opposed to standard deviation [false]')
    p.add_option('-q', '--mapq', type='int', default=20, help='minimum mapping quality to be used [20]')
    p.add_option('-a', '--annots', help='File with names and coordinates for regions to include in plot [None]')
    p.add_option('-b', '--bams', help='File with a bunch of bam file locations. Can be used instead of providing bams as arguments [None]')
    p.add_option('-o', '--out', help='Base name for output files, both a pdf plot and the values used in the plot [None]')
    
    
    opts, args = p.parse_args()

    #Step through each bam file and pull out the coverage at each site
    all_covs=[]
    
    if opts.bams:
        for line in open(opts.bams, 'r'):
            #all_covs is a list containing dictionaries: keys(1)=ref seq names, values(1) = dictionary: keys(2)=ref positions, values(2)=[forward_read_count, reverse_read_count]
            all_covs.append(get_cov(line.strip(), opts))
    
    if args:
        for each in args:
            #all_covs is a list containing dictionaries: keys(1)=ref seq names, values(1) = dictionary: keys(2)=ref positions, values(2)=[forward_read_count, reverse_read_count]
            all_covs.append(get_cov(each, opts))
    
    #For output of avg covs per bam
    #Headers assume that reference sequence is positive-sense
#    print "Gene\t%s\t%s" % ("\t".join(["Pos-%s"  % ('_'.join(x.split('_')[:3])) for x in args]), "\t".join(["Neg-%s"  % ('_'.join(x.split('_')[:3])) for x in args]))
    
    #Read in annotation info
    annots=read_annots(opts.annots)
    
    #Create a plot for each reference contig
    for ref in all_covs[0].keys():
        if ref in annots:
            x,y_pos,yerr_pos,y_neg,yerr_neg = info_for_plot([x[ref] for x in all_covs], annots[ref], ref, opts)
            plot_cov('%s_%s.pdf' % (opts.out, ref), x,y_pos,yerr_pos,y_neg,yerr_neg, [x[0] for x in annots])
###-----------------End of main()--------------------------->>>

def info_for_plot(covs, annots, ref, opts):
    x_coord=[]
    y_pos=[]
    yerr_pos=[]
    y_neg=[]
    yerr_neg=[]
    
    #For output txt file
    fout=open('%s_%s.txt' % (opts.out, ref), 'w')
    #Write header for file that will contain the data used to make the plot
    if opts.sem: fout.write("Region\tPos\tRnorm\tRsem\tFnorm\tFsem\n")
    else: fout.write("Region\tPos\tRnorm\tRstd\tFnorm\tFstd\n")
    
    #Lists that will hold the average coverages for the different genes
    pos_norm_bybam=[]
    neg_norm_bybam=[]

    #Step through each bam
    for bam in covs:
        #Step through each gene in the annotation file
        pos_avg_covs, neg_avg_covs = calc_covs(bam, annots)
        #Normalize by average average coverage across all ORFs
        pos_norm_bybam.append([x/np.mean(pos_avg_covs) for x in pos_avg_covs])
        neg_norm_bybam.append([x/np.mean(neg_avg_covs) for x in neg_avg_covs])

#To get x coordinates for plot
    for gene, start, stop in annots:
        x_coord.append((start+stop)/2)

    #Step through each gene
    for index in range(len(annots)):
        y_pos.append(np.mean([x[index] for x in pos_norm_bybam]))
        y_neg.append(np.mean([x[index] for x in neg_norm_bybam]))
        if opts.sem:
            yerr_pos.append(scipy.stats.sem([x[index] for x in pos_norm_bybam]))
            yerr_neg.append(scipy.stats.sem([x[index] for x in neg_norm_bybam]))
        else:
            yerr_pos.append(np.std([x[index] for x in pos_norm_bybam]))
            yerr_neg.append(np.std([x[index] for x in neg_norm_bybam]))
        #Write info to an output txt file
        fout.write("%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\n" % (annots[index][0], x[index], y_pos[index], yerr_pos[index], y_neg[index], yerr_neg[index]))

#    print "%s\t%s\t%s" % (gene, "\t".join([str(x) for x in pos_covs]), "\t".join([str(x) for x in neg_covs]))

    return x_coord, y_pos, yerr_pos, y_neg, yerr_neg

#Returns a list of avg cov values
def calc_covs(each, annots):
    pos_covs=[]
    neg_covs=[]
    #Step through each annotated gene
    for gene, start, stop in annots:
        #***Assuming that the reference was positive sense
        #Calculates and appends the average positive strand coverage across this one ORF
        pos_covs.append(avg_cov(each, int(start), int(stop), 1))
        #Calculates and appends the average neagtive strand coverage across this one ORF
        neg_covs.append(avg_cov(each, int(start), int(stop), 0))
    return pos_covs, neg_covs

def avg_cov(each, start, stop, index):
    all_b_covs=[]
    for b in range(start, stop+1):
        if b in each: all_b_covs.append(each[b][index])
        #Have this to account for bases with no coverage not being included in the pileup
        else: all_b_covs.append(0)
    return np.mean(all_b_covs)

def plot_cov(filename, x,y_pos,yerr_pos,y_neg,yerr_neg, names):
    
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(x, y_pos, 'r-', x, y_neg, 'g-')
    ax1.errorbar(x, y_pos, yerr=yerr_pos, fmt='o', color='r')
    ax1.errorbar(x, y_neg, yerr=yerr_neg, fmt='o', color='g')
    ax1.set_xlabel('Genes')
    ax1.set_ylabel('Normalized coverage')
    ax1.set_xticks(x)
    ax1.set_xticklabels(names)
    fig.savefig(filename)
    fig.clf()
    
    #plt.plot(x, y_pos, 'r-', x, y_neg, 'g-')
    #plt.errorbar(x, y_pos, yerr=yerr_pos)
    #plt.errorbar(x, y_neg, yerr=yerr_neg)

    #plt.axis([0, max(x), 0, max(info['f_cov']+ info['r_cov'])])
    #plt.legend(prop=fontP)
    #plt.savefig(filename)
    #plt.clf()
        

#Just raw counts being output, no normalization
def get_cov(afile, opts):
    bam = pysam.AlignmentFile(afile)
    cov_info={}
    #Step through each reference sequence in the bam file
    #print dir(bam)
    for ref in bam.references:
        cov_info[ref]={}
        #Step through each position for the reference in the pileup
        for pos in bam.pileup(ref):
            # # of forward and reverse strand reads covering this position
            f_count, r_count = get_counts(pos.pileups, opts)
            #I may be one off here, NEED TO CHECK
            cov_info[ref][pos.reference_pos] = [f_count, r_count]
    return cov_info

#Adjusted in v4 so that all of the annotations are linked to specific references
def read_annots(afile):
    info={}
    for line in open(afile, 'r'):
        cols=line.strip().split('\t')
        if cols[0] not in info: info[cols[0]] = []
        info[cols[0]].append([cols[1], int(cols[2]), int(cols[3])])
    return info
###------------->>>

def average(list):
    return sum(list)/len(list)

def get_counts(pileups, opts):
    f_count=0
    r_count=0
    for read in pileups:
        if read.alignment.mapping_quality >= opts.mapq and not read.alignment.is_unmapped:
            if (not opts.useR2 and read.alignment.is_read1) or (opts.useR2 and not read.alignment.is_read1):
                if read.alignment.is_reverse: r_count+=1
                else: f_count+=1
    return f_count, r_count

if __name__ == "__main__":
    main()
