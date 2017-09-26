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

#import matplotlib as mpl
#from matplotlib import pylab as plt
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#import numpy as np
#typeface='helvetica neue lt std'
#mpl.rcParams['font.family']=typeface
#mpl.rcParams['font.size']=22

#Changed option to using R2 instead of using both reads
#In version 3.0, changed so that multiple bams can be supplied as arguments. A separate plot will be made for each bam. 
def main():

    #To parse command line
    usage = "usage: %prog [options] bam1 [bam2 ...]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
#    p.add_option('-b', '--bam', help='input bam file. Must be indexed [None]')
    p.add_option('--useR2', default=False, action='store_true', help='Use this flag if you want use R2 instead of R1 [false]')    
    p.add_option('-q', '--mapq', type='int', default=20, help='minimum mapping quality to be used [20]')
    p.add_option('-s', '--smooth', help='Use this option if you want to smooth the coverage plots. Specify window size and slide, comma separated [None]')
    
    opts, args = p.parse_args()

    for bam in args:
        opts.bam=bam
        cov_info = get_cov(opts)
        plot_cov(cov_info, opts)
###-----------------End of main()--------------------------->>>


def plot_cov(cov_info, opts):
    for ref, info in cov_info.iteritems():
        filename = '%s_%s_strandcov.pdf' % (opts.bam, ref)
        plt.plot(info['pos'], info['f_cov'], 'r-', info['pos'], info['r_cov'], 'g-',  info['pos'], [0]* len(info['pos']), 'k-')
        plt.axis([min(info['pos']), max(info['pos']), 0, max(info['f_cov']+ info['r_cov'])])
        #plt.plot(info['pos'], info['f_cov'], 'r-', info['pos'], [-1*x for x in info['r_cov']], 'g-',  info['pos'], [0]* len(info['pos']), 'k-')
        #plt.axis([min(info['pos']), max(info['pos']), min([-1*x for x in info['r_cov']]), max(info['f_cov'])])
        plt.legend(prop=fontP)
        plt.savefig(filename)
        plt.clf()
        

def get_cov(opts):
    bam = pysam.AlignmentFile(opts.bam)
    cov_info={}
    #Step through each reference sequence in the bam file
    #print dir(bam)
    for ref in bam.references:
        cov_info[ref]={'f_cov':[], 'r_cov':[], 'pos':[]}
        #Step through each position for the reference in the pileup
        for pos in bam.pileup(ref):
            f_count, r_count = get_counts(pos.pileups, opts)
            cov_info[ref]['pos'].append(pos.reference_pos+1)
            cov_info[ref]['f_cov'].append(f_count)
            cov_info[ref]['r_cov'].append(r_count)
    if opts.smooth: cov_info = smooth(cov_info, [int(x) for x in opts.smooth.split(',')])
    return cov_info
###------------->>>

def average(list):
    return sum(list)/len(list)

def smooth(cov_info, smooth_params):
    win,step = smooth_params
    new_cov_info={}
    for ref, info in cov_info.iteritems():
        new_pos=[]
        new_f=[]
        new_r=[]
        start=0
        while start<=len(info['pos'])-win:
            new_pos.append(average(info['pos'][start:start+win]))
            new_f.append(average(info['f_cov'][start:start+win]))
            new_r.append(average(info['r_cov'][start:start+win]))
            start+=step
        new_cov_info[ref] = {'f_cov':new_f, 'r_cov':new_r, 'pos':new_pos}
    return new_cov_info

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
