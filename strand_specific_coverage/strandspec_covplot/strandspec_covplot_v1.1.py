#!/usr/bin/env python

# By Jason Ladner

from __future__ import division
import sys, optparse, os, pysam
import numpy as np
#For plotting
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as fp
fontP = fp()

#In version 1.1, added option to utilize unpaired reads (not flagged as read1 or read2)

#!#!#! From pysam website: "Coordinates in pysam are always 0-based (following the python convention). SAM text files use 1-based coordinates."

#For future, add an optional average coverage threshold for inclusion in combo plot

def main():

    #To parse command line
    usage = "usage: %prog [options] bam1 [bam2 ...]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
#    p.add_option('-b', '--bam', help='input bam file. Must be indexed [None]')
    p.add_option('-o', '--out', help='Base name for output files [None]')
    p.add_option('-q', '--mapq', type='int', default=20, help='minimum mapping quality to be used [20]')
    p.add_option('-s', '--smooth', help='Use this option if you want to smooth the coverage plots. Specify window size and slide, comma separated [None]')
    p.add_option('--useR2', default=False, action='store_true', help='Use this flag if you want use R2 instead of R1 [false]')    
    p.add_option('--useUnpaired', default=False, action='store_true', help='Use this flag if you want use reads not specified as R1 or R2 [false]')    
    p.add_option('--noIndivPlots', default=False, action='store_true',  help='Use this flag if you do not want to create plots for individual bams [false]')
#    p.add_option('--plotRev', default=False, action='store_true', help='Use this flag if you want to plot coverage for reads mapped to the reverse strand, instead of the forward strand. For RNA Access, reverse reads represent the strand of the reference genome. [false]')    
#    p.add_option('--dontNorm', default=False, action='store_true', help='Use this flag if you do not want to normalize coverage at each position by the average coverage across the genome/segment [false]')    
    
    opts, args = p.parse_args()

    fo_raw = open("%s_rawcov.txt" % (opts.out), 'w')
    fo_raw.write("File\tReference\tPosition(1-based)\tForCov\tRevCov\n")
    fo_norm = open("%s_normcov.txt" % (opts.out), 'w')
    fo_norm.write("File\tReference\tPosition(1-based)\tNormForCov\tNormRevCov\n")


    all_norm = {}
    for bam in args:
        print bam
        opts.bam=bam
        cov_info, norm_info = get_cov(opts)
        all_norm[bam] = norm_info
        for ref in cov_info:
            print ref
            for i, p in enumerate(cov_info[ref]["pos"]):
                fo_raw.write("%s\t%s\t%d\t%d\t%d\n" % (bam, ref, p, cov_info[ref]["f_cov"][i], cov_info[ref]["r_cov"][i]))
                fo_norm.write("%s\t%s\t%d\t%.4f\t%.4f\n" % (bam, ref, p, norm_info[ref]["f_cov"][i], norm_info[ref]["r_cov"][i]))
        if not opts.noIndivPlots:
            plot_cov(cov_info, 'rawcov', opts)
            print "Raw plot done for %s, %s" % (bam, ref)
            plot_cov(norm_info, 'normcov', opts)
            print "Norm plot done for %s, %s" % (bam, ref)
    
    #To troubleshoot bug
    #!#!#!#! Increasing the size of the sliding window solved the isue I was having, but still don't understand what the problem was
#    for ref in all_norm[bam].keys():
#        print ref, len(all_norm[bam][ref]['pos']), len(all_norm[bam][ref]['f_cov']), len(all_norm[bam][ref]['r_cov'])
    
    #For the combined plot, calculate the average and stdev for each window/position
    comb_norm = {}
    for ref in all_norm[bam].keys():
        comb_norm[ref]={'f_cov':np.array([]),'f_std':np.array([]), 'r_cov':np.array([]), 'r_std':np.array([]),'pos':np.array([])}
        for i, p in enumerate(all_norm[bam][ref]['pos']):
            fthis=[]
            rthis=[]
            for k, v in all_norm.iteritems():
                fthis.append(v[ref]['f_cov'][i])
                rthis.append(v[ref]['r_cov'][i])
            comb_norm[ref]['f_cov'] = np.append(np.mean(fthis), comb_norm[ref]['f_cov'])
            comb_norm[ref]['f_std'] = np.append(np.std(fthis), comb_norm[ref]['f_std'])
            comb_norm[ref]['r_cov'] = np.append(np.mean(rthis), comb_norm[ref]['r_cov'])
            comb_norm[ref]['r_std'] = np.append(np.std(rthis), comb_norm[ref]['r_std'])
            comb_norm[ref]['pos'] = np.append(p, comb_norm[ref]['pos'])
    print "Done combining"

    plot_cov_std(comb_norm, 'combo', opts)
    
    fo_raw.close()
    fo_norm.close()
###-----------------End of main()--------------------------->>>

#Making separate plots for the reverse and forward strands
def plot_cov_std(cov_info, outstr, opts):
    for ref, info in cov_info.iteritems():
        #forward strand
        filename = '%s_%s_For_%s.pdf' % (opts.out, ref, outstr)
#        print info['pos']
#        print info['f_cov']
        
        plt.plot(info['pos'], info['f_cov'], '-', color="#da7c30")
        plt.fill_between(info['pos'], info['f_cov']+info['f_std'], info['f_cov']-info['f_std'], facecolor="#da7c30", alpha=0.5)
        plt.axis([min(info['pos']), max(info['pos']), 0, max(info['f_cov'])+ max(info['f_std'])])
        #plt.plot(info['pos'], info['f_cov'], 'r-', info['pos'], [-1*x for x in info['r_cov']], 'g-',  info['pos'], [0]* len(info['pos']), 'k-')
        #plt.axis([min(info['pos']), max(info['pos']), min([-1*x for x in info['r_cov']]), max(info['f_cov'])])
        plt.legend(prop=fontP)
        plt.savefig(filename)
        plt.clf()

        #reverse strand
        filename = '%s_%s_Rev_%s.pdf' % (opts.out, ref, outstr)
        plt.plot(info['pos'], info['r_cov'], '-', color="#396ab1")
        plt.fill_between(info['pos'], info['r_cov']+info['r_std'], info['r_cov']-info['r_std'], facecolor="#396ab1", alpha=0.5)
        plt.axis([min(info['pos']), max(info['pos']), 0, max(info['r_cov'])+ max(info['r_std'])])
        #plt.plot(info['pos'], info['f_cov'], 'r-', info['pos'], [-1*x for x in info['r_cov']], 'g-',  info['pos'], [0]* len(info['pos']), 'k-')
        #plt.axis([min(info['pos']), max(info['pos']), min([-1*x for x in info['r_cov']]), max(info['f_cov'])])
        plt.legend(prop=fontP)
        plt.savefig(filename)
        plt.clf()
        

def plot_cov(cov_info, outstr, opts):
    for ref, info in cov_info.iteritems():
        filename = '%s_%s_%s.pdf' % (opts.bam, ref, outstr)
        plt.plot(info['pos'], info['f_cov'], 'r-', info['pos'], info['r_cov'], 'g-',  info['pos'], [0]* len(info['pos']), 'k-')
        plt.axis([min(info['pos']), max(info['pos']), 0, max(info['f_cov']+ info['r_cov'])])
        #plt.plot(info['pos'], info['f_cov'], 'r-', info['pos'], [-1*x for x in info['r_cov']], 'g-',  info['pos'], [0]* len(info['pos']), 'k-')
        #plt.axis([min(info['pos']), max(info['pos']), min([-1*x for x in info['r_cov']]), max(info['f_cov'])])
        plt.legend(prop=fontP)
        plt.savefig(filename)
        plt.clf()
        

def get_cov(opts):
    bam = pysam.AlignmentFile(opts.bam)
    
    #Make dict with read length info
    #!#!#! This could potentially cause a problem is the sam/bam file does not include a header with info on reference lengths
    #!#!#! I should probably update this section to include a warning and work around if this is the case
    reflen_dict={}
    for index, r in enumerate(bam.references):
        reflen_dict[r]=bam.lengths[index]

    #To hold coverage information
    cov_info={}
    norm_info={}

    #Step through each reference sequence in the bam file
    #print dir(bam)
    for ref in bam.references:
        cov_info[ref]={'f_cov':[], 'r_cov':[], 'pos':[]}
        #Step through each position for the reference in the pileup
        poscount=0
        for pos in bam.pileup(ref):
            poscount+=1
            #To account for bases at the beginning of the ref with 0 coverage
            if poscount == 1 and poscount != pos.reference_pos+1:
                for i in range(1, pos.reference_pos+1):
                    cov_info[ref]['pos'].append(i)
                    cov_info[ref]['f_cov'].append(0)
                    cov_info[ref]['r_cov'].append(0)

            f_count, r_count = get_counts(pos.pileups, opts)
            cov_info[ref]['pos'].append(pos.reference_pos+1)
            cov_info[ref]['f_cov'].append(f_count)
            cov_info[ref]['r_cov'].append(r_count)
    
        #To account for bases at the beginning of the ref with 0 coverage
        if poscount+1 != reflen_dict[ref]:
            if poscount+1 > reflen_dict[ref]: print "!!!Something is not right!!! Last pileup position (%d) should not be larger than the length of the reference (%d)" % (poscount+1, reflen_dict[ref])
            else: 
                for i in range(poscount+2, reflen_dict[ref]+1):
                    cov_info[ref]['pos'].append(i)
                    cov_info[ref]['f_cov'].append(0)
                    cov_info[ref]['r_cov'].append(0)
        
        #Create version where the coverages are normalized by the average coverage across the reference sequence
        norm_info[ref]={'f_cov':[], 'r_cov':[], 'pos':[]}
        f_avg = sum(cov_info[ref]['f_cov'])/len(cov_info[ref]['f_cov'])
        r_avg = sum(cov_info[ref]['r_cov'])/len(cov_info[ref]['r_cov'])
        for i, p in enumerate(cov_info[ref]['pos']):
            norm_info[ref]['pos'].append(p)
            norm_info[ref]['f_cov'].append(cov_info[ref]['f_cov'][i]/f_avg)
            norm_info[ref]['r_cov'].append(cov_info[ref]['r_cov'][i]/r_avg)

    if opts.smooth: 
        cov_info = smooth(cov_info, [int(x) for x in opts.smooth.split(',')])
        norm_info = smooth(norm_info, [int(x) for x in opts.smooth.split(',')])
    return cov_info, norm_info
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
            if (not opts.useR2 and read.alignment.is_read1) or (opts.useR2 and not read.alignment.is_read1) or (opts.useUnpaired and not read.alignment.is_read1 and not read.alignment.is_read2):
                if read.alignment.is_reverse: r_count+=1
                else: f_count+=1
    return f_count, r_count

if __name__ == "__main__":
    main()
