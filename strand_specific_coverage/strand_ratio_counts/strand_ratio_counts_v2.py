#!/usr/bin/env python

# By Jason Ladner

from __future__ import division
import sys, optparse, os, pysam

#For plotting
#import matplotlib
#matplotlib.use('PDF')
#import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties as fp
#fontP = fp()
#Numpy for stdev
#import numpy as np

#In version 2, added ability to limit calculations to a subset region of the genome

#Provide with a text file that has a head and three tab-delimited columns
    #Only the third column is used and should be bam file locations.
    #This script will calculate the ratio of reads mapped to the various strands and will add this info to the table in the ouput

def main():

    #To parse command line
    usage = "usage: %prog [options] name1,bamnamefile1 [name2,bamnamefile2 ...]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
    p.add_option('-i', '--inp', help='Name for input file [None]')
    p.add_option('-o', '--out', help='Name for output file [None]')
#    p.add_option('-m', '--mapQual', default=20, type='int', help='Minimum mapping quality for a read to be used [20]')
    #These are optional, only if you want to only look at reads that overlap a certain region
    p.add_option('-b', '--beg', type='int', help='1st base in range')
    p.add_option('-e', '--end', type='int', help='last base in range')
    p.add_option('-v', '--overlap', type='int', default=10, help='required overlap for reads to be inlcuded [10]')

    opts, args = p.parse_args()

    #open file for output
    fout = open(opts.out, 'w')
    
    linecount=0
    for line in open(opts.inp, 'r'):
        linecount+=1
        if linecount==1: fout.write("%s\tPropForw\tPropRev\tStrandRatio\tNumFor\tNumRev\n" % (line.strip()))
        else:
            cols=line.strip().split('\t')
            proprev, propforw, ratio, forw, rev = parse_strand_indiv(cols[2], opts)
            fout.write("%s\t%.4f\t%.4f\t%.4f\t%d\t%d\n" % ("\t".join(cols), propforw, proprev, ratio, forw, rev))
    fout.close()

    #Make figure
    #fig = plt.figure()
    #ax1 = fig.add_subplot(1,1,1)
    #ax1.boxplot(values)
    #ax1.set_xlabel('Sample Types')
    #ax1.set_ylabel('# Positive sense reads/#Negative sense reads')
    #ax1.set_xticks(range(1,len(names)+1))
    #ax1.set_xticklabels(names)
    #fig.savefig(opts.out)
    #fig.clf()


###-----------------End of main()--------------------------->>>

#def parse_regions(coord_file):
#    counts_dict={}
#    for line in open(coord_file, 'r'):
#        cols=line.strip().split('\t')
#        counts_dict[tuple(cols)]=0
#    return counts_dict

def parse_strand_indiv(each, opts):
    ratios=[]
    sambam = pysam.AlignmentFile(each)
    num_forward=0
    num_reverse=0
    #Step through each read in the sam file
    for read in sambam:
        #Check to make sure the read is mapped and forward strand, which should be the first mate in the file
        #Checks whether the read overlaps with the 
        if read.is_proper_pair and read.is_read1: 
            if not opts.beg:
                if read.is_reverse: num_reverse+=1
                elif read.mate_is_reverse: num_forward+=1
                else: print "PROBLEM: neither read in pair is mapped to reverese strand"
            elif read.get_overlap(opts.beg, opts.end)>=opts.overlap:
                if read.is_reverse: num_reverse+=1
                elif read.mate_is_reverse: num_forward+=1
                else: print "PROBLEM: neither read in pair is mapped to reverese strand"
    if num_reverse==0 and num_forward==0: 
        print "NO reads covering specified region for %s!!!!!!!" % (each)
        return -99, -99, -99, num_forward, num_reverse
    elif num_forward==0: 
        print "NO forward strand reads for %s!!!!!!!!!!!" % (each)
        return num_reverse/(num_reverse+num_forward), num_forward/(num_reverse+num_forward), -99, num_forward, num_reverse
    elif num_reverse==0: 
        print "NO reverse strand reads for %s!!!!!!!!!!!" % (each)
    return num_reverse/(num_reverse+num_forward), num_forward/(num_reverse+num_forward), num_reverse/num_forward, num_forward, num_reverse
###------------->>>

if __name__ == "__main__":
    main()
