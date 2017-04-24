#!/usr/bin/env python

# By Jason Ladner

#This script looks for single reads that are chimeric, with different portions mapping to different regions

#In v2, I tweaked the way I was making the plot in matplot lib to try to make it easy to create a legend
#In v3, Fixed bug related to reads with doubly aligned bases. Only counted the left most alignments
#   **** Known issue related with mutations (I think) near breakpoints, that lead to unaligned, internal bases. This could be accounted for, but it's unclear whether bases should be added to left or right.
#In v3.1, Slightly changed the way the sizes for the legend circles is calculated to avoid non-integers

from __future__ import division
import sys, optparse, os, pysam, math
import numpy as np
from scipy.stats import gaussian_kde

#For plotting
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as fp
fontP = fp()


def main():

    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
    
    p.add_option('-b', '--bam', help='input sam or bam file. [None, REQ]')
    p.add_option('-o', '--out', default="chimerareads", help='String to add to bam names for output. [chimerareads]')
    p.add_option('-m', '--minPerc', type='float', default=0.99, help='Minumum percent of query read that must be aligned to the reference, when considering all alignemnts. [0.99]')
    
    opts, args = p.parse_args()

    SA_info = {}
    all_aligned = {}
    
    bam = pysam.AlignmentFile(opts.bam)
    #Step through all of the reads
    for read in bam.fetch():
        #To keep track of all mapped reads
        if not read.is_unmapped: 
            if read.reference_name not in all_aligned: all_aligned[read.reference_name]={}
            all_aligned[read.reference_name][read.query_name]=''
        
        #Check to make sure that 1) the read is mapped and 2) the read has a supplementary alignment
        if not read.is_unmapped and read.has_tag("SA"):
            #Add read name to dict if not already present
            if read.query_name not in SA_info: SA_info[read.query_name]=[]
            
            #Checks for hard clipped bases, and corrects pairs and length for this, if present
            if 5 in [x[0] for x in read.cigartuples]:
                corr_pairs, corr_length = corr_hardclipped(read)
                #Add info about this alignment
                SA_info[read.query_name].append([read.reference_name, corr_length, corr_pairs, read.cigarstring])
                
            else:
                #Add info about this alignment
                SA_info[read.query_name].append([read.reference_name, read.infer_query_length(), read.get_aligned_pairs(matches_only=True), read.cigarstring])
    
    goodcount=0
    fout_del=open("%s_dels.txt" % opts.out, 'w')
    fout_dup=open("%s_dups.txt" % opts.out, 'w')
    fout_del.write("ReadName\tType\tRefName\tDelLength\tDelLeft\tDelRight\tReadLength\tLeftLen\tRightLen\tOverlap\n")
    fout_dup.write("ReadName\tType\tRefName\tDupLength\tDupLeft\tDupRight\tReadLength\tLeftLen\tRightLen\n")
    for read, info in SA_info.iteritems():
        qsites_aligned, qwithdups = get_aligned([x[2] for x in info], 0)
        #rsites_aligned, rwithdups = get_aligned([x[2] for x in info], 1)
        
        #Check to make sure that there are only two alignments and that they are both to the same reference
        if len([x[0] for x in info]) == 2 and len(set([x[0] for x in info])) == 1:
            #Check to make sure that the minPerc of the read is aligned, when considering both alignments
            if len(qsites_aligned)/max([x[1] for x in info]) >= opts.minPerc:
                goodcount+=1
                missing = missing_bases(info[0][2], info[1][2])
                if missing: fout_del.write("%s\tDeletion\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (read, info[0][0], len(missing), min(missing)+1, max(missing)+1, max([x[1] for x in info]), len(info[0][2]), len(info[1][2]), len(qwithdups) - len(qsites_aligned)))
                else:
                    ####!!!! Not currently using this portion
                    ####!!!! Do I need to adjust?!?!?!
                    common = common_bases([x[1] for x in info[0][2]], [x[1] for x in info[1][2]])
                    fout_dup.write("%s\tDuplication\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (read, info[0][0], len(common), min(common)+1, max(common)+1, max([x[1] for x in info]), len(info[0][2]), len(info[1][2])))
    fout_del.close()
    fout_dup.close()

    #Make plots for deletions
    del_info = make_plots("%s_dels.txt" % opts.out, opts)
    #Write out info about deletions
    fout = open("%s_delsinfo.txt" % opts.out, 'w')
    fout.write("Ref\t#UniqueChim\t#ReadsChim\tAvgReadsPerChim\tStdReadsPerChim\tPropMappedChim\n")
    for ref, info in del_info.iteritems():
        fout.write("%s\t%d\t%d\t%.4f\t%.4f\t%.4f\n" % (ref, len(info.keys()), sum(info.values()), np.mean(info.values()), np.std(info.values()), sum(info.values())/len(all_aligned[ref])))
    fout.close()
    
    
#    dup_info = make_plots("%s_dups.txt" % opts.out, opts)


###-----------------End of main()--------------------------->>>

def missing_bases(pairs1, pairs2):
    ref1, ref2 = [x[1] for x in pairs1], [x[1] for x in pairs2]
    full=range(min(ref1+ref2), max(ref1+ref2)+1)
    covered = get_covered(pairs1, pairs2)
    missing = list(set(full).difference(set(covered)))
    return missing

#Currently using this to correct for query bases that are multiply mapped to the reference
#Only the left most mapping is included in the output
def get_covered(pairs1, pairs2):
    overlap = set([x[0] for x in pairs1]).intersection(set([x[0] for x in pairs2]))
    if overlap:
        d1 = dict(pairs1)
        d2 = dict(pairs2)
        for q, r in d2.iteritems():
            if q not in d1: d1[q]=r
            elif r<d1[q]: d1[q]=r
        return sorted(d1.values())

    else: return sorted(list(set([x[1] for x in pairs1]+[x[1] for x in pairs2])))

def common_bases(cov1, cov2):
    return list(set(cov1).intersection(set(cov2)))


def corr_hardclipped(read):
    length = read.infer_query_length()
    for t, l in read.cigartuples:
        if t==5: length+=l
    pairs = read.get_aligned_pairs(matches_only=True)
    if read.cigartuples[0][0]==5:
        corr_pairs=[]
        for q, r in pairs:
            corr_pairs.append((q+read.cigartuples[0][1],r))
        return corr_pairs, length
    return pairs,length

def get_aligned(tups, index):
    combo=[]
    for each in tups:
        combo+=[x[index] for x in each]
    return list(set(combo)), combo


###---For plotting

def make_plots(info_file, opts):
    
    plot_info = {}

# Parse info from file
    linecount=0
    for line in open(info_file, 'r'):
        linecount+=1
        if linecount>1:
            cols = line.strip().split('\t')
            #Check to see if reference name is already a key in the dictionary
            if cols[2] not in plot_info: plot_info[cols[2]]={}
            #Check to see if the exact same deletion is already a key in a sub-dictionary, add a count
            if (cols[4],cols[5]) not in plot_info[cols[2]]: plot_info[cols[2]][(cols[4],cols[5])]=1
            else: plot_info[cols[2]][(cols[4],cols[5])]+=1

    for ref, counts in plot_info.iteritems():
        freq_plot(ref, counts, info_file, opts)
    
    #Return info dict
    return plot_info

def freq_plot(ref, counts, infofile, opts):

    #Prep data for plot
    start = []
    stop = []
    freq = []


    for pos, f in counts.iteritems():
        start.append(int(pos[0]))
        stop.append(int(pos[1]))
        freq.append(f)

    #To color points by density
    xy = np.vstack([stop,start])
#    print xy
    z = gaussian_kde(xy)(xy)
    #Sort the points by density so that the densest points are plotted last
    idx = z.argsort()
    stop, start, freq, z = np.asarray(stop)[idx], np.asarray(start)[idx], np.asarray(freq)[idx], z[idx]


    fig, ax = plt.subplots()
    cax = ax.scatter(stop, start, s=freq*10, c=z, alpha=0.5)
    cbar = fig.colorbar(cax, ticks=[min(z), (min(z)+max(z))/2, max(z)])
    cbar.ax.set_yticklabels(['Low', 'Medium', 'High'])
#    lgnd = fig.legend(handles = [cax], labels = ['test'], loc="upper left", scatterpoints=2, columnspacing=20)
    
    ax.set_xlabel('End')
    ax.set_ylabel('Start')
#    ax1.axhline(y=0, xmin=0, xmax=10000, color='k')
    ax.plot([0, max(stop)], [0, max(stop)], ls='--', c='0.3')
    ax.set_xlim([0,max(stop)])
    ax.set_ylim([0,max(stop)])
    ## plot where the genes are
#    if opts.orfs:
#        for x,i in enumerate(opts.orfs):
#            w=10
#            y=100 + i[-1] * w 
##            print i[0],y,i[1]-i[0]
#            #print (i[1]-i[0])%3,i[2]
#            ax1.arrow(i[0],y,i[1]-i[0],0.0,alpha=0.4,head_width=w,width=w,head_length=0,length_includes_head=True,facecolor='gray',edgecolor='k',zorder=1)
#            ax1.text(np.mean([i[1],i[0]]),y,'%s'%(i[2]),va='center',ha='center',size=10,zorder=2)
#    ax1.set_xticks(x)
#    ax1.set_xticklabels(names)

#For legend
    gll = plt.scatter([],[], s=1*10, marker='o', color='#555555', alpha=0.5)
    gl = plt.scatter([],[], s=mid4leg(max(freq))*10, marker='o', color='#555555', alpha=0.5)
    ga = plt.scatter([],[], s=max(freq)*10, marker='o', color='#555555', alpha=0.5)

    fig.legend((gll,gl,ga), (str(1), "%.0f" % mid4leg(max(freq)), str(max(freq))), scatterpoints=1,
       loc=(0.15, 0.70), ncol=1, fontsize=8, labelspacing=2)
    fig.savefig("%s-%s_freqs.pdf" % (infofile, ref))
    fig.clf()

def mid4leg(maxfreq):
    if maxfreq%2==0: return maxfreq/2
    else: return (maxfreq+1)/2



###---May not be using below this

# will cut fasta name off at the first whitespace
def read_fasta_lists_simple_names(file):
    fin = open(file, 'r')
    count=0

    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:].split()[0])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)

    return names, seqs


#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()

def read_fasta_dict_simple_names(file):
    names, seqs = read_fasta_lists_simple_names(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

###------------END of functions used in building new consensus from pileup----------------------------


        
        
###------------->>>

if __name__ == "__main__":
    main()

