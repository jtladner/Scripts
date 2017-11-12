#!/usr/bin/env python

from __future__ import division
import optparse, os
from subprocess import Popen, PIPE
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
#import itertools as it

#This script reads in a vcf file with genotype info and counts the number of parsimony informative SNPs and indels

def main():
    usage = '%prog [options] subset1, subset2'
    p = optparse.OptionParser()
    p.add_option('-v', '--vcf',  help='vcf file. [None, REQ]')
    p.add_option('-w', '--win',  type='int', default=200, help='Window size to check for enrichment. [200]')
    p.add_option('-n', '--num',  type='int', default=4, help='Min number of changes for enrichment. [4]')
#    p.add_option('-o', '--out',  help='base name for output files. [None, REQ]')
    
    opts, args = p.parse_args()
    
    find_hotspots(opts)
###------------------------------

def find_hotspots(opts):
    ct_dict={}
    #Step through the vcf file and get info for T-C transtions
    for line in open(opts.vcf, 'r'):
        if line.startswith('#') and not line.startswith('##'):
            cols=line.strip().split('\t')
            ids=cols[9:]
        elif not line.startswith('#'):
            cols=line.strip().split('\t')
            #If this is a T to C transition
            if cols[3]=='T' and cols[4]=='C':
                c_indices=[i for i,x in enumerate(cols[9:]) if x=='C']
                for i in c_indices:
                    if ids[i] not in ct_dict: ct_dict[ids[i]]=[int(cols[1])]
                    else: ct_dict[ids[i]].append(int(cols[1]))
    #Step through seqs and look for windows of T to C mutation enrichment
    for seq, muts in ct_dict.iteritems():
        clusters=find_clusters(muts, opts)
        if clusters:
            print "%s\t%s" % (seq, "\t".join([",".join([str(y) for y in x]) for x in clusters]))


def find_clusters(mut_list, opts):
    clusters=[]
    mut_list.sort()
    i=0
    while i<=len(mut_list)-opts.num:
        y=opts.num-1
        if (mut_list[i+y]-mut_list[i])<opts.win:
            while (mut_list[i+y]-mut_list[i])<opts.win:
                y+=1
                if i+y >= len(mut_list): break
            clusters.append(mut_list[i:i+y])
        i+=1
    return clusters
                            




#### Not currently using the functions below here

def parse_variant(var_str):
    if var_str.startswith("ANN="):
        var_str=var_str.split("ANN=")[1]
        var_types=[x.split('|')[1] for x in var_str.split(',')]
        var_locs=[x.split('|')[3] for x in var_str.split(',')]
        if len(var_types)==1: return var_types[0], var_locs[0]
        else:
            dom_index=get_dom_index(var_types)
            return var_types[dom_index], var_locs[dom_index]
    else: return "intergenic_region", "terminus"

def get_dom_index(var_list):
    if 'stop_gained' in var_list: return var_list.index('stop_gained')
    if 'frameshift_variant' in var_list: return var_list.index('frameshift_variant')
    elif 'missense_variant' in var_list: return var_list.index('missense_variant')
    elif 'synonymous_variant' in var_list: return var_list.index('synonymous_variant')
    else: print var_list
    
def read_fasta_dict(file):
    names, seqs = read_fasta_lists(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict
    
# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
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

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

