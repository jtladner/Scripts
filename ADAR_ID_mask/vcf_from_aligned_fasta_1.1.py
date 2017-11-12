#!/usr/bin/env python

from __future__ import division
import optparse, os, time
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#This script reads in an aligned fasta file and reports sequence-level and alignment-level stats regarding variable sites and missing data
#Version 2.0 added support for stats on seqs with just one fasta
#Version 2.1 added support for unaligned fastas

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-f', '--fasta',  help='Aligned fasta. [None]')
    p.add_option('-r', '--ref',  help='Name of sequence in the alignment to use as a reference. [None]')
    p.add_option('-o', '--out',  help='Name for output vcf file. [None]')
    p.add_option('-g', '--gap', default='-', help='Character used for gaps. [-]')

    opts, args = p.parse_args()
    
    make_vcf(opts)
        
#----------------------End of main()

def make_vcf(opts):

    #Read in seqs
    names, seqs = read_fasta_lists(opts.fasta)
    refseq=seqs[names.index(opts.ref)]
    #If the sequences are unaligned (not all the same length)
    if len(set([len(x) for x in seqs]))>1:
        print 'Sequences must all be the same length!!!'
        return
                
    else:
        fout=open(opts.out, 'w')
        fout.write('##fileformat=VCFv4.2\n##fileDate=%s\n##source=vcf_from_aligned_fasta_1.0.py\n##reference=%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (time.strftime("%m-%d-%Y"), opts.ref, '\t'.join(names)))
        #Get info for variable positions
        idxs, alts, positions = align_stats(seqs, refseq, opts)
        #Write out info to vcf file 
        for list_idx, var_idx in enumerate(idxs):
            for a in alts[list_idx]:
                fout.write('%s\t%d\t.\t%s\t%s\t.\t.\t.\t.\t%s\n' % (opts.ref, positions[list_idx], refseq[var_idx], a, genos(seqs, var_idx, refseq[var_idx].upper())))
    fout.close()

def genos(seqs, var_idx, ref_base):
    these_genos=[]
    bases=[x[var_idx].upper() for x in seqs]
    for b in bases:
        if b==ref_base: these_genos.append('.')
        else: these_genos.append(b)
    return '\t'.join(these_genos)

def seq_stats_info(s):
    var_i=[]
    for i in range(len(s)):
        if s[i].upper()=='N': var_i.append(i)
    return var_i
        
def perc(num, den):
    return num/den*100

def align_stats(seqs, ref, opts):
    positions=[]
    idxs=[]
    alts=[]
    ref_count=0
    for i in range(len(ref)):
        if ref[i] != opts.gap: ref_count+=1
#        else: print 'Gap at %d' % (i+1)
        uniq_bases=set([x[i].upper() for x in seqs]).difference(set(['N']))
        if len(uniq_bases) > 1: 
            idxs.append(i)
            alts.append(list(uniq_bases.difference(set([ref[i].upper()]))))
            positions.append(ref_count)
    return idxs, alts, positions
    
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

