#!/usr/bin/env python

from __future__ import division
import optparse, os

#This script creates a nexus for input to BEAST from an aligned fasta file
#All sites are included in a single charset

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-f', '--fasta',  help='Aligned fasta. [None]')
    p.add_option('-c', '--charset', default='snps', help='Name for charset.[snps]')
    p.add_option('-o', '--out',  help='Name for output nexus file. [None]')
    opts, args = p.parse_args()
    
    make_beast_nexus(opts)
        
#----------------------End of main()

def make_beast_nexus(opts):
    fout=open(opts.out, 'w')

    #Read in seqs
    names, seqs = read_fasta_lists(opts.fasta)
        
    
    fout.write("#NEXUS\n[File created using beast_nexus_prep_simple.py using %s]\n\nBEGIN TAXA;\n" % (opts.fasta))
    fout.write("DIMENSIONS NTAX=%d;\n\nTAXLABELS\n%s\n;\n\nEND;\n" % (len(names), '\n'.join(names)))
    fout.write("BEGIN CHARACTERS;\nDIMENSIONS NCHAR=%d;\nFORMAT DATATYPE=DNA MISSING=N GAP=-;\nMATRIX\n\n%s\n;\n\nEND;\n\n" % (len(seqs[0]), '\n'.join(['%s %s' % (names[x], seqs[x]) for x in range(len(names))])))    
    fout.write("BEGIN ASSUMPTIONS;\n\tcharset %s = 1-%d;\nend;\n" % (opts.charset, len(seqs[0])))


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

###------------------------------------->>>>    

if __name__ == "__main__":
    main()
