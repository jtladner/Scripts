#!/usr/bin/env python3

import sys

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

###-------->>>

if __name__=='__main__':
    names, seqs = read_fasta_lists(sys.argv[1])
    for index in range(len(seqs)):
        print ("%s\t%d" % (names[index], len(seqs[index])))
