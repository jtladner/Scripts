#!/usr/bin/env python

from __future__ import division
import optparse, os

#This script uses an aligned fasta file and a tab deliminted file containing CDS coordinates to create a nexus input for BEAST

#In version 2.0, added a flag to throw if you only want coding sequence to be included in the nexus file

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-f', '--fasta',  help='Aligned fasta. [None]')
    p.add_option('-c', '--coords', help='Tab delimited file with coordinates of CDS. Should have at least 3 tab delimited columns. The first is not used, will probably have some sort of CDS name. The next two have start and stop base positions.[None]')
    p.add_option('-o', '--out',  help='Name for output nexus file. [None]')
    p.add_option('--onlyCDS', default=False, action="store_true",  help='Use this flag if you only want coding regions to be included in the output nexus file. [None]')
    opts, args = p.parse_args()
    
    make_beast_nexus(opts)
        
#----------------------End of main()

def make_beast_nexus(opts):
    fout=open(opts.out, 'w')

    #Read in seqs
    names, seqs = read_fasta_lists(opts.fasta)
    #Get coding coordinates
    coding_coords=get_coords(opts.coords)
    
    #Make concatenated coding seqs
    coding_seqs=['']*len(seqs)
    for start, end in coding_coords:
        for i in range(len(seqs)):
            coding_seqs[i]+=seqs[i][start-1:end]
    
    
    if opts.onlyCDS:
        fout.write("#NEXUS\n[File created using beast_nexus_prep.py using %s and %s]\n\nBEGIN TAXA;\n" % (opts.fasta, opts.coords))
        fout.write("DIMENSIONS NTAX=%d;\n\nTAXLABELS\n%s\n;\n\nEND;\n" % (len(names), '\n'.join(names)))
        fout.write("BEGIN CHARACTERS;\nDIMENSIONS NCHAR=%d;\nFORMAT DATATYPE=DNA MISSING=N GAP=-;\nMATRIX\n\n%s\n;\n\nEND;\n\n" % (len(coding_seqs[0]), '\n'.join(['%s %s' % (names[x], coding_seqs[x]) for x in range(len(names))])))    
        fout.write("BEGIN ASSUMPTIONS;\n\tcharset coding = 1-%d;\nend;\n" % (len(coding_seqs[0])))

    else:
        #Get non-coding coordinates
        noncoding_coords=extrap_noncoding(coding_coords, len(seqs[0]))
    
        #Make concatenated noncoding seqs
        noncoding_seqs=['']*len(seqs)
        for start, end in noncoding_coords:
            for i in range(len(seqs)):
                noncoding_seqs[i]+=seqs[i][start-1:end]
    
        concat_seqs=[coding_seqs[i]+noncoding_seqs[i] for i in range(len(seqs))]
    
        coding_start=1
        coding_end=len(coding_seqs[0])
        noncoding_start=coding_end+1
        noncoding_end=len(concat_seqs[0])
    
        fout.write("#NEXUS\n[File created using beast_nexus_prep.py using %s and %s]\n\nBEGIN TAXA;\n" % (opts.fasta, opts.coords))
        fout.write("DIMENSIONS NTAX=%d;\n\nTAXLABELS\n%s\n;\n\nEND;\n" % (len(names), '\n'.join(names)))
        fout.write("BEGIN CHARACTERS;\nDIMENSIONS NCHAR=%d;\nFORMAT DATATYPE=DNA MISSING=N GAP=-;\nMATRIX\n\n%s\n;\n\nEND;\n\n" % (len(concat_seqs[0]), '\n'.join(['%s %s' % (names[x], concat_seqs[x]) for x in range(len(names))])))    
        fout.write("BEGIN ASSUMPTIONS;\n\tcharset coding = %d-%d;\n\tcharset noncoding = %d-%d;\nend;\n" % (coding_start, coding_end, noncoding_start, noncoding_end ))



    fout.close()
    
def extrap_noncoding(coding_coords, seq_len):
    non_coords=[]
    #To handle noncoding at the very beginning of the sequence
    if coding_coords[0][0] != 1:
        non_coords.append((1,coding_coords[0][0]-1))    
        
    #To handle noncoding regions in between coding seqs
    coding_sorted=sorted(coding_coords[:])
    for i in range(len(coding_sorted[:-1])):
        if coding_sorted[i+1][0]-coding_sorted[i][1]>0:
            non_coords.append((coding_sorted[i][1]+1,coding_sorted[i+1][0]-1))

    #To handle non-coding at the very end of the sequence
    if coding_coords[-1][1] != seq_len:
        non_coords.append((coding_coords[-1][1]+1, seq_len))

    print non_coords
    return non_coords

def get_coords(c_file):
    fin=open(c_file, 'r')
    coords=[]
    for line in fin:
        cols=line.strip().split('\t')
        coords.append((int(cols[1]), int(cols[2])))
    return coords

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
