#!/usr/bin/env python

#This script reads a file that contains info for a number of subsequences that you want to output to a new file
#The new seqs are output to stdout
#Each line of the file specified after the command line should include: fasta_file, seq_name, base_start, base_stop

import argparse, os
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import inout as io             #Available at https://github.com/jtladner/Modules

def main():

    arg_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument( '-f', '--fileHeader', default = "File", help = "Header used in input for the column specifying the path to the fasta file containing the sequence of interest.")
    arg_parser.add_argument( '-s', '--seqHeader', default = "SeqName", help = "Header used in input for the column specifying the name of the parent sequence of interest.")
    arg_parser.add_argument( '-t', '--startHeader', default = "Start", help = "Header used in input for the column specifying the start position of interest (1-indexed).")
    arg_parser.add_argument( '-p', '--stopHeader', default = "Stop", help = "Header used in input for the column specifying the stop position of interest (1-indexed).")
    arg_parser.add_argument( '-r', '--revcompHeader', default = "RevComp", help = "Header used in input for the column specifying whether (1) ot not (0) to reverse complement the sequence.")
    arg_parser.add_argument('--removeDash', default = False, action="store_true", help = "Use this flag is you want to remove dashes from the sub-sequences after extraction.")

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-i', '--input', help = "Tab-delimited file. One row should be included per sequence to extract. There should be 5 columns per line: path to fasta file containing sequence, parent sequence name, 1-index start position, 1-index stop position, whether the sequence should be reverse complemented (1) or not (0)", required=True )
    reqArgs.add_argument( '-o', '--out', help = "Name for output fasta file", required=True )

    args = arg_parser.parse_args()
    
    # Read in info for seqs to extract
    infoD = io.fileDictFull(args.input, delim="\t", valType="str")
    
    # Create dictionary for output fasta file
    outD = {}
    
    for i, fn in enumerate(infoD[args.fileHeader]):
        seqName, seq = get_seq(fn, infoD[args.seqHeader][i], infoD[args.startHeader][i], infoD[args.stopHeader][i])
        if infoD[args.revcompHeader][i] == "1":
            seq = rev_comp(seq)
            
        
        #Add sequence to output dictionary
        if args.removeDash:
            outD[seqName] = seq.replace("-", "")
        else:
            outD[seqName] = seq
    
    #Write out fasta file
    ft.write_fasta_dict(outD, args.out)

###--- End of main() ------->>>


def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join([complement.get(base, base) for base in reversed(seq)])
    return reverse_complement


def get_seq(file, seq_name, start, stop):
    seq_dict = ft.read_fasta_dict_upper(file)
    seq_of_interest = seq_dict[seq_name][int(start)-1:int(stop)]
    return '%s_%d_%d %s' % (seq_name, int(start), int(stop), os.path.basename(file)) , seq_of_interest

###--- End of functions ------->>>

if __name__ == '__main__':
    main()