#!/usr/bin/env python

# This script reads a fasta file with one or more sequences
# And outputs a fasta containing all the same seqs, but truncated to a specific portion of the input seqs

import argparse, os
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import inout as io             #Available at https://github.com/jtladner/Modules

def main():

    arg_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument( '-s', '--startCoord', default = 0, type=int, help = "0-indexed coordinate for truncated alignment to begin.")
    arg_parser.add_argument( '-e', '--endCoord', default = 1000000000, type=int, help = "0-indexed coordinate for truncated alignment to end (inclusive).")

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-i', '--input', help = "Fasta file to be truncated", required=True )
    reqArgs.add_argument( '-o', '--out', help = "Name for output, truncated fasta file", required=True )

    args = arg_parser.parse_args()
    
    # Read in fasta alignment
    fD = ft.read_fasta_dict_upper(args.input)
    
    # Truncate seqences
    outD = {n:s[args.startCoord:args.endCoord+1] for n,s in fD.items()}
    
    #Write out fasta file
    ft.write_fasta_dict(outD, args.out)

###--- End of main() ------->>>


###--- End of functions ------->>>

if __name__ == '__main__':
    main()