#!/usr/bin/env python3

#This script reads a file that contains info for a number of subsequences that you want to output to a new file
#The new seqs are output to stdout
#Each line of the file specified after the command line should include: fasta_file, seq_name, base_start, base_stop

import argparse, os
import fastatools as ft        #Available at https://github.com/jtladner/Modules
import inout as io             #Available at https://github.com/jtladner/Modules

def main():

    arg_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument('--removeDash', default = False, action="store_true", help = "Use this flag is you want to remove dashes from the sub-sequences after extraction.")

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-i', '--input', help = "Fasta file from which you want to extract subsequences", required=True )
    reqArgs.add_argument( '-t', '--start', type=int, help = "The start position of interest (1-indexed).", required=True)
    reqArgs.add_argument( '-p', '--stop', type=int, help = "The stop position of interest (1-indexed).", required=True)
    reqArgs.add_argument( '-o', '--out', help = "Name for output fasta file", required=True )

    args = arg_parser.parse_args()
    
    # Read in info for seqs to extract
    fastaD = ft.read_fasta_dict(args.input)
    
    # Create dictionary for output fasta file
    outD = {}
    
    for n,s in fastaD.items():
        seqName = f'{n}_{args.start}_{args.stop}'
        seq = s[args.start-1:args.stop]
        
        #Add sequence to output dictionary
        if args.removeDash:
            outD[seqName] = seq.replace("-", "")
        else:
            outD[seqName] = seq
    
    #Write out fasta file
    ft.write_fasta_dict(outD, args.out)

###--- End of main() ------->>>

if __name__ == '__main__':
    main()
