#!/usr/bin/env python

import argparse, os, random
import fastatools as ft        #Available at https://github.com/jtladner/Modules

from collections import defaultdict


def main():

    arg_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    arg_parser.add_argument( '-r', '--reps', help = "Number of replicate datasets to generate for each level of divergence.", default = 1, type = int )

    reqArgs = arg_parser.add_argument_group('Required Arguments')
    reqArgs.add_argument( '-i', '--input', help = "Fasta file contianing the protein sequence to downsample.", required=True )
    reqArgs.add_argument( '-n', '--num', help = "Size(s) of downsampled datasets. Can be a comma-delimited list of integers", required=True )
    reqArgs.add_argument( '-o', '--output', help = "Directory name for output files. Will be created, if it doesn't already exist", required=True )

    args = arg_parser.parse_args()
    
    # Generate output directory
    if not os.path.isdir(args.output):
        os.mkdir(args.output)
    else:
        print("Warning: %s already exists!" % (args.output))
    
    # Read in fasta file to downsample
    names, seqs = ft.read_fasta_lists(args.input)
    
    # Extract file basename
    bName = ".".join(os.path.basename(args.input).split(".")[:-1])
    
    # Step through each dataset size
    sizes = [int(x) for x in args.num.split(",")]
    
    for s in sizes:
        sCount=0
        while sCount<args.reps:
            indexes = random.choices(range(len(names)), k=s)
            ft.write_fasta([names[i] for i in indexes], [seqs[i] for i in indexes], "%s/%s_n%04d-%03d.fasta" % (args.output, bName, s, sCount))
            sCount+=1
    
#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

