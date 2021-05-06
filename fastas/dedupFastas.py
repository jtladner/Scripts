#!/usr/bin/env python

import argparse, glob, os, sys
#import inout as io             #Available at https://github.com/jtladner/Modules
import fastatools as ft        #Available at https://github.com/jtladner/Modules

from collections import defaultdict

# This script removes identical sequences from fasta files
# Identity must be at both the name and sequence level to be removed

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("fastas", help="Fasta files to dedup", nargs='+')
    parser.add_argument("--prepend", default="dedup", help="String to add to the beginning of deduped files.")

    args = parser.parse_args()
    
    
    for each in args.fastas:

        # Read in fasta file
        names, seqs = ft.read_fasta_lists(each)
    
        # Convert to dictionary with keys = names, and values = lists of seqs
        fD = defaultdict(list)
        for i, n in enumerate(names):
            fD[n].append(seqs[i])
        
        newN = []
        newS = []
        
        for n, sL in fD.items():
            if len(set(sL)) == 1:
                newN.append(n)
                newS.append(sL[0])
            else:
                print(n)
                for s in sL:
                    newN.append(n)
                    newS.append(s)
        
        ft.write_fasta(newN, newS, "%s_%s" % (args.prepend, each))
    
#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

