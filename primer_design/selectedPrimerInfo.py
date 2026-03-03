#! /usr/bin/env python3

import argparse, os

import fastatools as ft	 # Available at: https://github.com/jtladner/Modules
import seqtools as st	   # Available at: https://github.com/jtladner/Modules
import inout as io		  # Available at: https://github.com/jtladner/Modules
from collections import defaultdict

# import numpy as np
# import pandas as pd
#import primer3


# The purpose of this script is to parse an output file from j-codehop (exported primers) and prepare it for input into primersFromCodeHopCores.py

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-f', '--fasta',  help='Aligned fasta file containing target seqs of interest as well as selected primer sequences.', required=True)
	reqArgs.add_argument('-p', '--pairs',  help='Tab-delimited text file with selected primer pairs (one per line). One column should have the header "Forward" and another "Reverse". Names of the primers provided in this file should match the names of the primer sequences in the provided input fasta file.', required=True)
	reqArgs.add_argument('-o', '--outName',  help='Base name for the output files.', required=True)

	args = parser.parse_args()

	# Read in aligned fasta file
	fD = ft.read_fasta_dict(args.fasta, simple=True)
	
	# Read in primer pairs
	pD = io.fileDictHeaderLists(args.pairs, "Forward", "Reverse")
	
	# Extract all primer sequence names
	allPrimerNames = []
	for k, vL in pD.items():
		allPrimerNames.append(k)
		allPrimerNames+=vL
	allPrimerNames = sorted(list(set(allPrimerNames)))
# 	print(allPrimerNames)
	
	# Initiate dictionary to hold masked sequences
	oD = {}
	# Initiate dictionary to hold prop match info for each primer and target
	mD=defaultdict(dict)
	
	# To store target sequence names
	targetNames=[]
	
	# Step through each of the primer sequences
	for pn in allPrimerNames:
		ps = fD[pn]
		oD[pn] = ps.replace("-", "")
		# Initiate a dictionary to link primer positions to primer bases
		psD={}
		for i, b in enumerate(ps):
			if b != "-":
				psD[i]=b
		
		# Step through all input sequences
		for n,s in fD.items():
			# If this is a target sequence
			if n not in allPrimerNames:
				targetNames.append(n)
				oD[f"{n} {pn}"]=""
				matches=0
				for pos, primeBase in psD.items():
					if s[pos] in st.iupacStringD[primeBase]:
						oD[f"{n} {pn}"]+="."
						matches+=1
					else:
						oD[f"{n} {pn}"]+=s[pos]
				mD[pn][n] = (matches, matches/len(psD))

	# Write out file with masked sequences
	ft.write_fasta_dict(oD, f"{args.outName}_maskedPrimers.fasta")
	
	# Simplify targets to remove repeats
	targetNames = sorted(list(set(targetNames)))
	
	# Write out primer centric files with info on # and proportion of matching bases
	with open(f"{args.outName}_propMatch.tsv", "w") as fout_prop:
		with open(f"{args.outName}_numMatch.tsv", "w") as fout_num:
			pnS = "\t".join(mD.keys())
			fout_prop.write(f"Sequence\t{pnS}\n")
			fout_num.write(f"Sequence\t{pnS}\n")
			for n in targetNames:
				propS = "\t".join([f"{mD[pn][n][1]:.3f}" for pn in mD])
				fout_prop.write(f"{n}\t{propS}\n")

				numS = "\t".join([str(mD[pn][n][0]) for pn in mD])
				fout_num.write(f"{n}\t{numS}\n")

	# Write out file focused on primer pairs
	with open(f"{args.outName}_propMatch_primerPairs.tsv", "w") as fout:
		tnS = "\t\t".join(targetNames)
		fout.write(f"F_primer\tR_primer\t{tnS}\n")
		for f,rL in pD.items():
			for r in rL:
				fout.write(f"{f}\t{r}\t")
				for tn in targetNames:
					fout.write(f"{mD[f][tn][1]:.3f}\t{mD[r][tn][1]:.3f}\t")
				fout.write("\n")
	


#####-------------------->>>


if __name__ == '__main__':
	main()
