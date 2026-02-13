#! /usr/bin/env python3

import argparse, os

import fastatools as ft	 # Available at: https://github.com/jtladner/Modules
import seqtools as st	   # Available at: https://github.com/jtladner/Modules
import inout as io		  # Available at: https://github.com/jtladner/Modules

import numpy as np
import pandas as pd
from collections import defaultdict
import primer3


# The purpose of this script is to parse an output file from j-codehop (exported primers) and prepare it for input into primersFromCodeHopCores.py

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-i', '--inName',  help='Input csv file with info for primers exported from j-codehop', required=True)
	reqArgs.add_argument('-o', '--outName',  help='Name fo output tsv file that can be used as input for primersFromCodeHopCores.py', required=True)

	# Arguments below here will NOT be used if '--simpleConsensus' flag is used.
# 	parser.add_argument('--maxMeltTemp', default=70, type=float, help='Target melting temperature for full primer, core + clamp.')
# 	parser.add_argument('--maxClamp', default=30, type=int, help='Minimum clamp length.')
# 	parser.add_argument('--gcMin', default=40, type=float, help='Minimum target GC content for full primer, core + clamp.')
# 	parser.add_argument('--gcMax', default=60, type=float, help='Minimum target GC content for full primer, core + clamp.')
	

	args = parser.parse_args()
	
	# Open up output file for writing
	with open(args.outName, "w") as fout:
		fout.write(f"Name\tCore\tDirection\tCoreStart\tCoreEnd\n")

		# Step through input file
		with open(args.inName, "r") as fin:
			lc=0
			for line in fin:
				lc+=1
				cols = line.rstrip("\n").split(", ")
				if lc==1:
					headerD = {head:i for i,head in enumerate(cols)}
				else:
					if len(cols)==len(headerD):
						pName=cols[headerD["Primer Name"]].replace(" ", "-")
						direc=cols[headerD['Direction']]
						coreLen = int(cols[headerD["Core length (NT(AA))"]].split("(")[0])
						coreSeq = cols[headerD["Primer Sequence 5'-3'"]][-coreLen:]
						locL = [int(i) for i in cols[headerD["Primer Location (NT)"]].split("-")]
						if direc=="forward":
							fout.write(f"{pName}\t{coreSeq}\t{direc}\t{locL[1]-coreLen+1}\t{locL[1]}\n")
						elif direc=="reverse":
							fout.write(f"{pName}\t{coreSeq}\t{direc}\t{locL[0]}\t{locL[0]+coreLen-1}\n")
						else:
							print(f"Error: {direc} is NOT a recognized direction, expecting either 'forward' or 'reverse'")
							
					else:
						print(f"Ternimating processing of file with this line: {line}")
						break
		


#####-------------------->>>


if __name__ == '__main__':
	main()
