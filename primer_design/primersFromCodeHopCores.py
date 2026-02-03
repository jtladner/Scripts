#! /usr/bin/env python3

import argparse, os

import fastatools as ft	 # Available at: https://github.com/jtladner/Modules
import seqtools as st	   # Available at: https://github.com/jtladner/Modules
import inout as io		  # Available at: https://github.com/jtladner/Modules

import numpy as np
import pandas as pd
from collections import defaultdict
import primer3



def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-n', '--ntFasta', help='Fasta file with aligned target nucleotide sequences. The alignment coordinated should correspond perfectly to those of teh aa-level alignment given to j-codehop, but it does not have to include all of the same sequences considered in the codehop core selection. You can focus here on a subset of high priority samples for the design of the clamp region.')
	reqArgs.add_argument('-i', '--inName',  help='Input tsv file with info for the codehop core regions selected by j-codehop, or similar. This file should contain at least 4 named columns: Core (core nt sequence), CoreStart (1-based, inclusive), CoreEnd (1-based, inclusive), Direction (forward or reverse). A Name column can be optionally provided. Column order does not matter, but headers must be exactly as specified', required=True)
	reqArgs.add_argument('-o', '--outName',  help='Input tsv file with info for the codehop core regions selected by j-codehop, or similar. This file should contain at least 4 named columns: Core (core nt sequence), CoreStart (1-based, inclusive), CoreEnd (1-based, inclusive), Direction (forward or reverse). A Name column can be optionally provided. Column order does not matter, but headers must be exactly as specified', required=True)

	parser.add_argument('-m', '--meltTemp', default=65, type=float, help='Target melting temperature for full primer, core + clamp.')
	parser.add_argument('-c', '--minClamp', default=8, type=int, help='Minimum clamp length.')

	args = parser.parse_args()
	
	# Read in nucleotide alignment
	fD = ft.read_fasta_dict(args.ntFasta)
		
	# Read in core info
	inDF = pd.read_csv(args.inName, sep="\t", header=0, keep_default_na=False, low_memory=False)
	
	# Check for expected column headers
	if "Core" not in inDF.columns:
		print("File provided with -i flag must contain a column labeled 'Core'. Check your input file.")
		return None
	elif "CoreStart" not in inDF.columns:
		print("File provided with -i flag must contain a column labeled 'CoreStart'. Check your input file.")
		return None
	elif "CoreEnd" not in inDF.columns:
		print("File provided with -i flag must contain a column labeled 'CoreEnd'. Check your input file.")
		return None
	elif "Direction" not in inDF.columns:
		print("File provided with -i flag must contain a column labeled 'Direction'. Check your input file.")
		return None
	
	# Open output file for writing and add header
	fout = open(args.outName, "w")
	fout.write("Name\tLength\tCoreLength\tClampLength\tDegeneracy\tDirection\tStart\tEnd\tSequence\tRevComp\tAvgTm\tMinTm\tMaxTm\tAvgCoreMis\tMinCoreMis\tMaxCoreMis\tAvgClampMis\tMinClampMis\tMaxClampMis\n")
	
	# Step through each core sequence
	for i,row in inDF.iterrows():
		
		# Check that you are finding matches to the provided core sequence at the provided location
		# This is a good sanity check. Though, you won't necessarily see a perfect match to all sequences
		coreSeq = row["Core"].upper()
		if row["Direction"].strip().upper() == "FORWARD" or row["Direction"].strip().upper() == "F":
			numMatched,coreMis = checkCoreMatches(coreSeq, {n:s[row["CoreStart"]-1:row["CoreEnd"]] for n,s in fD.items()})
			i=row["CoreStart"]-1
			allCores = st.expand_degenerate_seq(coreSeq)
			direction="F"
		elif row["Direction"].strip().upper() == "REVERSE" or row["Direction"].strip().upper() == "R":
			revComp = st.revCompDNA(row["Core"])
			numMatched,coreMis = checkCoreMatches(revComp, {n:s[row["CoreStart"]-1:row["CoreEnd"]] for n,s in fD.items()})
			i=row["CoreEnd"]
			allCores = st.expand_degenerate_seq(revComp)
			direction="R"
		
		print(f"{numMatched}/{len(fD)} sequences perfectly match core: {row['Core'].upper()}")
		
		primer, clamp, mismatches = designClamp(coreSeq, i, direction, fD, args.meltTemp, clampSize=args.minClamp)
		allPrimers = st.expand_degenerate_seq(primer)
		tmL = [primer3.calc_tm(ps) for ps in allPrimers]
		
		if "Name" not in inDF.columns:
			primerName = row["Name"]
		else:
			primerName = f"{st.translate(allCores[0])}-{direction}"
		
		if direction=="F":
			start=row["CoreEnd"]-len(primer)+1
			end=row["CoreEnd"]
		else:
			start=row["CoreStart"]
			end=row["CoreStart"]+len(primer)-1
		
		fout.write(f"{primerName}\t{len(primer)}\t{len(coreSeq)}\t{len(clamp)}\t{len(allPrimers)}\t{direction}\t{start}\t{end}\t{primer}\t{st.revCompDNA(primer)}\t{np.mean(tmL):.2f}\t{min(tmL):.2f}\t{max(tmL):.2f}\t{np.mean(coreMis):.2f}\t{min(coreMis):.2f}\t{max(coreMis):.2f}\t{np.mean(mismatches):.2f}\t{min(mismatches)}\t{max(mismatches)}\n")


#####-------------------->>>

def checkCoreMatches(coreSeq, fD):
	expandedSeqs = {k:"" for k in st.expand_degenerate_seq(coreSeq)}
	mismatches=[]
	perfects=0
	for s in fD.values():
		if s in expandedSeqs:
			mismatches.append(0)
			perfects+=1
		else:
			l=[]
			for every in expandedSeqs:
				l.append(sum([1 for i in range(len(every)) if s[i] != every[i]]))
			mismatches.append(min(l))
	
	return perfects, mismatches

def avgTm(ambigSeq):
	tmL = [primer3.calc_tm(ps) for ps in st.expand_degenerate_seq(ambigSeq)]
	return np.mean(tmL)

# clampSize provided to function is just the starting size. Final size will be determined by target Tm
def designClamp(coreSeq, i, direction, fD, meltTemp, clampSize=10):
	if direction=="F":
		availD = {n:s[:i].replace("-", "") for n,s in fD.items()}
	elif direction=="R":
		availD = {n:st.revCompDNA(s[i:].replace("-", "")) for n,s in fD.items()}
	else:
		print(f"Direction must be either F or R, {direction} is NOT a valid option.")
		return None

	toAddD = {n:s[-clampSize:] for n,s in availD.items()}
	toAddS = st.consensus(list(toAddD.values()), noAmbig=True)
	primer = toAddS+coreSeq
	while avgTm(primer)<meltTemp:
		prevLen = len(primer)
		clampSize+=1
		toAddD = {n:s[-clampSize:] for n,s in availD.items()}
		toAddS = st.consensus(list(toAddD.values()), noAmbig=True)
		primer = toAddS+coreSeq
		if not len(primer)>prevLen:
			print("Not enough clamp bases available to meeting melt temp threshold.")
			break
	
	# For final clamp, count the number of mismatches between the clamp and each sequence
	mismatches=[]
	for n,s in toAddD.items():
		mis = sum([1 for i in range(len(toAddS)) if s[i] != toAddS[i]])
		mismatches.append(mis)
	
	return primer, toAddS, mismatches


if __name__ == '__main__':
	main()
