#! /usr/bin/env python3

import argparse, os
from pathlib import Path

import fastatools as ft	 # Available at: https://github.com/jtladner/Modules
import seqtools as st	   # Available at: https://github.com/jtladner/Modules
import inout as io		  # Available at: https://github.com/jtladner/Modules

import numpy as np
import itertools as it
# import scipy.stats as stats
import pandas as pd
from collections import defaultdict
import primer3

# This scripts parses the output from primersFromCodeHopCores.py and finds pairs of primers that match specified criteria

def main():

	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	reqArgs = parser.add_argument_group('required arguments')
	# Note that these arguments are added directly to the new argument group "reqArgs", not "parser" 
	reqArgs.add_argument('-n', '--ntFasta', help='Fasta file with aligned target nucleotide sequences. The alignment coordinates should correspond perfectly to those of teh aa-level alignment given to j-codehop, but it does not have to include all of the same sequences considered in the codehop core selection. You can focus here on a subset of high priority samples for the design of the clamp region.')
	reqArgs.add_argument('-i', '--inName',  help='Input tsv file with individual forward and reverse primers designed using primersFromCodeHopCores.py', required=True)
	reqArgs.add_argument('-o', '--outName',  help='Name for output file with primer pairs', required=True)

	# Optional arguments
# 	parser.add_argument('-l', '--log', help='Name for log file about the filtering process. If not provided, a name will be constructed from the output file name.')
	parser.add_argument('-m', '--minAmpLen', default=300, type=float, help='Minimum for average amplicon length, for a pair to be considered. Amplicon length is calculated excluding the primers.')
	parser.add_argument('-x', '--maxAmpLen', default=700, type=float, help='Maximum for average amplicon length, for a pair to be considered. Amplicon length is calculated excluding the primers.')
	parser.add_argument('--meltTempDiff', default=5, type=float, help='Maximum average melt temp diff for forward and reverse primers.')
	parser.add_argument('--structureThresh', default=45, type=float, help='Maximum threshold for heterodimer Tm.')
	parser.add_argument('--maxPerAmp', default=10, type=int, help='Maximum number of primer pairs to report per amplicon. The top scoring pairs will be reported.')
	parser.add_argument('--scoreCriteria', help='User can provide a tab-delim file with at least two columns: "Factor" and "Weight". Factor column is used to indicate which criteria from the output primer pair file should be used for scoring. Weight should be a floating point value indicating the relative weight to give each criteria. The higher the weight, the more important the factor will be indetermining the total score. If not specified, default criteria will be used.')
	parser.add_argument('--applyStructureThresh_toAll', default=False, action="store_true", help='Use this flag is you want the structure threshold to be applied retrospectively to the hairpin and homodimer average values calculated previously. This allows the user to apply a stricter threshold, to limit the number of possible primer pairs.')
	
	args = parser.parse_args()

	# Read in scoring criteria, if provided
	if args.scoreCriteria:
		scoreCritD = io.fileDictHeader(args.scoreCriteria, "Factor", "Weight")
	else:
		scoreCritD = {"F_AvgClampMis":1, "R_AvgClampMis":1, "F_AvgHairpinTm":0.5, "R_AvgHairpinTm":0.5, "AvgTmDiff":0.25, "F_AvgHomoDimTm":0.25, "R_AvgHomoDimTm":0.25, "AvgHeterodimerTm":0.25}
		print("Using default scoring criteria\n", scoreCritD)
	
	# Read in nucleotide alignment
	fD = ft.read_fasta_dict(args.ntFasta)
		
	# Read in individual primer info
	inDF = pd.read_csv(args.inName, sep="\t", header=0, keep_default_na=False, low_memory=False)
	
	# Check for expected column headers
	for headName in ["Name", "Direction", "Start", "End", "Sequence", "AvgTm"]:
		if headName not in inDF.columns:
			print(f"File provided with -i flag must contain a column labeled '{headName}'. Check your input file.")
			return None
		
	# Forward primer dictionary
	forPrimerD = defaultdict(list)
	# Reverse primer dictionary
	revPrimerD = defaultdict(list)
	
	# Step through each primer and collect info into two dicts, one for forward primers and one for reverse
	for i,row in inDF.iterrows():
		if row['Direction']=="F":
			forPrimerD[(row['Name'],row['End'])].append(row)
		elif row['Direction']=="R":
			revPrimerD[(row['Name'], row['Start'])].append(row)
	
	# Initiate dictionary to save info for output
	outD = defaultdict(list)
	
	# Step through each forward primer
	for ftup, foptL in forPrimerD.items():
		# Step through each forward primer
		for rtup, roptL in revPrimerD.items():
			
			# Check whether average amplicon length fits specified range
			ampStart = foptL[0]["End"]
			ampEnd = roptL[0]["Start"]
			ampSeqs = [s[ampStart-1:ampEnd].replace("-", "") for s in fD.values()]
			ampLen = [len(s) for s in ampSeqs]
			avgAmpLen = np.mean(ampLen)
			if args.minAmpLen <= avgAmpLen <= args.maxAmpLen:
				
				for eachF in foptL:
					if not args.applyStructureThresh_toAll or max([eachF["AvgHomoDimTm"], eachF["AvgHairpinTm"],]) <= args.structureThresh:
						allPrimersF = st.expand_degenerate_seq(eachF["Sequence"])
						for eachR in roptL:
							if not args.applyStructureThresh_toAll or max([eachR["AvgHomoDimTm"], eachR["AvgHairpinTm"],]) <= args.structureThresh:
								allPrimersR = st.expand_degenerate_seq(eachR["Sequence"])
								hdL = []
								for fs in allPrimersF:
									hdL+=[primer3.calc_heterodimer_tm(fs, rs) for rs in allPrimersR]
								if np.mean(hdL) < args.structureThresh:
									avgTmDiff = abs(eachF['AvgTm']-eachR['AvgTm'])
									if avgTmDiff <= args.meltTempDiff:
										ampID = (ampStart, ampEnd)
										outD[ampID].append({
											"F_name":eachF['Name'], "R_name":eachR['Name'],
											"AvgAmpLength":avgAmpLen, "AvgTmDiff":avgTmDiff,
											"F_Length":eachF['Length'], "F_CoreLength":eachF['CoreLength'], "F_ClampLength":eachF['ClampLength'],
											"R_Length":eachR['Length'], "R_CoreLength":eachR['CoreLength'], "R_ClampLength":eachR['ClampLength'],
											"F_Degeneracy":eachF['Degeneracy'], "R_Degeneracy":eachR['Degeneracy'],
											"F_Start":eachF['Start'], "F_End":eachF['End'],
											"R_Start":eachR['Start'], "R_End":eachR['End'],
											"F_seq":eachF['Sequence'], "R_seq":eachR['Sequence'],
											"F_RevComp":eachF['RevComp'], "R_RevComp":eachR['RevComp'],
											"AvgHeterodimerTm":np.mean(hdL), "MinHeterodimerTm":min(hdL), "MaxHeterodimerTm":max(hdL),
											"F_AvgCoreMis":eachF['AvgCoreMis'], "R_AvgCoreMis":eachR['AvgCoreMis'],
											"F_AvgClampMis":eachF['AvgClampMis'], "R_AvgClampMis":eachR['AvgClampMis'],
											"F_AvgTm":eachF['AvgTm'], "R_AvgTm":eachR['AvgTm'],
											"F_AvgHomoDimTm":eachF['AvgHomoDimTm'], "R_AvgHomoDimTm":eachR['AvgHomoDimTm'],
											"F_AvgHairpinTm":eachF['AvgHairpinTm'], "R_AvgHairpinTm":eachR['AvgHairpinTm']
										})

	# Open output file for writing primer pair details and add header
	fout = open(f"{args.outName}_primerPairs.tsv", "w")
	outCats = ["F_name", "R_name", "Score", "F_AvgClampMis", "R_AvgClampMis", "F_AvgCoreMis", "R_AvgCoreMis", "F_AvgHomoDimTm", "R_AvgHomoDimTm", "F_AvgHairpinTm", "R_AvgHairpinTm", "AvgTmDiff", "AvgHeterodimerTm", "AvgAmpLength", "F_Degeneracy", "R_Degeneracy", "F_Length", "F_CoreLength", "F_ClampLength", "R_Length", "R_CoreLength", "R_ClampLength", "F_Start", "F_End", "R_Start", "R_End", "F_seq", "R_seq", "F_RevComp", "R_RevComp", "F_AvgTm", "R_AvgTm"]
	outCatsStr = "\t".join(outCats)
	fout.write(f"{outCatsStr}\n")

	# Open output file for writing amplicon details and add header
	fout_amp = open(f"{args.outName}_amplicons.tsv", "w")
	outAmpCats = ["F_name", "R_name", "Start", "End", "AvgAmpLength", "AvgPercID", "F_AvgCoreMis", "R_AvgCoreMis", "F_Degeneracy", "R_Degeneracy", "F_CoreLength", "R_CoreLength"]
	outAmpCatsStr = "\t".join(outAmpCats)
	fout_amp.write(f"{outAmpCatsStr}\n")


	# When multiple primer pairs were found per amplicon, generate relative scores based to prioritize among them
	for pairID, infoDL in outD.items():
		
		# Write out summary info for target amplicon (a given pair of F and R cores)
		aD = {c:infoDL[0][c] for c in outAmpCats if c in infoDL[0]}
		aD["Start"] = infoDL[0]["F_End"]
		aD["End"] = infoDL[0]["R_Start"]
		ampSeqs = [s[aD["Start"]-1:aD["End"]] for s in fD.values()]
		percID_L = []
		for s1,s2 in it.combinations(ampSeqs, 2):
			percID_L.append(st.percID_nt(s1, s2))
		aD["AvgPercID"] = np.mean(percID_L)
		outAmpStr = "\t".join([str(aD[c]) for c in outAmpCats])
		fout_amp.write(f"{outAmpStr}\n")
		
		
		if len(infoDL)>1:
			percD={}
			for crit in scoreCritD:
				percD[crit] = {}
				critL = [d[crit] for d in infoDL]
				for val in sorted(list(set(critL))):
#					percD[crit][val] = stats.percentileofscore(critL, val)
					if len(set(critL)) == 1:
						percD[crit][val] = 0
					else:
						percD[crit][val] = (val-min(critL))/(max(critL)-min(critL))
			
			for pD in infoDL:
				score=0
				for crit, vD in percD.items():
					score += vD[pD[crit]] * scoreCritD[crit]
				pD["Score"] = score
			
			# Sort by score and output the top (lowest) scoring pairs 
			sortedL = sorted([(pD["Score"], pD["F_AvgClampMis"], pD["R_AvgClampMis"],i) for i, pD in enumerate(infoDL)])
			for each in sortedL[:args.maxPerAmp]:
				outStr = "\t".join([str(infoDL[each[-1]][c]) for c in outCats])
				fout.write(f"{outStr}\n")


		else:
			for pD in infoDL:
				infoDL[i]["Score"] = score
				outStr = "\t".join([str(pD[c]) for c in outCats])
				fout.write(f"{outStr}\n")

#####-------------------->>>


if __name__ == '__main__':
	main()
