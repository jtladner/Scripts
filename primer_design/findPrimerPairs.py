#! /usr/bin/env python3

import argparse, os
from pathlib import Path

import fastatools as ft	 # Available at: https://github.com/jtladner/Modules
import seqtools as st	   # Available at: https://github.com/jtladner/Modules
import inout as io		  # Available at: https://github.com/jtladner/Modules

import numpy as np
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
	parser.add_argument('-m', '--minAmpLen', default=300, type=float, help='Minimum for average amplicon length, for a pair to be considered.')
	parser.add_argument('-x', '--maxAmpLen', default=700, type=float, help='Maximum for average amplicon length, for a pair to be considered.')
	parser.add_argument('--meltTempDiff', default=5, type=float, help='Maximum average melt temp diff for forward and reverse primers.')
	parser.add_argument('--structureThresh', default=45, type=float, help='Maximum threshold for heterodimer Tm.')
	parser.add_argument('--max2report', default=5, type=int, help='Maximum number of primers to include in output, per input core.')
	parser.add_argument('--minProp', default=0.2, type=float, help='Minimum proportion of base in reference sequences to be considered for clamp inclusion.')
	parser.add_argument('--max2consider', default=100000, type=int, help='Maximum full length clamp sequences to consdier. --maxClamp will be dynamically adjusted to meet this amount.')
	

	args = parser.parse_args()

	# Construct log file name, if not provided
# 	if not args.log:
# 		dirName = os.path.dirname(args.outName)
# 		base_noExt = Path(args.outName).stem
# 		args.log = f"{dirName}/{base_noExt}_Log.tsv"
	
	# Read in nucleotide alignment
	fD = ft.read_fasta_dict(args.ntFasta)
		
	# Read in individual primer info
	inDF = pd.read_csv(args.inName, sep="\t", header=0, keep_default_na=False, low_memory=False)
	
	# Check for expected column headers
	for headName in ["Name", "Direction", "Start", "End", "Sequence", "AvgTm"]:
		if headName not in inDF.columns:
			print(f"File provided with -i flag must contain a column labeled '{headName}'. Check your input file.")
			return None
	
	# Open output file for writing and add header
	fout = open(args.outName, "w")
	fout.write("F_name\tR_name\tAvgAmpLength\tF_Degeneracy\tR_Degeneracy\tF_Start\tF_End\tR_Start\tR_End\tF_seq\tR_seq\tF_RevComp\tR_RevComp\tAvgHeterodimerTm\tMinHeterodimerTm\tMaxHeterodimerTm\tF_AvgCoreMis\t\tR_AvgCoreMis\tF_AvgClampMis\tR_AvgClampMis\tF_AvgTm\tR_AvgTm\tF_AvgHomoDimTm\tR_AvgHomoDimTm\tF_AvgHairpinTm\tR_AvgHairpinTm\n")

	
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
		
		propMatched = numMatched/len(fD)
		print(f"{numMatched}/{len(fD)} sequences perfectly match core: {row['Core'].upper()}")
		
		# Determine name to be used in outputs
		if "Name" not in inDF.columns:
			primerName = f"{st.translate(allCores[0])}-{direction}"
		else:
			primerName = row["Name"]
		
		
		if args.simpleConsensus:
			primerL, clampL, mismatchesL = designClamp(coreSeq, i, direction, fD, args.meltTemp, clampSize=args.minClamp)
		else:
			primerL, clampL, mismatchesL, logD = designClamps(coreSeq, i, direction, fD, args)
			flog.write(f"{primerName}\t{propMatched:.3f}\t{logD['maxLen']}\t{logD['start']}\t{logD['gc']}\t{logD['tm']}\t{logD['hairpin']}\t{logD['homodimer']}\n")


		
		for num,primer in enumerate(primerL):
			allPrimers = st.expand_degenerate_seq(primer)
			tmL = [primer3.calc_tm(ps) for ps in allPrimers]
			
			if direction=="F":
				start=row["CoreEnd"]-len(primer)+1
				end=row["CoreEnd"]
			else:
				start=row["CoreStart"]
				end=row["CoreStart"]+len(primer)-1
			
			# Calculate hairpin Tm
			hpL = [primer3.calc_hairpin_tm(ps) for ps in allPrimers]
			# Calculate homodimer Tm
			hmL = [primer3.calc_homodimer_tm(ps) for ps in allPrimers]
			
			fout.write(f"{primerName}\t{len(primer)}\t{len(coreSeq)}\t{len(clampL[num])}\t{len(allPrimers)}\t{direction}\t{start}\t{end}\t{np.mean(hpL):.2f}\t{min(hpL):.2f}\t{max(hpL):.2f}\t{np.mean(hmL):.2f}\t{min(hmL):.2f}\t{max(hmL):.2f}\t{primer}\t{st.revCompDNA(primer)}\t{np.mean(tmL):.2f}\t{min(tmL):.2f}\t{max(tmL):.2f}\t{np.mean(coreMis):.2f}\t{min(coreMis):.2f}\t{max(coreMis):.2f}\t{np.mean(mismatchesL[num]):.2f}\t{min(mismatchesL[num])}\t{max(mismatchesL[num])}\n")

	fout.close()
	flog.close()
	
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

def avgGC(ambigSeq):
	gcL = [st.calcGC(ps) for ps in st.expand_degenerate_seq(ambigSeq)]
	return np.mean(gcL)

def maxHairpin(ambigSeq):
	hpL = [primer3.calc_hairpin_tm(ps) for ps in st.expand_degenerate_seq(ambigSeq)]
	return max(hpL)

def maxHomodimer(ambigSeq):
	hmL = [primer3.calc_homodimer_tm(ps) for ps in st.expand_degenerate_seq(ambigSeq)]
	return max(hmL)


# clampSize provided to function is just the starting size. Final size will be determined by target Tm
# Just using a simple consensus approach
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
	
	return [primer], [toAddS], [mismatches]

# This function will generate and consider multiple clamps, attempting to select one with optimal characteristics
def designClamps(coreSeq, i, direction, fD, args):
	logD = {k:0 for k in ["maxLen", "start", "gc", "tm", "hairpin", "homodimer"]}
	
	if direction=="F":
		availD = {n:s[:i].replace("-", "") for n,s in fD.items()}
	elif direction=="R":
		availD = {n:st.revCompDNA(s[i:].replace("-", "")) for n,s in fD.items()}
	else:
		print(f"Direction must be either F or R, {direction} is NOT a valid option.")
		return None

	toAddD = {n:s[-args.maxClamp:] for n,s in availD.items()}
	toAddS_ambig = st.consensus_minProp(list(toAddD.values()), args.minProp)
	toAddSL = st.expand_degenerate_seq(toAddS_ambig)
	if len(toAddSL)>args.max2consider:
		newMax = args.maxClamp
		while len(toAddSL)>args.max2consider:
			newMax-=1
			toAddD = {n:s[-newMax:] for n,s in availD.items()}
			toAddS_ambig = st.consensus_minProp(list(toAddD.values()), args.minProp)
			toAddSL = st.expand_degenerate_seq(toAddS_ambig)
		#print(f"Max clamp length adjusted to {newMax}")
		logD["maxLen"] = newMax
	else:
		logD["maxLen"] = args.maxClamp
	
#	print(f"Considering {len(toAddSL)} clamp sequences with max clamp length.")
	for cls in toAddSL[:]:
		while len(cls)>args.minClamp:
			cls=cls[1:]
			toAddSL.append(cls)
	toAddSL = list(set(toAddSL))
# 	print(f"Considering {len(toAddSL)} clamp sequences. Avg Length = {avglen(toAddSL)}.")
	logD["start"] = len(toAddSL)
	
	
	# GC content filter
	primerL = [cl+coreSeq for cl in toAddSL]
	gcL = [avgGC(ps) for ps in primerL]
	passGC = [toAddSL[i] for i,gc in enumerate(gcL) if args.gcMin<=gc<=args.gcMax]
	if (len(passGC)/len(toAddSL))<0.05:
		print(f"Skipping GC filter because losing more than 95% of clamp options. Average clamp GC content = {np.mean(gcL)}.")
	else:
		toAddSL = passGC[:]
		logD["gc"] = len(toAddSL)
# 	print(f"{len(toAddSL)} sequences within GC target range. Avg Length = {avglen(toAddSL)}.")

	# Tm filter
	primerL = [cl+coreSeq for cl in toAddSL]
	tmL = [avgTm(ps) for ps in primerL]
	passTM = [toAddSL[i] for i,tm in enumerate(tmL) if args.meltTemp<=tm<=args.maxMeltTemp]
	if len(passTM)==0:
		print(f"Skipping Tm filter because The number of clamp sequences because none of the options meet specified melting temperature range. You may want to consider adjusting your clamp size or target melting temperature.")
	else:
		toAddSL = passTM[:]
		logD["tm"] = len(toAddSL)
# 	print(f"{len(toAddSL)} sequences passed Tm check. Avg Length = {avglen(toAddSL)}.")
	
	# Hairpin filter
	primerL = [cl+coreSeq for cl in toAddSL]
	hpL = [maxHairpin(ps) for ps in primerL]
	passHP = [toAddSL[i] for i,tm in enumerate(hpL) if tm<=args.structureThresh]
	if len(passHP)==0:
		print(f"Skipping hairpin filter because none of the options are below specified threshold.")
	else:
		toAddSL = passHP[:]
		logD["hairpin"] = len(toAddSL)
# 	print(f"{len(toAddSL)} sequences passed Hairpin check. Avg Length = {avglen(toAddSL)}.")
	
	# Homodimer filter
	primerL = [cl+coreSeq for cl in toAddSL]
	hmL = [maxHomodimer(ps) for ps in primerL]
	passHM = [toAddSL[i] for i,tm in enumerate(hmL) if tm<=args.structureThresh]
	if len(passHM)==0:
		print(f"Skipping homodimer filter because none of the options are below specified threshold.")
	else:
		toAddSL = passHM[:]
		logD["homodimer"] = len(toAddSL)
# 	print(f"{len(toAddSL)} sequences passed Homodimer check. Avg Length = {avglen(toAddSL)}.")
	
	# For filtered clamp options, count the number of mismatches between the clamp and each sequence
	maxMM=defaultdict(list)
	avgMM=defaultdict(list)
	mmD={}
	ammD={}
	for cs in toAddSL:
		mismatches=[]
		for n,s in toAddD.items():
			mis = sum([1 for i in range(len(cs)) if s[-i] != cs[-i]])
			mismatches.append(mis)
		mmD[cs] = mismatches
		maxMM[max(mismatches)].append(cs)
		avgMM[np.mean(mismatches)].append(cs)
		ammD[cs] = np.mean(mismatches)
	
	# Initial selection based on the maximum number of mismatches between the clamp and the targets
	minMaxMM = min(maxMM.keys())
	output = maxMM[minMaxMM]
	
	# If a focus on max mismatches alone returns more options than the user requested
	if len(output)>args.max2report:
		sortedSeqs = sorted([(ammD[cs],cs) for cs in output])
		output = [x[1] for x in sortedSeqs[:args.max2report]]
	# If a focus on max mismatches alone returns fewer options than the user requested
	elif len(output)<args.max2report:
		extra=[]
		thisMax = minMaxMM
		needed=args.max2report-len(output)
		while len(extra)<needed:
			thisMax+=1
			for k,v in maxMM.items():
				if k<=thisMax and k != minMaxMM:
					extra+=v
		if len(extra)>needed:
			sortedSeqs = sorted([(ammD[cs],cs) for cs in extra])
			extra = [x[1] for x in sortedSeqs[:needed]]

		output+=extra 

	primerL = [cl+coreSeq for cl in output]

	return primerL, output, [mmD[cs] for cs in output], logD

def avglen(sL):
	lenL = [len(s) for s in sL]
	return np.mean(lenL)


if __name__ == '__main__':
	main()
