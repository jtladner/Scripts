#!/usr/bin/env python

# By Jason Ladner

from __future__ import division
import sys, optparse, os, pysam, itertools

#Identifies individual reads covering 2 positions of interest and records the combined genotypes

#In version 2.0, starting saving read name info in order to handle R1 and R2 in a smarter way, but this could lead to high memory usage for large datasets
    #If both reads cover positions of interest, and they agree about the genotypes, only one is counted
    #If they disagree, neither is counted
    
def main():

    #To parse command line
    usage = "usage: %prog [options] sam1 [sam2 ...]"
    p = optparse.OptionParser(usage)
    
    p.add_option('-i', '--input', help='Text file with info for the position pairs to be examined. One line per pair, tab-delimited with the following columns: pos1, ref1, alt1, pos2, ref2, alt2 [None, OPT]')
    p.add_option('-o', '--out', help='Name for output file, required if using the -i option. [None, OPT/REQ]')
#    p.add_option('-s', '--sam', help='input sam file, or bam file. Should be sorted by query name [None, REQ]')
    p.add_option('-q', '--minQual', default=20, type='int', help='Minimum base quals to be used in a geno [20]')
    p.add_option('-m', '--mapQual', default=20, type='int', help='Minimum mapping quality for a read to be used [20]')
#    p.add_option('-b', '--buffer', default=2, type='int', help='Number of bases to use as a buffer. If two insert locations are within this distance they are considered the same [2]')
    p.add_option('-1', '--first', type='int', help='position in the reference of the left-most base to phase')
    p.add_option('-2', '--second', type='int', help='position in the reference of the right-most base to phase')
    p.add_option('--R1only', default = False, action = "store_true", help='Use this flag to only consider R1')
    p.add_option('--R2only', default = False, action = "store_true", help='Use this flag to only consider R2')

    opts, args = p.parse_args()
    
    if opts.R1only and opts.R2only: print "--R1only and --R2only cannot be used together"
    
    #Batch mode
    elif opts.input:
        fout = open(opts.out, 'w')
        fout.write("File\tPosition 1\tPosition 2\tRef\Ref\tRef\Alt\tAlt\Ref\tAlt\Alt\n")
        for each in args:
            opts.sam = each
            for line in open(opts.input, 'r'):
                cols = line.strip().split("\t")
                hap_tup_dict = phase(int(cols[0]), int(cols[3]), opts)
                exp_geno_tups = make_exp_geno_tups(cols)
                fout.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\n" % (each, cols[0], cols[3], len(set(hap_tup_dict.get(exp_geno_tups[0], []))), len(set(hap_tup_dict.get(exp_geno_tups[1], []))), len(set(hap_tup_dict.get(exp_geno_tups[2], []))), len(set(hap_tup_dict.get(exp_geno_tups[3], [])))))
    
    #Old way of just processing a single pair specified on the command line
    else:
        for each in args:
            opts.sam = each
            hap_tup_dict = phase(opts.first, opts.second, opts)

###-----------------End of main()--------------------------->>>

def make_exp_geno_tups(cols):
    return [(cols[1],cols[4]),(cols[1],cols[5]),(cols[2],cols[4]),(cols[2],cols[5])]

def phase(pos1, pos2, opts):
    hap_tup_dict={}
    mutants={}
    sam = pysam.Samfile(opts.sam)
    #Step through each read in the sam file
    for read in sam:
        #Check to make sure the read is mapped and meets the mapping qulaity threshold
        if not read.is_unmapped and read.mapping_quality >= opts.mapQual:
            #Read type checks
            if (opts.R1only and read.is_read1) or (opts.R2only and read.is_read2) or (not opts.R1only and not opts.R2only):
                #Reference start and end are 0-based, but end points to one past the last base
                if read.reference_start >= read.reference_end: print "!!!!End is NOT larger than start: start=%d, end=%d" % (read.reference_start, read.reference_end)
                if max(pos1-1, pos2-1) < read.reference_end and min(pos1-1, pos2-1) >= read.reference_start:
                        read_ref_dict, ref_read_dict = get_paired_base_dicts(read)
                        first_geno, first_quals = read_base_at_ref_pos(read, read_ref_dict, ref_read_dict, pos1-1)
                        second_geno, second_quals = read_base_at_ref_pos(read, read_ref_dict, ref_read_dict, pos2-1)
                        
                        if min(first_quals) >= opts.minQual and min(second_quals) >= opts.minQual:
                            if (first_geno, second_geno) not in hap_tup_dict: hap_tup_dict[(first_geno, second_geno)] = [read.query_name]
                            else: hap_tup_dict[(first_geno, second_geno)].append(read.query_name)

    hap_tup_dict = rmv_inconsistent_pairs(hap_tup_dict)

    print opts.sam
    for key, val in hap_tup_dict.iteritems():
        print key, len(val), len(set(val))
    
    return hap_tup_dict
###------------->>>

def rmv_inconsistent_pairs(hap_tup_dict):
    for a,b in itertools.combinations(hap_tup_dict.keys(), 2):
        ovlp = set(hap_tup_dict[a]).intersection(set(hap_tup_dict[b]))
        if ovlp:
            for each in ovlp:
                print each, a, b
                hap_tup_dict[a].remove(each)
                hap_tup_dict[b].remove(each)
    return hap_tup_dict


def get_paired_base_dicts(read):
    aligned_base_tuples = read.get_aligned_pairs()
    read_ref_dict=dict(aligned_base_tuples)
#    print read.get_aligned_pairs(with_seq=True)
    ref_read_dict=dict((y, x) for x, y in aligned_base_tuples)
    return read_ref_dict, ref_read_dict

def read_base_at_ref_pos(read, read_ref_dict, ref_read_dict, ref_pos):
    quals=[]
    read_pos = ref_read_dict[ref_pos]
    if not read_pos: return '-', [1000]
    else:
        geno=read.query_sequence[read_pos]
        quals.append(read.query_qualities[read_pos])
        n=1
        while not read_ref_dict[read_pos-n]:
            geno = read.query_sequence[read_pos-n] + geno
            quals.append(read.query_qualities[read_pos])
            n+=1
        return geno, quals
    

if __name__ == "__main__":
    main()
