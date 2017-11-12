#!/usr/bin/env python

from __future__ import division
import optparse, os, time
from subprocess import Popen, PIPE
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

#This script reads in an aligned fasta file, creates a vcf, identifies TC substitution clusters and then creates a new fasta file with these types of changes reverted to reference

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()

#From fasta to vcf script
    p.add_option('-f', '--fasta',  help='Aligned fasta. [None, REQ]')
    p.add_option('-r', '--ref',  help='Name of sequence in the alignment to use as a reference. [None, REQ]')
    p.add_option('-o', '--out',  help='Base name for output files. [None, OPT]')
    p.add_option('-g', '--gap', default='-', help='Character used for gaps. [-]')


#From ID TtoC clusters script
    p.add_option('-v', '--vcf',  help='vcf file. [None, OPT]')
    p.add_option('-w', '--win',  type='int', default=200, help='Window size to check for enrichment. [200]')
    p.add_option('-n', '--num',  type='int', default=4, help='Min number of changes for enrichment. [4]')

    opts, args = p.parse_args()
    
    #Set opts.out using fasta file name, if not specified on the command line
    if not opts.out:
        opts.out = '.'.join(opts.fasta.split('.')[:-1])
    
    if not opts.vcf:
        opts.vcf = make_vcf(opts)
    
    #Get dict with entry for each seq with a TC cluster
    #Keys are seq names, values are lists of lists of cluster positions from the vcf
    cluster_dict = get_cluster_info(opts)

    #Make a new fasta where these potential ADAR-edited sites have been reverted back to Ts
    reverse_ADAR(cluster_dict, opts)
    
    #Get context (from the reference sequence) for the putative ADAR-edited sites
    get_context(cluster_dict, opts)
#----------------------End of main()

def get_context(cluster_dict, opts):
    #Get ref seq
    fasta_dict = read_fasta_dict(opts.fasta)
    refseq = fasta_dict[opts.ref]
    #Dict to move between vcf pos and alignment pos
    vcf_to_align_pos_dict = connect_vcf_to_align(fasta_dict[opts.ref], opts.gap)

    #Make empty lists to hold 5' and 3' bases
    prime5=[]
    prime3=[]

    #Get list of unique sites with ADAR editing
    uniq_sites = only_uniq(cluster_dict.values())
    for site in uniq_sites:
        prime5.append(refseq[vcf_to_align_pos_dict[site-1]].upper())
        prime3.append(refseq[vcf_to_align_pos_dict[site+1]].upper())

    #Print stats
    print_stats(prime5, '5 prime')
    print_stats(prime3, '3 prime')
        

def print_stats(thelist, name):
    print '%s:' % name
    for each in set(thelist):
        print "\t%s: %d (%.2f)" % (each, thelist.count(each), thelist.count(each)/len(thelist)*100)
    
def only_uniq(clusters):
    master_list=[]
    for each in clusters:
        for single in each: master_list+=single
    return list(set(master_list))

def reverse_ADAR(cluster_dict, opts):
    fasta_dict = read_fasta_dict(opts.fasta)
    vcf_to_align_pos_dict = connect_vcf_to_align(fasta_dict[opts.ref], opts.gap)

    #Step through each seq with a TC cluster and replace the seq in the seqs li
    for seq in cluster_dict.keys():
        fasta_dict[seq] = replace_ADAR_Cs(fasta_dict[seq], cluster_dict[seq], vcf_to_align_pos_dict)

    #Write new fasta
    write_fasta(fasta_dict.keys(), fasta_dict.values(), "%s_ADARreversed.fasta" % (opts.out))


def replace_ADAR_Cs(seq, clusters, pos_dict):
    for clust in clusters:
        for pos in clust:
            if seq[pos_dict[pos]].upper() != 'C': print "PROBLEM!: %d/%d supposed to be C, actually %s" % (pos, pos_dict[pos], seq[pos_dict[pos]])
            seq = seq[:pos_dict[pos]] + 'T' + seq[pos_dict[pos]+1:]
    return seq


def connect_vcf_to_align(ref_seq, gap_char):
    pos_dict={}
    ref_count=0
    for index, base in enumerate(ref_seq):
        if base != gap_char:
            ref_count+=1
            pos_dict[ref_count]=index
    return pos_dict

def get_cluster_info(opts):
    tc_dict={}
    #Step through the vcf file and get info for T-C transtions
    for line in open(opts.vcf, 'r'):
        #If header row
        if line.startswith('#') and not line.startswith('##'):
            cols=line.strip().split('\t')
            ids=cols[9:]
        #If data row
        elif not line.startswith('#'):
            #print 'data'
            cols=line.strip().split('\t')
            #If this is a T to C transition
            #print cols[3], cols[4]
            if cols[3].upper()=='T' and cols[4].upper()=='C':
                #print 'TC'
                c_indices=[i for i,x in enumerate(cols[9:]) if x.upper()=='C']
                for i in c_indices:
                    if ids[i] not in tc_dict: tc_dict[ids[i]]=[int(cols[1])]
                    else: tc_dict[ids[i]].append(int(cols[1]))
    #Step through seqs and look for windows of T to C mutation enrichment
    cluster_dict={}
    
    #Open file for output
    fout = open("%s_clusters.txt" % opts.out, 'w')
    
    #print len(tc_dict)
    for seq, muts in tc_dict.iteritems():
        #print seq, len(muts)
        clusters=find_clusters(muts, opts)
        if clusters:
            fout.write("%s\t%s\n" % (seq, "\t".join([",".join([str(y) for y in sorted(x)]) for x in clusters])))
            cluster_dict[seq]=clusters
    
    fout.close()
    
    return cluster_dict

#-------START ----- Functions from the fasta to vcf script

def make_vcf(opts):

    #Read in seqs
    names, seqs = read_fasta_lists(opts.fasta)
    refseq=seqs[names.index(opts.ref)]
    #If the sequences are unaligned (not all the same length)
    if len(set([len(x) for x in seqs]))>1:
        print 'Sequences must all be the same length!!!'
        return

    else:
        vcf_name = "%s_%s.vcf" % (opts.out, opts.ref)
        fout=open(vcf_name, 'w')
        fout.write('##fileformat=VCFv4.2\n##fileDate=%s\n##source=ID_revert_TC_clusters.py\n##reference=%s\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (time.strftime("%m-%d-%Y"), opts.ref, '\t'.join(names)))
        #Get info for variable positions
        idxs, alts, positions = align_stats(seqs, refseq, opts)
        #Write out info to vcf file 
        for list_idx, var_idx in enumerate(idxs):
            for a in alts[list_idx]:
                fout.write('%s\t%d\t.\t%s\t%s\t.\t.\t.\t.\t%s\n' % (opts.ref, positions[list_idx], refseq[var_idx], a, genos(seqs, var_idx, refseq[var_idx].upper())))
    fout.close()

    return vcf_name

def genos(seqs, var_idx, ref_base):
    these_genos=[]
    bases=[x[var_idx].upper() for x in seqs]
    for b in bases:
        if b==ref_base: these_genos.append('.')
        else: these_genos.append(b)
    return '\t'.join(these_genos)

def seq_stats_info(s):
    var_i=[]
    for i in range(len(s)):
        if s[i].upper()=='N': var_i.append(i)
    return var_i
        
def perc(num, den):
    return num/den*100

def align_stats(seqs, ref, opts):
    positions=[]
    idxs=[]
    alts=[]
    ref_count=0
    for i in range(len(ref)):
        if ref[i] != opts.gap: ref_count+=1
#        else: print 'Gap at %d' % (i+1)
        uniq_bases=set([x[i].upper() for x in seqs]).difference(set(['N']))
        if len(uniq_bases) > 1: 
            idxs.append(i)
            alts.append(list(uniq_bases.difference(set([ref[i].upper()]))))
            positions.append(ref_count)
    return idxs, alts, positions
    
# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()

#-------END ----- Functions from the fasta to vcf script



#-------START ----- Functions from the ID TtoC clusters script

def find_clusters(mut_list, opts):
    clusters=[]
    mut_list.sort()
    i=0
    while i<=len(mut_list)-opts.num:
        y=opts.num-1
        if (mut_list[i+y]-mut_list[i])<opts.win:
            while (mut_list[i+y]-mut_list[i])<opts.win:
                y+=1
                if i+y >= len(mut_list): break
            clusters.append(mut_list[i:i+y])
        i+=1
    return condense(clusters)

def condense(clusters):
    out=[]
    while len(clusters)>0:
        first, rest = clusters[0], clusters[1:]
        first = set(first)
        
        lenf = -1
        while len(first)>lenf:
            lenf = len(first)
            
            rest2=[]
            for r in rest:
                if len(first.intersection(set(r)))>0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest=rest2
        out.append(list(first))
        clusters = rest
    return out

#### Not currently using the functions below here

def parse_variant(var_str):
    if var_str.startswith("ANN="):
        var_str=var_str.split("ANN=")[1]
        var_types=[x.split('|')[1] for x in var_str.split(',')]
        var_locs=[x.split('|')[3] for x in var_str.split(',')]
        if len(var_types)==1: return var_types[0], var_locs[0]
        else:
            dom_index=get_dom_index(var_types)
            return var_types[dom_index], var_locs[dom_index]
    else: return "intergenic_region", "terminus"

def get_dom_index(var_list):
    if 'stop_gained' in var_list: return var_list.index('stop_gained')
    if 'frameshift_variant' in var_list: return var_list.index('frameshift_variant')
    elif 'missense_variant' in var_list: return var_list.index('missense_variant')
    elif 'synonymous_variant' in var_list: return var_list.index('synonymous_variant')
    else: print var_list
    
def read_fasta_dict(file):
    names, seqs = read_fasta_lists(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

#-------END ----- Functions from the ID TtoC clusters script


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

