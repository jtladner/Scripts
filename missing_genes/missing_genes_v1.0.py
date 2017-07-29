#!/usr/bin/env python

from __future__ import division
import sys, optparse, os
from subprocess import Popen, PIPE

def main():
    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)

    #Input/output files
    p.add_option('-u', '--unpaired', help='Fastq files with unpaired reads to be mapped to reference. For this purpose, all files should be treated as unpaired. You can combine the mapping of multiple fastqs by entering them together as a comma separated list [None]')
    p.add_option('-r', '--ref', help='Fasta file that was used as reference in mapping. [None, Required]')
    p.add_option('-o', '--out', default='unknownID', help='Base name for output files. Output files will be written to the current working directory. [unknownID]')

#    p.add_option('-l', '--min_length', type='int', default=10, help='Minimum region length to plot [10]')
#    p.add_option('-g', '--genes', help='Gene info for the reference strain')
    p.add_option('--mapQ', type='int', default=20, help='Minimum mapping quality for a read to be used in the pileup generation for the non-repeat version [20]')
    p.add_option('--maxCov', type='int', default=500, help='Max per base coverage to be used in the pileup generation for the final quality check [500]')
    p.add_option('--recursThresh', type='int', default=50000, help='Minimum level of coverage required to change consensus. [3000]')
    p.add_option('--offset', type='int', default=33, help='Base quality offset used in the pileup. I believe the default in samtools is Sanger(33) [33]')
    p.add_option('--baseQual', type='int', default=20, help='Minimum base quality for a base to be counted when looking at coverage. [20]')
    p.add_option('--procs', type='int', default=16, help='Number of processors to use in multi-threaded portions [16]')
    #p.add_option('--minOvl', type='int', default=20, help='Minimum overlap between missing regions in normal and mapq pileup for the region to be reported [20]')
    p.add_option('--scoreMin', default='L,0,-0.12', help='Minimum score for a good alignment in bowtie2. For 100bp reads, L,0,0=0mismatches, L,0,-0.06=1, L,0,-0.12=2, L,0,-0.18=3, L,0,-0.24=4, L,0,-0.30=5, L,0,-0.6=10. [L,0,-0.12]')

    #Optional starting places
    p.add_option('-b', '--bam', help='Bam file from previous run. **Optional starting place. [None]')
    p.add_option('-s', '--sam', help='Sam file from previous run. **Optional starting place. [None]')
    p.add_option('-p', '--std_pile', help='Basic pileup from previous run. **Optional starting place. [None]')
    p.add_option('--qual_pile', help='Pileup from previous run that only contains high quality mapped reads. **Optional starting place. [None]')
    p.add_option('--miss_regions', help='Missing regions file from previous run. **Optional starting place. This automatically throws the --justGenes flag. [None]')
    p.add_option('-g', '--genes',  help='Gene info for the reference strain [None, REQD]')
    p.add_option('--justGenes', action='store_true', default=False, help='Use this flag if you already have the missing regions file and you just want to find the genes that correspond')
    
    #To determine if a gene is missing
    p.add_option('--maxZeroPerc', type='float', default=0.4, help='If a gene has a larger percent of bases with zero coverage than this, it is considered missing [0.4]')
    p.add_option('--minAvgCov', type='float', default=0.01, help='If a gene has a lower percent of average coverage than this, it is considered missing [0.01]')
    
    opts, args = p.parse_args()

    #Print out the command that was provided
    print '  '.join(sys.argv)

    #If both pileups are provided, automatically set justGenes = True
    if opts.std_pile and opts.qual_pile: opts.justGenes=True

    #Increase recursion limit before starting script
    sys.setrecursionlimit(opts.recursThresh)

    #This is probably not necessary in the current version, but shouldn't hurt
    #Change all filenames to their absolute path versions
    if opts.ref: opts.ref=os.path.abspath(opts.ref)
    if opts.unpaired: opts.unpaired=comma_sep_abs_path(opts.unpaired)

    #Create pileups, based on whatever starting material you have
    if not opts.justGenes:
        if not (opts.std_pile and opts.qual_pile):
            if not opts.bam:
                if not opts.sam:
                    if not opts.unpaired: 
                        print 'Insufficient starting material for calculation of missing regions, must provide either fastq (-u), sam (-s), bam (-b) or pileups (-p)!'
                        return
                    opts.sam = bowtie2_index_align(opts)
                opts.bam = make_sort_bam(opts)
            opts.std_pile, opts.qual_pile = make_pileups(opts)


    #If gene info is provided, then will write a file containing all of the genes that are at least partially missing
    if opts.genes and opts.std_pile and opts.qual_pile:
        pos_count_dict_dict_std, avg_cov_std = get_pos_counts(opts.std_pile, opts)
        pos_count_dict_dict_qual, avg_cov_qual = get_pos_counts(opts.qual_pile, opts)
        write_missing_repeat_genes(pos_count_dict_dict_std, pos_count_dict_dict_qual, avg_cov_qual, opts)
    else: print  'Could not write missing genes, you did not supply gene info and/or enough information to make the pileup' 
    
# End of main() --------------------------------------------------------------------------------

def write_missing_repeat_genes(pos_count_dict_dict_std, pos_count_dict_dict_qual, avg_cov, opts):
    gene_info_dict = make_gid(opts.genes)
    fout = open('%s_missing_genes.txt' % opts.out, 'w')
    fout_all = open('%s_all_genes.txt' % opts.out, 'w')
    fout_rep = open('%s_repeat_in_ref_genes.txt' % opts.out, 'w')
    for outfile in [fout, fout_rep]:
        outfile.write('Org_ID\tchrom_ID\tGene\tPerc_Zero_Cov\tCov_Perc_Avg\n')
    fout_all.write('Org_ID\tchrom_ID\tGene\tPerc_Zero_Cov\tCov_Perc_Avg\tQual_Perc_Zero_Cov\tQual_Cov_Perc_Avg\n')
    
    for gene, info in gene_info_dict.iteritems():
        num_zero_bases=0
        num_zero_bases_qual=0
        base_covs=[]
        base_covs_qual=[]
        
        gene_bases=range(int(info[2]), int(info[3])+1)
        
        #Get info from standard pileup
        for chrom in pos_count_dict_dict_std:
            if info[0] in chrom:
                for pos in gene_bases:
                    print chrom, len(pos_count_dict_dict_std[chrom]), pos
                    this_cov=pos_count_dict_dict_std[chrom][pos]
                    base_covs.append(this_cov)
                    if this_cov == 0: num_zero_bases+=1
        
        #Get info from quality controlled pileup            
        for chrom in pos_count_dict_dict_qual:
            if info[0] in chrom:
                for pos in gene_bases:
                    this_cov=pos_count_dict_dict_qual[chrom][pos]
                    base_covs_qual.append(this_cov)
                    if this_cov == 0: num_zero_bases_qual+=1
                    
        #Just a check to make sure everything is working right
        if len(base_covs_qual) == len(base_covs) == len(gene_bases):
            perc_zero = num_zero_bases/len(gene_bases)
            perc_avg_cov = (sum(base_covs)/len(base_covs))/avg_cov
            
            perc_zero_qual = num_zero_bases_qual/len(gene_bases)
            perc_avg_cov_qual = (sum(base_covs_qual)/len(base_covs_qual))/avg_cov
            
            fout_all.write('%s\t%s\t%s\t%s\t%.2f\t%s\t%.2f\n' % (gene, info[0], info[1], perc_zero, perc_avg_cov, perc_zero_qual, perc_avg_cov_qual))
            #If it seems to be missing in the std pileup
            if (perc_zero >= opts.maxZeroPerc or perc_avg_cov < opts.minAvgCov):
                fout.write('%s\t%s\t%s\t%s\t%.2f\n' % (gene, info[0], info[1], perc_zero, perc_avg_cov))
            #If it only seems to be missing in the quality controlled pileup
            elif (perc_zero_qual >= opts.maxZeroPerc or perc_avg_cov_qual < opts.minAvgCov):
                fout_rep.write('%s\t%s\t%s\t%s\t%.2f\n' % (gene, info[0], info[1], perc_zero_qual, perc_avg_cov_qual))
        
        else: print 'There is a problem with %s, either returning info for multiple chromosomes or missing some bases' % gene
        
    fout.close()
    fout_all.close()
    fout_rep.close()

#Keys = Org's gene ID
#Values are lists with [chrom, gene_name, start, end]
def make_gid(gene_file):
    fin = open(gene_file, 'r')
    gid={}
    linecount=0
    for line in fin:
        linecount+=1
        if linecount>1:
            cols = line.strip().split('\t')
            gid[cols[0]] = cols[1:]
    fin.close()
    return gid


##!!! Expecting fastq quality scores to be 33-offset
def bowtie2_index_align(opts):

    #Make index for reference
    cmd='bowtie2-build %s ref_index' % (opts.ref)
    os.popen(cmd)

    sam_name='%s_aligned.sam' % opts.out

    cmd='bowtie2 -p %d --ignore-quals --phred33 -N 0 --score-min %s -x ref_index -U %s -S %s' % (opts.procs, opts.scoreMin, opts.unpaired, sam_name)
    print cmd
    unpaired_bowtie=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    unpaired_bowtie.wait()

    return sam_name

def make_sort_bam(opts):
    #Make bam from sam, and sort the bam
    sam_prefix='.'.join(opts.sam.split('.')[:-1])
    bam_name='%s_sorted.bam' % sam_prefix
    cmd='samtools view -bS %s | samtools sort -m 800000000 - %s' % (opts.sam, sam_prefix + '_sorted')
    make_bam=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_bam.wait()
    return bam_name    

def make_pileups(opts):
    #Make new faidx for current reference
    cmd='samtools faidx %s' % opts.ref
    faidx=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    faidx.wait()

    #Make pileup from bam
    pile_name='%s_mapQ0.pile' % ('.'.join(opts.bam.split('.')[:-1]))
    cmd='samtools mpileup -f %s -q 0 -d %d  %s >%s' % (opts.ref, opts.maxCov, opts.bam, pile_name)
    std_pileup=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    std_pileup.wait()

    #Make pileup from bam
    qual_pile_name='%s_mapQ%d.pile' % ('.'.join(opts.bam.split('.')[:-1]), opts.mapQ)
    cmd='samtools mpileup -f %s -q %d -d %d  %s >%s' % (opts.ref, opts.mapQ, opts.maxCov, opts.bam, qual_pile_name)
    mapq_pileup=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    mapq_pileup.wait()

    return pile_name, qual_pile_name


def make_fasta_pos_dict_simple_names(ref):
    chrom_pos_dict_dict={}
    names, seqs = read_fasta_lists_simple(ref)
    for index in range(len(names)):
        chrom_pos_dict_dict[names[index]] = {}
        for base in range(1,len(seqs[index])+1): 
            chrom_pos_dict_dict[names[index]][base]=0
    return chrom_pos_dict_dict


def get_pos_counts(pileup, opts):
    #Make a dictionary from the reference fasta
    ref_pos_dict_dict=make_fasta_pos_dict_simple_names(opts.ref)
    non_zero_cov_counts=[]
    
    fin=open(pileup, 'r')
    for line in fin:
        chrom, pos, bases, indels, ends, quals = pileup_info(line)
        #Removes info for bases with quality lower than the specified minimum threshold
        good_bases, good_quals = quality_filter(bases, quals, opts)
        ref_pos_dict_dict[chrom][pos]=len(good_bases)
        non_zero_cov_counts.append(len(good_bases))
    #Just in case a whole contig doesn't have any coverage
#    if not non_zero_cov_counts: non_zero_cov_counts=[1]
    return ref_pos_dict_dict, sum(non_zero_cov_counts)/len(non_zero_cov_counts)


def quality_filter(bases, quals, opts):
    new_bases=[]
    new_quals=[]
    for index in range(len(bases)):
        if ord(quals[index])-opts.offset >= opts.baseQual:
            new_bases.append(bases[index])
            new_quals.append(quals[index])
    return new_bases, new_quals


def pileup_info(line):
    cols=line.strip().split('\t')
    chrom=cols[0]
    pos=int(cols[1])
    num_reads=int(cols[3])
    bases, indels, ends = parse_bases(cols[4].upper(), num_reads)
    quals=cols[5]
    return chrom, pos, bases, indels, ends, quals


def parse_bases(raw_bases, num_reads):
    #If any reads begin at this position, remove these markers and the associated mapping quality info
    if '^' in raw_bases:
        bases=remove_beg_info(raw_bases)
    else: bases=raw_bases
    #Count the number of reads ending and then remove the $ charcters
    ends=bases.count('$')
    bases=bases.replace('$', '')

    #Get info about indels following this position and remove this info form bases
    indels={}           #This line is needed to ensure an empty dict named indels for strings without indel present
    if '-' in bases or '+' in bases:
        bases, indels = indel_info(bases, indels)

    if len(bases)==num_reads: return bases, indels, ends
    else:
        print 'Problem!!! Exp bases: %d, Output bases: %d, Input str: %s, Output str: %s' % (num_reads, len(bases), raw_bases, bases)
        return #This will cause script to stop because three arguments are expected to be returned

def remove_beg_info(bases):
    while '^' in bases:
        index=bases.index('^')
        bases=bases[:index]+bases[index+2:]
    return bases

def indel_info(bases, ind_dict):
    ind_dict={}
    if '-' in bases:
        del_info=[]
        new_bases, del_info = characterize_indel(bases, '-', del_info)
        ind_dict['delete']=del_info
        bases=new_bases
    if '+' in bases:
        ins_info=[]
        new_bases, ins_info = characterize_indel(bases, '+', ins_info)
        ind_dict['insert']=ins_info
        bases=new_bases
    return bases, ind_dict

def characterize_indel(bases, type_char, ind_info):
    start_index=bases.index(type_char)
    digitcount=0
    for each in bases[start_index+1:]:
        if each.isdigit(): digitcount+=1
        else: break
    indel_bases=int(bases[start_index+1:start_index+1+digitcount])
    indel_str=bases[start_index+1+digitcount:start_index+1+digitcount+indel_bases]

    ind_info.append(indel_str)
    bases=bases[:start_index]+bases[start_index+1+digitcount+indel_bases:]
    #Recursively go through this process if the string has more indel info
    if type_char in bases: bases, ind_info = characterize_indel(bases, type_char, ind_info)

    return bases, ind_info



def comma_sep_abs_path(string):
    file_list=string.split(',')
    abs_paths=[os.path.abspath(x) for x in file_list]
    return ','.join(abs_paths)

def get_all_subkeys(dd):
    all_keys=[]
    for dict in dd.values():
        all_keys+=dict.keys()
    return list(set(all_keys))


def make_gen_info_file(reference):
    new_filename = '%s_infofile_for_coverage.txt' % "".join(reference.split('.')[:-1])
    chrom_lengths={}
    fout = open(new_filename, 'w')
    names, seqs = read_fasta_lists(reference)
    for index in range(len(names)):
        fout.write('%s\t%s\n' % (names[index], len(seqs[index])))
        chrom_lengths[names[index].split()[0]] = len(seqs[index])
    fout.close()
    return new_filename, chrom_lengths

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (
def read_fasta_lists_simple(file):
        fin = open(file, 'r')
    count=0
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:].split()[0])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
                                                                                                                                                                                                                                        
    return names, seqs

def sam2bam(samfile):
    sam_prefix=samfile.split('.')[0]
    cmd='samtools view -bS %s | samtools sort -m 800000000 - %s' % (samfile, sam_prefix + '_sorted')
    print cmd
    os.popen(cmd)
    return '%s_sorted.bam' % sam_prefix
                                                                                                                                                                                                                                                
                                                
###----------------->>>

if __name__ == '__main__':
    main()


