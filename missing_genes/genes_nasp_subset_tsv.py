#!/usr/bin/env python

from __future__ import division
import sys, optparse, os
from subprocess import Popen, PIPE

def main():
    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)

    #Input/output files
#    p.add_option('-r', '--ref', help='Fasta file that was used as reference in mapping. [None, Required]')
    p.add_option('-g', '--gbk',  help='GenBank file for reference [None, REQ]')
    p.add_option('-o', '--out', help='Name for output file. [None, REQ]')
    p.add_option('-t', '--tsv', help='Subset NASP .tsv file containing only positions of interest. [None]')
    p.add_option('--propCov', type='float', default=0.4, help='If a gene has this proportion, or larger, of bases in the subset tsv, it is considered "covered" [0.4]')
    p.add_option('--recursThresh', type='int', default=50000, help='To reset the recursion threshold. [5000]')
    
    opts, args = p.parse_args()

    #Print out the command that was provided
    print '  '.join(sys.argv)

    #Increase recursion limit before starting script
    sys.setrecursionlimit(opts.recursThresh)

    #This is probably not necessary in the current version, but shouldn't hurt
    #Change all filenames to their absolute path versions
#    if opts.ref: opts.ref=os.path.abspath(opts.ref)

    if opts.gbk and opts.tsv:
        tsv_sites = parse_tsv(opts.tsv)
        gene_info_dict = parse_genbank(opts.gbk)
        covered_genes = []
        for gene, info in gene_info_dict.iteritems():
            contig = info[0]
            if contig in tsv_sites:
                ovrlp = tsv_sites[contig].intersection(info[2])
                if len(ovrlp)/len(info[2])>=opts.propCov:
                    covered_genes.append((info[0], gene, info[1], len(ovrlp)/len(info[2])))
        fout = open(opts.out, "w")
        for each in sorted(covered_genes):
            fout.write("%s\t%s\t%s\t%.4f\n" % (each[0], each[1], each[2], each[3]))
            
    else: print  'Could not determine covered genes, you must supply gene info and a subset tsv file' 
    
# End of main() --------------------------------------------------------------------------------

def parse_tsv(tsv):
    linecount=0
    tsv_dict = {}
    for line in open(tsv, "r"):
        linecount+=1
        if linecount>1:
            cols = line.strip().split()
            contig, pos = cols[0].split("::")
            
            #!!!!Likely a special case scenario
            contig = contig.split('|')[1]
            
            if contig not in tsv_dict: tsv_dict[contig] = set([int(pos)])
            else: tsv_dict[contig].add(int(pos))
    return tsv_dict
            
def parse_genbank(gbk):
    contig = ""
    atfeats = 0
    count=0
    info={x:'' for x in ['ltag', 'prod', 'pos']}
    info_dict={}
    ltag_split = 0
    product_split = 0
    inbetwn = 1
    for line in open(gbk,"r"):
        if line.startswith("LOCUS"):
            atfeats = 0
            contig = line.strip().split()[1]
        if atfeats:
            if ltag_split:
                if line.strip().endswith('"'):
                    info['ltag'] = info['ltag'] + line.strip()[:-1]
                    ltag_split = 0
                else:
                    info['ltag'] = info['ltag'] + line.strip()
            
            elif product_split:
                if line.strip().endswith('"'):
                    info['prod'] = info['prod'] + line.strip()[:-1]
                    product_split = 0
                else:
                    info['prod'] = info['prod'] + line.strip()
            elif line.startswith("     source"):
                length = line.strip().split()[1].split("..")[1]
                if contig.endswith(length): contig = contig[:-len(length)]
            elif line.startswith("     gene"):
                inbetwn = 1
            elif line.startswith("     CDS"):
                inbetwn = 0
                count+=1
                #Save info from previous CDS
                if count>1:
                    info_dict[info['ltag']] = [contig, info['prod'], info['pos']]
                    info={x:'' for x in ['ltag', 'prod', 'pos']}
                coords = line.strip().split()[1]
                if coords.startswith("complement"): coords = coords[11:-1]
                if coords.startswith("join"):
                    info['pos'] = set([])
                    splitcoords = coords[5:-1].split(',')
                    for each in splitcoords:
                        thisstart, thisstop = [int(x.strip("<").strip(">")) for x in each.split("..")]
                        info['pos'].update(set(range(thisstart, thisstop+1)))
                else:
                    start, end = [int(x.strip("<").strip(">")) for x in coords.split("..")]
                    info['pos'] = set(range(start, end+1))
            elif not inbetwn:
                if line.strip().startswith("/locus_tag"):
                    if line.strip().endswith('"'): info['ltag'] = line.strip().split('=')[1][1:-1]
                    else:
                        info['ltag'] = line.strip().split('=')[1][1:]
                        ltag_split = 1
                elif line.strip().startswith("/product"):  
                    if line.strip().endswith('"'): 
#                        print line.strip().split('=')[1][1:-1]
                        info['prod'] = line.strip().split('=')[1][1:-1]
                    else:
                        info['prod'] = line.strip().split('=')[1][1:]
                        product_split = 1
        else:
            if line.startswith("FEATURES"):
                atfeats=1
    return info_dict

def write_missing_repeat_genes(pos_count_dict_dict_std, pos_count_dict_dict_qual, avg_cov, opts):
    fout = open('%s_missing_genes.txt' % opts.out, 'w')
    fout_all = open('%s_all_genes.txt' % opts.out, 'w')
    fout_rep = open('%s_repeat_in_ref_genes.txt' % opts.out, 'w')
    for outfile in [fout, fout_rep]:
        outfile.write('Org_ID\tchrom_ID\tGene\tPerc_Zero_Cov\tCov_Perc_Avg\n')
    fout_all.write('Org_ID\tchrom_ID\tGene\tPerc_Zero_Cov\tCov_Perc_Avg\tQual_Perc_Zero_Cov\tQual_Cov_Perc_Avg\n')
    
        
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



def pos_from_tsv(opts):
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

                                                
###----------------->>>

if __name__ == '__main__':
    main()


