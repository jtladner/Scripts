#!/usr/bin/env python

# By Jason Ladner

###!!! Requires installations of samtools, bedtools(genomeCoverageBed) and bowtie2, both of which must reside within your $PATH

from __future__ import division
import sys, optparse, os, datetime
from subprocess import Popen, PIPE
#For plotting
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties as fp
fontP = fp()

def main():

    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
    p.add_option('-u', '--unpaired', help='Fastq file with unpaired reads to be mapped to reference [None]')
    p.add_option('-1', '--paired1', help='Fastq file with read 1 of paired reads [None]')
    p.add_option('-2', '--paired2', help='Fastq file with read 2 of paired reads [None]')
    p.add_option('-r', '--ref', help='Fasta file that was used as reference in mapping. [None, Required]')
    p.add_option('-o', '--out', default='new_consensus', help='Base name for output files. Output files will be written to the current working directory. [new_consensus]')
    
    #General
    p.add_option('--recursThresh', type='int', default=50000, help='Minimum level of coverage required to change consensus. [3000]')
    p.add_option('--temp', default='./temp', help='Name for temporary working directory. This will be created at the begining of the script and then removed at the end [./temp]')
    p.add_option('--procs', type='int', default=8, help='Number of processors to use in multi-threaded portions [8]')
    p.add_option('--Nbeg', type='int', default=10000, help='# of Ns to add to the beginning of each contig [10000]')
    p.add_option('--Nend', type='int', default=10000, help='# of Ns to add to the end of each contig [10000]')    
    
    #For samtools
    p.add_option('--mapQ', type='int', default=10, help='Minimum mapping quality for a read to be used in the pileup generation [10]')
    p.add_option('--maxCov', type='int', default=300, help='Max per base coverage to be used in the pileup generation [300]')
    #This step allows for a large amount of time saved when a small number of the reads actually map to the reference. Should probably be used for viral samples, but not bacterial.
    p.add_option('--NOonlyMapped', default=False, action='store_true', help='Use this flag to turn off the part of the script that creates a smaller sam with only mapped reads prior to sorting')
    
    #For bowtie2
    p.add_option('--nceil', default='0,0.6',  help='Function to supply as --n-ceil parameter in bowtie2. [0,0.6]')
    p.add_option('--minIns', type='int', default=0, help='Minimum insert size for a pair to be concordant [0]')
    p.add_option('--maxIns', type='int', default=500, help='Maximum insert size for a pair to be concordant [500]')
    
    
    #Make new consensus
    p.add_option('--propThresh', type='float', default=0.5, help='Requires that a greater proportion of the reads (than this threshold) must support an alternative for the consensus to be changed. [0.5]')
    p.add_option('--covThresh', type='int', default=3, help='Minimum level of coverage required to change consensus. [3]')
    p.add_option('--offset', type='int', default=33, help='Base quality offset used in the pileup. I believe the default in samtools is Sanger(33) [33]')
    p.add_option('--baseQual', type='int', default=20, help='Minimum base quality for a base to be used in the evalution of the consensus. [20]')
    
    #Specifies whether to use CAP3 to try to join contigs after best consensus has been built
    #If CAP3 is run, bowtie_buildout will be rerun on the new contigs, if any changes are made
    p.add_option('--NOcap3Join', default=False, action='store_true', help='Use this flag to turn off the CAP3 contig joining step')
    
    #For final quality check
    p.add_option('--NoReqSameBase', default=False, action='store_true', help='Use this flag to turn off the default behavior, which requires at least --finalCov coverage of a single base to be in the final assembly')
    p.add_option('--finalCov', type='int', default=5, help='Minimum level of coverage required to keep a base in the final assembly. If a base is below this threshold, it is changed to an "N" [5]')
    p.add_option('--finalMapQ', type='int', default=20, help='Minimum mapping quality for a read to be used in the pileup generation for the final quality check [20]')
    p.add_option('--finalMaxCov', type='int', default=500, help='Max per base coverage to be used in the pileup generation for the final quality check [500]')
    p.add_option('--Nfinal', type='int', default=15, help='# of Ns to add to the beginning and end of each contig prior to the final mapping. This allows for easier mapping at the ends of the contigs.[15]')
    
    opts, args = p.parse_args()
    
    #Increase recursion limit before starting script
    sys.setrecursionlimit(opts.recursThresh)
    
    #Open a file to write summary info to
    opts.fout_sum = open('%s_summary.txt' % opts.out, 'w')
    #Write current date and time to this summary file
    opts.fout_sum.write('Date Run:\n\t%s\n' % str(datetime.datetime.now()))
    #Write command run to this summary file
    opts.fout_sum.write('Command_run (from within %s):\n\t%s\n\n' % (os.getcwd(), '  '.join(sys.argv)))
    
    #Change all filenames to their absolute path versions
    opts.ref=os.path.abspath(opts.ref)
    if opts.unpaired: opts.unpaired=comma_sep_abs_path(opts.unpaired)
    if opts.paired1: opts.paired1=comma_sep_abs_path(opts.paired1)
    if opts.paired2: opts.paired2=comma_sep_abs_path(opts.paired2)
    
    #Save current working directory
    opts.startDir=os.getcwd()

    #Run the bowtie buildout process
    new_cons_fasta=bowtie_build(opts)
    
    #Join contigs with CAP3, unless NOcap3Join flag is used
    if not opts.NOcap3Join:
        cmd='cap3 %s' % (new_cons_fasta)
        cap3_join=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        cap3_join.wait()
        
        #Clean extra files
        for ext in ['.cap.contigs.links', '.cap.contigs.qual', '.cap.info', '.cap.ace']: 
            if os.path.isfile('%s%s' % (new_cons_fasta, ext)):
                os.remove('%s%s' % (new_cons_fasta, ext))
        
        #If CAP3 joined some of the contigs, then create new consensus and rerun bowtie buildout
        if cap3_improved(new_cons_fasta):
            #Join CAP3 'conitgs' and 'singlets' into a new fasta file, and save this as the current reference
            opts.ref=join_cap3_contigs(new_cons_fasta)
            opts.out='cap3_%s' % opts.out
            #Rerun bowtie buildout on the new consensus 
            final_cons_fasta=bowtie_build(opts)
            #Summarize both changlog files
            summarize_changelog(opts, '%s.changelog' % opts.out, '%s.changelog' % opts.out[5:])
        
        #Otherwise the pre-Cap3 contigs are the final contigs
        #And there is only one changelog file to summarize
        else: 
            final_cons_fasta=new_cons_fasta
            summarize_changelog(opts, '%s.changelog' % opts.out)
        
        #Delete CAP3 contig and singlet files
        for ext in ['.cap.contigs', '.cap.singlets']: 
            if os.path.isfile('%s%s' % (new_cons_fasta, ext)):
                os.remove('%s%s' % (new_cons_fasta, ext))

    #If using NOcap3Join, sets the fasta from the bowtie buildout as the final fasta
    else: 
        final_cons_fasta=new_cons_fasta
        summarize_changelog(opts, '%s.changelog' % opts.out)

    #Finally, map the reads onto the new reference, create coverage plots and note any regions with zero coverage or not enough good coverage to be verified
    coverage_check(final_cons_fasta, opts)

###-----------------End of main()--------------------------->>>

#Look through changelog and summarize results in the summary file
def summarize_changelog(opts, *args):
    count_dict_dict={}
    contig_list=[]
    
    for cl in args:
        fin=open(cl, 'r')
        for line in fin:
            #Extract info from this line
            contig, pos, ref, change, perc, cov=line.split('\t')
            #Add empty count dict if contig is not already represented
            if contig not in count_dict_dict: 
                count_dict_dict[contig]={'num_N':0, 'num_other':0, 'num_ins':0, 'num_del': 0}
                contig_list.append(contig)
            #Add the appropriate count for this line
            if ref=='-': count_dict_dict[contig]['num_ins']+=1
            elif change=='-': count_dict_dict[contig]['num_del']+=1
            elif ref=='N': count_dict_dict[contig]['num_N']+=1
            else: count_dict_dict[contig]['num_other']+=1
    #Write results
    opts.fout_sum.write('\nSummary of changes made to seqs in buildout process:\n')
    for con in contig_list:
        opts.fout_sum.write('\t%s\n\t\tInsertions: %d\n\t\tDeletions: %d\n\t\tSNP change from N: %d\n\t\tSNP change from Non-N: %d\n' % (con, count_dict_dict[con]['num_ins'], count_dict_dict[con]['num_del'], count_dict_dict[con]['num_N'], count_dict_dict[con]['num_other']))
         
            
def comma_sep_abs_path(string):
    file_list=string.split(',')
    abs_paths=[os.path.abspath(x) for x in file_list]
    return ','.join(abs_paths)

###*******************Adds coverage functionality***************

def coverage_check(contig_fasta, opts):
    #Add Ns to the beginning and ends of contigs if specified. Ns_ref variable makes it easy to remove this temporary file at the end
    Ns_ref=add_Ns(contig_fasta, opts.Nfinal, opts.Nfinal)

    #remaps reads to final reference
    sam = bowtie2_index_align(Ns_ref, opts)
    #converts sam to bam and makes pileup
    if opts.NOonlyMapped: bam, pile = make_sort_bam_pileup(sam, Ns_ref, opts.finalMapQ, opts.finalMaxCov)
    else:
        mapped_sam = sam_just_mapped(sam)
        bam, pile = make_sort_bam_pileup(mapped_sam, Ns_ref, opts.finalMapQ, opts.finalMaxCov)
    
    #Quality coverage info
    updated_contigs = check_quality_coverage(pile, Ns_ref, opts) 
    
    #If you have paired data, run test for chimeric regions
    if opts.paired1 and opts.paired2:
        pair_cov_test(sam, opts)
    else: opts.fout_sum.write('\nNo paired data, cannot run test for chimeric contigs\n')
    
    #Absolute Coverage info
    cov_files, chrom_length_dict, genome_info_file = calc_covs([bam], Ns_ref)
    missing_dict=find_missing_regions(cov_files, chrom_length_dict)
    missing_region_files = write_missing_regions(missing_dict)
    plot_names=plot_coverage(cov_files)
    
    #Clean up junk, but keep bam
    for file in [genome_info_file, Ns_ref, sam, mapped_sam, 'ref_index.1.bt2', 'ref_index.2.bt2', 'ref_index.3.bt2', 'ref_index.4.bt2', 'ref_index.rev.1.bt2', 'ref_index.rev.2.bt2']: 
        if os.path.isfile(file):
            os.remove(file)
    
    #Rename some generically named files
    for file in [pile, bam]:
        if os.path.isfile(file):
            os.rename(file, '%s_%s' % (opts.out, file.split('/')[-1]))

###---------Function to make file noting the bases that do not have enough high quality coverage to be verified with the current settings

##!!!I need to change this so that it looks for the number of reference bases. 
##!!!Also, it should change any bad bases to 'N'
##!!!This is probably also the place to look for potentially ambiguous bases
def check_quality_coverage(pileup, ref, opts):
    #Make a dictionary from the reference fasta
    ref_dict=read_fasta_dict_simple_names(ref)
    #Turn each seq into a list 
    for name, seq in ref_dict.items():
        ref_dict[name]=list(seq)

    #Make two dicts with the same keys as ref_dict, but with values that are empty lists
    #These lists will hold all of the base positions that are present for a given contig in the pileup (pos_dict) and all of the positions with insufficient coverage (weak_dict)
    #The bases that are not present, have no coverage of good enough quality
    pos_dict={}
    weak_dict={}
    strong_dict={}
    #Will be used to keep track of the number of bases that lack the proper level of quality coverage, per contig
    num_bases_lack_qual_cov={}
    for name in ref_dict.keys():
        pos_dict[name]=[]
        weak_dict[name]=[]
        strong_dict[name]={opts.finalCov:0, opts.finalCov*2:0, opts.finalCov*5:0, opts.finalCov*10:0,50:0, 100:0, 500:0, 1000:0}
        num_bases_lack_qual_cov[name]=0
    
    #Open file to write out info for bases that do not have enough good quality coverage to be verified
    qual_cov_out='%s_bases_lack_qual_cov.txt' % opts.out
    fout=open(qual_cov_out, 'w')
    fout.write('Chromosome\tPosition\tTotalBases\tTotalRefBases\tGoodBases\tGoodRefBases\n')
    
    
    #Opens pileup file and starts stepping through each line
    fin=open(pileup, 'r')
    for line in fin:
        #Get info about this position
        chrom, pos, bases, indels, ends, quals = pileup_info(line)
        
        #Add position to pos_dict
        pos_dict[chrom].append(pos)
        
        #Extract only 'good' quality bases
        good_bases, good_quals = quality_filter(bases, quals, opts)
        
        #Get reference base counts
        ref_base_count=bases.count('.') + bases.count(',')
        good_ref_base_count = good_bases.count('.') + good_bases.count(',')
        
        #Add counts to strong_dict is coverage is above a certain level
        for num in set([opts.finalCov, opts.finalCov*2, opts.finalCov*5, opts.finalCov*10, 50, 100, 500, 1000]):
            if len(good_bases) >= num:
                strong_dict[chrom][num]+=1 
        
        #The way this line is currently written, the specified amount of coverage must be met when just looking at bases that match the reference
        if good_ref_base_count<opts.finalCov:
            fout.write('%s\t%s\t%d\t%d\t%d\t%d\n' % (chrom, pos, len(bases), ref_base_count, len(good_bases), good_ref_base_count))
            num_bases_lack_qual_cov[chrom]+=1
            weak_dict[chrom].append(pos)
            
    fout.close()
    
    #Change insufficiently supported bases to 'N's and write out new contig fasta 
    final_contig_names, final_contig_seqs, zero_coverage_base_count_dict = make_final_contigs(ref_dict, pos_dict, weak_dict)
    write_fasta(final_contig_names, final_contig_seqs, '%s_final_contigs.fasta' % opts.out)
    
    #Writing info ot a summary results file
    opts.fout_sum.write('\nLengths of final contigs:\n')
    for index, name in enumerate(final_contig_names):
        opts.fout_sum.write('\t%s\t%d\n' % (name, len(final_contig_seqs[index])))
    
    opts.fout_sum.write('\nNumber of bases with >0 but <%dx of >=BaseQual %d coverage (likely includes some of the %d Ns added to begs and ends):\n' % (opts.finalCov, opts.baseQual, opts.Nfinal))
    for chrom, count in num_bases_lack_qual_cov.items():
        opts.fout_sum.write('\t%s\t%d\n' % (chrom, count))
    
    opts.fout_sum.write('\nNumber of bases without any >=BaseQual %d coverage (likely includes some of the %d Ns added to begs and ends):\n' % (opts.baseQual, opts.Nfinal))
    for chrom, count in zero_coverage_base_count_dict.items():
        opts.fout_sum.write('\t%s\t%d\n' % (chrom, count))
    
    opts.fout_sum.write('\nNumber and Percent of contig bases with various levels of >=BaseQual %d  coverage:\n\t\tCov_Level\tNum_Bases\tPerc_Bases\n' % opts.baseQual)
    for chrom, count_dict in strong_dict.items():
        opts.fout_sum.write('\t%s:\n' % (chrom))
        chrom_length=len(final_contig_seqs[final_contig_names.index(chrom)])
        for num in [opts.finalCov, opts.finalCov*2, opts.finalCov*5, opts.finalCov*10, 50, 100, 500, 1000]:
            if chrom_length==0: opts.fout_sum.write('\t\t%d\t%d\t0.00\n' % (num, count_dict[num]))
            else: opts.fout_sum.write('\t\t%d\t%d\t%.2f\n' % (num, count_dict[num], count_dict[num]/chrom_length*100))
                     

#!!!Should add in the functionality of writing all the completely missing bases to a file
def make_final_contigs(ref_dict, pos_dict, weak_dict):
    new_dict={}
    zero_coverage_base_count_dict={}
    for contig in ref_dict.keys():
        #Copy full sequence into the new_dict
        new_dict[contig]=ref_dict[contig][:]
        
        #Find positions with zero coverage
        missing_pos_list= list(set(range(1, len(ref_dict[contig])+1)) - set(pos_dict[contig]))
        zero_coverage_base_count_dict[contig]=len(missing_pos_list)
        
        #Change any base without sufficient coverage to an 'N'
        for pos in weak_dict[contig] + missing_pos_list:
            new_dict[contig][pos-1]='N'
    
    #Turns seqs back into strings and trim the Ns from the Ends
    for contig, seq_list in new_dict.items():
        new_dict[contig]=''.join(seq_list).strip('N')
    
    return new_dict.keys(), new_dict.values(), zero_coverage_base_count_dict    
    
def calc_covs(bams, ref):
    gen_info_file, chrom_length_dict = make_gen_info_file(ref)
    covs=[]
    for b in bams:
        covs.append(cov_bed(b, gen_info_file))
    return covs, chrom_length_dict, gen_info_file

def make_gen_info_file(reference):
    new_filename = '%s_infofile_for_coverage.txt' % ".".join(reference.split('.')[:-1])
    chrom_lengths={}
    fout = open(new_filename, 'w')
    names, seqs = read_fasta_lists(reference)
    for index in range(len(names)):
        fout.write('%s\t%s\n' % (names[index], len(seqs[index])))
        chrom_lengths[names[index].split()[0]] = len(seqs[index])
    fout.close()
    return new_filename, chrom_lengths

        
def cov_bed(bamfile, gen_info_file):
    coverage_file = '%s.coverage' % ("".join(bamfile.split('.')[:-1]))
    cmd = 'genomeCoverageBed -ibam %s -d -g %s >%s' % (bamfile, gen_info_file, coverage_file)
    cov_proc=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    cov_proc.wait()
    return coverage_file

def find_missing_regions(cov_files, chrom_lengths):
    md={}                        #Will contain info for regions with zero coverage
    #Steps through each organism's coverage file
    for org in cov_files:
        md[org] = {}
        region=[0,0]
        fin = open(org, 'r')
        for line in fin:
            chrom, pos, cov = line.strip().split()
            #This is when you encounter the first base in a region of no coverage
            if cov == '0' and region[0] == 0: 
                region[0] = int(pos)
                region[1] = int(pos)
                if int(pos) == chrom_lengths[chrom]:
                    md[org][chrom] = md[org].get(chrom, [])
                    md[org][chrom].append(region)
                    region=[0,0]
                
                
            #This is when you hit the Nth no coverage base in a row
            elif cov == '0' and region[0] != 0: 
                region[1] = int(pos)
                if int(pos) == chrom_lengths[chrom]:
                    md[org][chrom] = md[org].get(chrom, [])
                    md[org][chrom].append(region)
                    region=[0,0]
            
            #This is when you have hit the end of a no coverage region
            elif cov != '0' and region[0] != 0:
                md[org][chrom] = md[org].get(chrom, [])
                md[org][chrom].append(region)
                region=[0,0]
    return md


def write_missing_regions(md):
        file_names=[]
        #Steps through each strain included in the analysis
        for strain, sd in md.items():
                outname = '%s_missing_regions.txt' % ('.'.join(strain.split('.')[:-1]))
                file_names.append(outname)
                fout = open(outname, 'w')
                for chrom, regions in sd.items():
                        if '|' in chrom:
                                chrom_out_name = chrom.split('|')[-2]
                        else: chrom_out_name = chrom
                        for start, end in regions: 
                                fout.write('%s\t%d\t%d\t%d\n' % (chrom_out_name, start, end, end-start+1))
                fout.close()
        return file_names
        
#Many of the functions in this script are written to allow the user to provide multiple files to process
#Even though for this purpose, we will only have one
def plot_coverage(cov_files):
    plot_names=[]
    for file in cov_files:
        pos_list_dict, cov_list_dict = extract_cov_info(file)
        for chrom in pos_list_dict.keys():
            filename = '%s.pdf' % chrom.replace('|', '_')
            plot_names.append(filename)
            plt.plot(pos_list_dict[chrom], cov_list_dict[chrom], 'r-', linewidth=3)
            plt.axis([1, max(pos_list_dict[chrom]), 0, max(cov_list_dict[chrom])])
            plt.legend(prop=fontP)
            plt.savefig(filename)
            plt.clf()
    return plot_names
            
def extract_cov_info(file):
    pos_dict={}
    cov_dict={}
    fin=open(file, 'r')
    for line in fin:
        chrom, pos, cov = line.strip().split()    
        if chrom not in pos_dict: 
            pos_dict[chrom]=[]
            cov_dict[chrom]=[]
        pos_dict[chrom].append(int(pos))
        cov_dict[chrom].append(int(cov))
    return pos_dict, cov_dict

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (
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


###******************End of coverage section*******************

def join_cap3_contigs(new_cons_fasta):
    sing_names, sing_seqs = read_fasta_lists_simple_names('%s.cap.singlets' % new_cons_fasta)
    joined_names, joined_seqs = read_fasta_lists_simple_names('%s.cap.contigs' % new_cons_fasta)
    #Changing the joined_names slightly to make sure they don't conflict with any of the sing_names
    joined_names=add_random(joined_names)
    all_names=sing_names+joined_names
    all_seqs=sing_seqs+joined_seqs
    cap3_cons_name='%s/cap3_%s' % ('/'.join(new_cons_fasta.split('/')[:-1]), new_cons_fasta.split('/')[-1])
    write_fasta(all_names, all_seqs, cap3_cons_name)
    return cap3_cons_name

def add_random(string_list, num_extra=6):
    import random, string
    new_list=[]
    for each in string_list:
        rand_bit='_'
        counter=0
        while counter < num_extra:
            counter+=1
            rand_bit+=random.choice(string.letters + string.digits)
        new_list.append(each + rand_bit)
    return new_list

def cap3_improved(new_cons_fasta):
    old_names, old_seqs = read_fasta_lists_simple_names(new_cons_fasta)
    new_sing_names, new_sing_seqs = read_fasta_lists_simple_names('%s.cap.singlets' % new_cons_fasta)
    if len(new_sing_names) >= len(old_names): return False
    else: return True 
    
def bowtie_build(opts):

    #Create temporary working direcory and move to that directory
    os.mkdir(opts.temp)
    os.chdir(opts.temp)

    #Add Ns to the beginning and ends of contigs if specified. Ns_ref variable makes it easy to remove this temporary file at the end
    #if opts.Nbeg or opts.Nend:
    Ns_ref=add_Ns(opts.ref, opts.Nbeg, opts.Nend)
    opts.ref=Ns_ref
    #else: Ns_ref='fake_name_that_hopefully_wont_match_a_real_file'
            
    #First alignment of reads to reference with bowtie2
    iter_num=1
    sam = bowtie2_index_align(opts.ref, opts)
    if opts.NOonlyMapped: bam, pileup = make_sort_bam_pileup(sam, opts.ref, opts.mapQ, opts.maxCov)
    else:
        mapped_sam = sam_just_mapped(sam)
        bam, pileup = make_sort_bam_pileup(mapped_sam, opts.ref, opts.mapQ, opts.maxCov)
    num_changes, current_consensus = make_consensus(pileup, opts.ref, opts)    
    print ('%d changes made in iteration %d' % (num_changes, iter_num), file=sys.stderr)
    
    #Continue with aditional iterations as long as changes continue to be made
    while num_changes:
        iter_num+=1
        sam = bowtie2_index_align(current_consensus, opts)
        if opts.NOonlyMapped: bam, pileup = make_sort_bam_pileup(sam, current_consensus, opts.mapQ, opts.maxCov)
        else:
            mapped_sam = sam_just_mapped(sam)
            bam, pileup = make_sort_bam_pileup(mapped_sam, current_consensus, opts.mapQ, opts.maxCov)
        num_changes, current_consensus = make_consensus(pileup, current_consensus, opts)
        print ('%d changes made in iteration %d' % (num_changes, iter_num), file=sys.stderr)
    
    #Create final files, leave working directory and clean it
    trim_Ns(current_consensus)
    new_cons_fasta = '%s/%s.fasta' % (opts.startDir, opts.out)
    os.rename(current_consensus, new_cons_fasta)
    os.rename('change_log.txt', '%s/%s.changelog' % (opts.startDir, opts.out))    
    
    for file in [Ns_ref, '%s.fai' % Ns_ref, sam, mapped_sam, bam, pileup, '%s.fai' % current_consensus, 'ref_index.1.bt2', 'ref_index.2.bt2', 'ref_index.3.bt2', 'ref_index.4.bt2', 'ref_index.rev.1.bt2', 'ref_index.rev.2.bt2',]: 
        if os.path.isfile(file):
            os.remove(file)
    os.chdir(opts.startDir)
    os.rmdir(opts.temp)
    
    #Returns the name of the new consensus fasta file. This should be an absolute path
    return new_cons_fasta

###-------------Start of functions to add and remove Ns from reference---------

def add_Ns(fasta, num_beg, num_end):
    new_fasta_name='%s_%dNbeg_%dNend.fasta' % ('.'.join(fasta.split('.')[:-1]), num_beg, num_end)
    names, seqs = read_fasta_lists_simple_names(fasta)
    seqs=[num_beg*'N'+x+num_end*'N' for x in seqs]
    write_fasta(names, seqs, new_fasta_name)
    return new_fasta_name

def trim_Ns(fasta):
    names, seqs = read_fasta_lists_simple_names(fasta)
    seqs=[x.strip('N') for x in seqs]
    write_fasta(names, seqs, fasta)
    
    
###-------------End of functions to add and remove Ns from reference---------

###----------Start of bowtie2 functions----------------------------------------

##!!! Expecting fastq quality scores to be 33-offset
def bowtie2_index_align(ref, opts):
    
    #Make index for reference
    cmd='bowtie2-build %s ref_index' % (ref)
    print(cmd)
    make_db=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_db.wait()
    
    
    sam_name='aligned.sam'
    
    #Align unpaired reads
    if opts.unpaired and opts.paired1 and opts.paired2:
        cmd='bowtie2 -p %d --phred33 -N 1 --n-ceil %s -x ref_index -U %s -1 %s -2 %s -S %s' % (opts.procs, opts.nceil, opts.unpaired, opts.paired1, opts.paired2, sam_name)
        print(cmd)
        unpaired_bowtie=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        unpaired_bowtie.wait()
    elif opts.unpaired:
        cmd='bowtie2 -p %d --phred33 -N 1 --n-ceil %s -x ref_index -U %s -S %s' % (opts.procs, opts.nceil, opts.unpaired, sam_name)
        print (cmd)
        unpaired_bowtie=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        unpaired_bowtie.wait()
    #Align paired reads    
    elif opts.paired1 and opts.paired2:
        cmd='bowtie2 -p %d --phred33 -N 1 -I %d -X %d --n-ceil %s -x ref_index -1 %s -2 %s -S %s' % (opts.procs, opts.minIns, opts.maxIns, opts.nceil, opts.paired1, opts.paired2, sam_name)
        paired_bowtie=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        paired_bowtie.wait()
    else: 
        print ('Must provide at least 1) an unpaired .fastq file (-u) or 2) two paired end fastq files (-1 and -2)!!!!')
        return
    
    return sam_name

###----------End of bowtie2 functions----------------------------------------

###----------Start of samtools functions----------------------------------------


def make_sort_bam_pileup(sam, ref, mapQ, maxCov, do_pile=True):
    #Make new faidx for current reference
    cmd='samtools faidx %s' % ref
    faidx=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    faidx.wait()
    
    #Make bam from sam, and sort the bam
    sam_prefix='.'.join(sam.split('.')[:-1])
    bam_name='%s_sorted.bam' % sam_prefix
    cmd=f"samtools view -bS {sam} | samtools sort -m 800000000 - -o {bam_name}"
    print(cmd)
    make_bam=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_bam.wait()
    
    if do_pile:
        #Make pileup from bam
        pile_name='%s.pile' % sam_prefix
        cmd='samtools mpileup -f %s -q %d -d %d  %s >%s' % (ref, mapQ, maxCov, bam_name, pile_name)
        print (cmd)
        make_pileup=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        make_pileup.wait() 
        return bam_name, pile_name
        
    else: return bam_name


#This function doesn't actually call samtools, but its related
#Creates a new samfile containing only the header and the information for reads that actually mapped to the referene.
#Uses bit 4 of the flag to check if the read is mapped. 
def sam_just_mapped(samfile):
    sam = open(samfile, 'r')
    new_name='%s_onlymapped.sam' % '.'.join(samfile.split('.')[:-1])
    new = open(new_name, 'w')
    for line in sam:
        if line[0] != '@':
            cols=line.split()
            if not int(cols[1]) & 0x4:new.write(line) 
        else: new.write(line)
    sam.close()
    new.close()
    return new_name
        
###----------End of samtools functions----------------------------------------        


###-----------Start of functions used in making consensus from pileup-----------

def quality_filter(bases, quals, opts):
    new_bases=[]
    new_quals=[]
    for index in range(len(bases)):
        if ord(quals[index])-opts.offset >= opts.baseQual:
            new_bases.append(bases[index])
            new_quals.append(quals[index])
    return new_bases, new_quals


def make_consensus(pileup, ref, opts):
    #Make a dictionary from the reference fasta
    ref_dict=read_fasta_dict_simple_names(ref)
    #Turn each seq into a list 
    for name, seq in ref_dict.items():
        ref_dict[name]=list(seq)
    
    #Open file to append info about changes that are made
    change_log='change_log.txt'
    fout=open(change_log, 'a')
    num_changes=0                    #Will keep track of the number of changes that are made to the consensus
    
    #Will hold positions that have been deleted and therefore the pileup positions for these will be skipped
    to_ignore={}
        
    fin=open(pileup, 'r')
    current_chrom=''                    #Need this so that we can reset base_base each time a new contig is encountered
    for line in fin:

        #Adding a check for whether there are actually any reads at this position
        if line.split("\t")[3]=="0":
            continue

        chrom, pos, bases, indels, ends, quals = pileup_info(line)
        #New step that removes info for bases with quality lower than the specified minimum threshold
        #Only using the quality score cutoff for the SNP changes, not indels
        good_bases, good_quals = quality_filter(bases, quals, opts)
        
        if chrom != current_chrom:
            base_base=-1                #This value changes with insertions and deletions to keep the position to index relationship correct, it starts at -1 b/c pos is 1-based and indices are 0-based
            current_chrom=chrom
            to_ignore={}        #!!!! I think I need to reset to_ignore!!!! This has been added in *_basequaly.py version
            
        #Make sure this position has not been deleted
        if pos not in to_ignore and len(good_bases)>0:
            #Change the base at a given position if it is more than the specified % of bases support an alternative 
            alt=majority_alt(good_bases, opts)
            if alt and len(good_bases)>=opts.covThresh: 
                num_changes+=1
                fout.write('%s\t%d\t%s\t%s\t%.2f\t%d\n' % (chrom, pos, ref_dict[chrom][pos+base_base], alt, good_bases.count(alt)/len(good_bases), len(good_bases)))
                ref_dict[chrom][pos+base_base]=alt
                
        
            #Insert or remove bases if more than the specified % of bases support an indel
            if indels:
                if len(indels)>1:
                    #Will only consider the type of structural change with the most support
                    if len(indels['insert'])>len(indels['delete']): del(indels['delete'])
                    elif len(indels['delete'])>len(indels['insert']): del(indels['insert'])
                    #If they are equal then neither has a chance of being greater than 50% of the whole
                    else: 
                        del(indels['delete'])
                        del(indels['insert'])
                #If some reads show insertions after this base
                if 'insert' in indels:
                    for each in set(indels['insert']):
                        if indels['insert'].count(each)/(len(bases)-ends)>opts.propThresh and (len(bases)-ends)>=opts.covThresh:
                            num_changes+=1
                            fout.write('%s\t%d\t-\t%s\t%.2f\t%d\n' % (chrom, pos, each, indels['insert'].count(each)/(len(bases)-ends), len(bases)-ends))
                            ref_dict[chrom].insert((pos+base_base+1), each)
                            base_base+=1                            #Need to increment this because all subsequent indices have been shifted forwards
                            break    
                #If some reads show deletions after this base
                elif 'delete' in indels:
                    for each in set(indels['delete']):        
                        if indels['delete'].count(each)/(len(bases)-ends)>opts.propThresh and (len(bases)-ends)>=opts.covThresh:
                            pos_to_ignore=pos
                            num_changes+=1
                            fout.write('%s\t%d\t%s\t-\t%.2f\t%d\n' % (chrom, pos, each, indels['delete'].count(each)/(len(bases)-ends), len(bases)-ends))
                            for position in each:
                                pos_to_ignore+=1
                                to_ignore[pos_to_ignore]=''
                                del(ref_dict[chrom][pos+base_base+1])
                                base_base-=1                            #Need to decrement this because all subsequent indices have been shifted backwards
                            break
    
    #Turn each seq into a string
    for name, seq in ref_dict.items():
        ref_dict[name]=''.join(seq)
    
    #Write new consensus fasta file
    new_cons_name = 'new_consensus.fasta'
    write_fasta(list(ref_dict.keys()), list(ref_dict.values()), new_cons_name)
    
    return num_changes, new_cons_name
                        
def majority_alt(bases, opts):
    for b in ['A', 'C', 'G', 'T']:
        if bases.count(b)/len(bases)>opts.propThresh: return b
    #If none of the alts are above the propThresh
    return False
    

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
    indels={}        #This line is needed to ensure an empty dict named indels for strings without indel present
    if '-' in bases or '+' in bases:
        bases, indels = indel_info(bases, indels)
    
    if len(bases)==num_reads: return bases, indels, ends
    else:
        print ('Problem!!! Exp bases: %d, Output bases: %d, Input str: %s, Output str: %s' % (num_reads, len(bases), raw_bases, bases))
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
#    print (bases)
    start_index=bases.index(type_char)
    digitcount=0
    for each in bases[start_index+1:]:
        if each.isdigit(): digitcount+=1
        else: break
    indel_bases=int(bases[start_index+1:start_index+1+digitcount])
#    print (indel_bases)
    indel_str=bases[start_index+1+digitcount:start_index+1+digitcount+indel_bases]
#    print (indel_str)
        
    ind_info.append(indel_str)
    bases=bases[:start_index]+bases[start_index+1+digitcount+indel_bases:]   
    #Recursively go through this process if the string has more indel info
    if type_char in bases: bases, ind_info = characterize_indel(bases, type_char, ind_info)
    
    return bases, ind_info

# will cut fasta name off at the first whitespace
def read_fasta_lists_simple_names(file):
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


#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()

def read_fasta_dict_simple_names(file):
    names, seqs = read_fasta_lists_simple_names(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

###------------END of functions used in building new consensus from pileup----------------------------

###-----------Start of functions used to check the assembly for potential chimeric sections----Can ONLY be run if you have paired data

####################################################################################################
# By Jason Ladner
#
#   !!Should actually be calculating the endpoint of the left read using the CIGAR string, 
#     but the length of the sequence should be good enough, for now at least. 
#
###################################################################################################

def pair_cov_test(sam_file, opts):
    
    opts.fout_sum.write( '\nTesting for potential chimeras using pairs:\n')
    
    ref_dict={}
    fin=open(sam_file, 'r')
    for line in fin:
        #For each line with info on a fastq read 
        if not line.startswith('@'):
            cols=line.strip().split('\t')
            #Checks to make sure that this read belongs to a pair, both mates are mapped and that they are 'properly aligned', according to the aligner
            #Also checks to make sure that this mate is the left-most in respect to the reference, this keeps us from double counting pairs
            if not int(cols[1]) & 0x4 and not int(cols[1]) & 0x8 and int(cols[1]) & 0x2 and int(cols[8])>0 and int(cols[1]) & 0x1:
                count_good_bases(ref_dict, cols[2], int(cols[3]), len(cols[9]), int(cols[7]))
                
        #For each line that has information about the reference sequences
        elif line.startswith('@SQ'): 
            cols=line.strip().split('\t')
            ref_name=cols[1][3:]
            length=int(cols[2][3:])
            ref_dict[ref_name]=[0]*length
                
    for ref, counts in ref_dict.items():
        start=0
        end=0
        for index, count in enumerate(counts):
            base=index+1
            if count==0: 
                if start==0: 
                    start=base
                    end=base
                else:
                    end=base
            #Don't use the nest three lines if you DO want to print out regions at the begining of contigs
            elif start==1:
                start=0
                end=0
                
            elif start!=0: 
                opts.fout_sum.write('\t%s\t%d\t%d\t%d\n' % (ref, start, end, end-start+1))
                start=0
                end=0

#Only use these two lines if you want to print out regions at the end of contigs
#        if start!=0:
#            print ('%s\t%d\t%d\t%d' % (ref, start, end, end-start+1))

def count_good_bases(ref_dict, ref_name, left_start, left_length, right_start):
    for i in range(left_start+left_length-2, right_start-2):
        ref_dict[ref_name][i]+=1
        
        
###------------->>>

if __name__ == "__main__":
    main()
