#!/usr/bin/env python

# By Jason Ladner

#This version is for paired end data and will run cutadapt on the reads prior to aligning
#Version 3.1 added prinseq filtering
#Version 3.4 added option for providing multiple comma separated bams along with multiple comma separated outnames
#Version 3.5 just has a small bug fix
#Version 3.6 added support for trimming multiple different sequences in the same way sispa seqs are trimmed.
#Version 4.0:Additional support for the automation of sunning multiple samples, especially for nextseq samples, for which multiple fastqs need to be concatenated
#Version 4.2: Fixed a problem with glob.glod. Added a sorted() command to make sure that the different lanes get concatenated in the same order for R1, R2, I1 and I2. 
#Version 4.4: Added option to use anomalous read pairs
#Version 4.5: Fixed a bug that was causing an issue with multi-base deletions. In previous versions, after deleting the first base, it would start to delete bases upstream instead of additional bases downstream from the first change.
#             Also changed the handling of indels to be more consistent with the handling of SNPs, so now you need to hit the minimum coverage just with reads that support the consensus change
#Version 5.0: Introduced an additional trimming step after cutadapt that removes an additional # bases after a sispa adaptor is trimmed to get rid of the random primer
#Version 5.1: Changed naming of sams, bams, pileups, etc. Also, cleaning up several unwanted files, like the prinseq output, only_mapped_sam, also converting the full sam to a bam and then removing the sam. 
#Version 5.2: Made the location of picard a variable that could be changed on the command line
#Version 5.3: Corrected prinseq command to be -trim_qual_right instead of -trim_right
#Version popgen: changed quality filtering and adaptor clipping defaults, and incorporated a step to call SNVs with freebayes. Includes '14' duplicate removal in prinseq
#Version popgen 1.2: Made picard and optional step
#Version popgen 1.3: Fixed problem with using -i with --Nbeg and/or --Nend

###!!! Requires installations of samtools and bowtie2, both of which must reside within your $PATH


from __future__ import division
import sys, optparse, os, glob
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np

def main():

    #To parse command line
    usage = "usage: %prog [options] bam1 [bam2 ...]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
    p.add_option('-i', '--input', help='Input file for running multiple samples from fastqs [None]')
    p.add_option('-u', '--unpaired', help='Fastq file with unpaired reads to be mapped to reference [None]')
    p.add_option('-1', '--paired1', help='Fastq file with read 1 of paired reads [None]')
    p.add_option('-2', '--paired2', help='Fastq file with read 2 of paired reads [None]')
    p.add_option('-r', '--ref', help='Fasta file that was used as reference in mapping. [None, Required]')
    p.add_option('-o', '--out', default='new_consensus', help='Base name for output files. Output files will be written to the current working directory. [new_consensus]')
    p.add_option('-s', '--sam', help='Optional starting place [None]')
    p.add_option('-b', '--bam', help='Optional starting place [None]')
    p.add_option('-p', '--pile', help='Optional starting place [None]')

    
    #General
    p.add_option('--recursThresh', type='int', default=50000, help='Minimum level of coverage required to change consensus. [3000]')
    p.add_option('--temp', help='Name for working directory that will hold internediate files. This is not currently being cleaned [based on -o]')
    p.add_option('--procs', type='int', default=2, help='Number of processors to use in multi-threaded portions [2]')
    p.add_option('--Nbeg', type='int', default=0, help='# of Ns to add to the beginning of each contig [0]')
    p.add_option('--Nend', type='int', default=0, help='# of Ns to add to the beginning of each contig [0]')    
    
    #For samtools
    p.add_option('--mapQ', type='int', default=10, help='Minimum mapping quality for a read to be used in the pileup generation [10]')
    p.add_option('--maxCov', type='int', default=300, help='Max per base coverage to be used in the pileup generation [300]')
    #This step allows for a large amount of time saved when a small number of the reads actually map to the reference. Should probably be used for viral samples, but not bacterial.
    p.add_option('--NOonlyMapped', default=False, action='store_true', help='Use this flag to turn off the part of the script that creates a smaller sam with only mapped reads prior to sorting')
    
    #For bowtie2
    p.add_option('--nceil', default='0,0.3',  help='Function to supply as --n-ceil parameter in bowtie2. [0,0.3]')
    p.add_option('--minIns', type='int', default=0, help='Minimum insert size for a pair to be concordant [0]')
    p.add_option('--maxIns', type='int', default=500, help='Maximum insert size for a pair to be concordant [500]')
    
    #Make new consensus
    p.add_option('--propThresh', type='float', default=0.5, help='Requires that a greater proportion of the reads (than this threshold) must support an alternative for the consensus to be changed. [0.5]')
    p.add_option('--covThresh', type='int', default=5, help='Minimum level of coverage required to avoid an N in the consensus. [5]')
    p.add_option('--offset', type='int', default=33, help='Base quality offset used in the pileup. I believe the default in samtools is Sanger(33) [33]')
    p.add_option('--baseQual', type='int', default=20, help='Minimum base quality for a base to be used in the evalution of the consensus. [20]')
    
    #Cutadapt and Prinseq
    #How to call programs
    p.add_option('--ca', default='cutadapt', help='How to call cutadapt. Need to use version >=1.2 [cutadapt]')
    p.add_option('--ps', default='prinseq-lite_0.20.3.pl', help='How to call prinseq. [prinseq-lite-0.20.3.pl]')
    p.add_option('--pc', default='~/programs/picard.jar', help='How to call picard.jar. [~/programs/picard.jar]')
#    p.add_option('--tr', default='~/programs/trimmomatic-0.33.jar', help='How to call trimmomatic. [~/programs/trimmomatic-0.33.jar]')
    p.add_option('--i1', default='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', help="Illumina adaptor to be trimmed from the 3' end of R1. [GATCGGAAGAGCACACGTCTGAACTCCAGTCAC]")
    p.add_option('--i2', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', help='Reverse complement of the Illumina adaptor that will be trimmed from the 3; end of R2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    p.add_option('-l', '--minLength', default=70, type='int', help='Minimum length to keep (after primer clipping) [70]')
    p.add_option('--trimQual', default=30, type='int', help='Minimum quality to keep at 3 prime end f seqs [30]')
    p.add_option('--trimType', default='min', help='Type of quality score calculation to use. Allowed options are min, mean, max and sum [min]')
    p.add_option('--trimWin', default=5, type='int', help='Size of window to use for quality trimming [5]')
    p.add_option('-q', '--meanQual', default=20, type='int', help='Minimum mean quality for a read to be kept [20]')
    p.add_option('--sispa', help='For Sispa , use GCCGGAGCTCTGCAGATATC. Add GCCGGAGCTCTGCAGATATCGGCCATTATGGCCGGG and GCCGGAGCTCTGCAGATATCAAAAAAAAAAAAA for Sispa RACE. For kikwit amplicons, add TGTAAAACGACGGCCAGT. Sispa adapter to remove from the seqs. Rev comp of this sequence will be removed from the 3 prime end [None]')
    #Use this flag if you are starting with reads, but you don't want to run cutadapt
    p.add_option('--indexFilt', type='int', help='Use this to specify the quality cutoff for index filtering. Currently only supported in combination with -i [None].')
    p.add_option('--useAnomPairs', default=False, action='store_true', help='Use this flag to use anomalous read pairs in pileup [False]')
    #Number of bases to hard clip from the beginning of R1, end of R2, to get rid of any potential random hexamers
    p.add_option('--hc', default=6, type='int', help='Number of bases to be hard clipped from the beginning of R1 and end of R2. This will be done with prinseq [6]')

    #To turn off different parts of the pipeline
#    p.add_option('--NOtrimmomatic', default=False, action='store_true', help='Use this flag to turn off the part of the script that runs trimmomatic')
    p.add_option('--NOcutadapt', default=False, action='store_true', help='Use this flag to turn off the part of the script that runs cutadapt')
    p.add_option('--NOprinseq', default=False, action='store_true', help='Use this flag to turn off the part of the script that runs prinseq')
    p.add_option('--NOpicard', default=False, action='store_true', help='Use this flag to turn off the part of the script that runs picard')


    opts, args = p.parse_args()

    opts.index1=''
    opts.index2=''
    
    #Increase recursion limit before starting script
    sys.setrecursionlimit(opts.recursThresh)
    
    #To allow for trimming multiple sispa like sequences (trimmed from both sides of reads), this function replaces the sequences with the string that will go into the cutadapt command
    if opts.sispa:
        opts.sispa=make_sispa_str(opts.sispa)
    

    #If providing multiple bams as args. Output names will be taken from bam names
    if args:
        for bamfile in args:
            opts.bam=bamfile
            opts.out='.'.join(bamfile.split('.')[:-1])
            opts.temp=None
            run_iteration(opts)


    #If you are running multiple samples from fastqs
    elif opts.input:
        finin = open(opts.input, 'r')
        start_ref=opts.ref
        for line in finin:
            #This resting of the ref is important if you are adding Ns onto the ends of the ref
            opts.ref=start_ref
            cols=line.strip().split('\t')
            opts.out=cols[0]
            print cols[1]
            opts.paired1=','.join(sorted(glob.glob('%s*R1*' % cols[1])))
            opts.paired2=','.join(sorted(glob.glob('%s*R2*' % cols[1])))
            if opts.indexFilt:
                opts.index1=','.join(sorted(glob.glob('%s*I1*' % cols[1])))
                opts.index2=','.join(sorted(glob.glob('%s*I2*' % cols[1])))
            if len(cols)>2: opts.ref=cols[2]
            print opts.paired1, opts.paired2, opts.ref
            run_iteration(opts)

    #If you are running multiple bam files
    elif ',' in opts.out:
        out_list=opts.out.split(',')
        if opts.bam:
            bam_list=opts.bam.split(',')
            if len(bam_list) == len(out_list):
                for i in range(len(out_list)):
                    opts.out=out_list[i]
                    opts.bam=bam_list[i]
                    opts.temp=None
                    run_iteration(opts)
            else: print 'Must provide the same number of input files and output names!!!!'
    
    #If you are running just one bam file or something else
    else: run_iteration(opts)

###-----------------End of main()--------------------------->>>

def make_sispa_str(seqs_comma):
    cut_str=''
    for seq in seqs_comma.split(','):
        cut_str+=' -g sispaF=%s -a sispaR=%s ' % (seq, rev_comp(seq))
    return cut_str

def concat_fastqs(fastq_str, newname):
    fastqs=fastq_str.split(',')
    cmd='zcat %s >%s' % (' '.join(fastqs), newname)
    print cmd
    concat=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    concat.wait()
    return newname


def run_iteration(opts):
    #Concatenate fastqs, if multiple
    #And set to absolute path
    if opts.paired1:
        if ',' in opts.paired1:
            opts.paired1=concat_fastqs(opts.paired1, '%s_R1.fastq' % (opts.out))
        opts.paired1=os.path.abspath(opts.paired1)
    if opts.paired2:
        if ',' in opts.paired2:
            opts.paired2=concat_fastqs(opts.paired2, '%s_R2.fastq' % (opts.out))
        opts.paired2=os.path.abspath(opts.paired2)
    
    if opts.index1:
        if ',' in opts.index1:
            opts.index1=concat_fastqs(opts.index1, '%s_I1.fastq' % (opts.out))
        opts.index1=os.path.abspath(opts.index1)
    if opts.index2:
        if ',' in opts.index2:
            opts.index2=concat_fastqs(opts.index2, '%s_I2.fastq' % (opts.out))
        opts.index2=os.path.abspath(opts.index2)
    
    #Index filter, if requested
    if opts.indexFilt:
        opts.paired1, opts.paired2 = filter_index(opts)
        opts.paired1=os.path.abspath(opts.paired1)
        opts.paired2=os.path.abspath(opts.paired2)

    #Change all filenames to their absolute path versions
    opts.ref=os.path.abspath(opts.ref)
    if opts.unpaired: opts.unpaired=os.path.abspath(opts.unpaired)
    if opts.pile: opts.pile=os.path.abspath(opts.pile)
    if opts.sam: opts.sam=os.path.abspath(opts.sam)
    if opts.bam: opts.bam=os.path.abspath(opts.bam)
    
    #If no name was provied for the temp directory, make a name for it based on the -o flag
    if not opts.temp: opts.temp='%s_interfiles' % opts.out

    #Save current working directory
    opts.startDir=os.getcwd()
    
    #Create consensus
    create_consensus(opts)
    
    #Added this line so that if running multiple analyses at once
    opts.temp=None


#Index filter scripts--------------------

def filter_index(opts):
    
    #Open new output files
    fout_r1 = open('%s_R1_Q%d.fastq' % (opts.out, opts.indexFilt), 'w')
    fout_r2 = open('%s_R2_Q%d.fastq' % (opts.out, opts.indexFilt), 'w')
#    fout_i1 = open('%s_Q%d.fastq' % (opts.I1.split('.')[0], opts.minQual), 'w')
#    fout_i2 = open('%s_Q%d.fastq' % (opts.I2.split('.')[0], opts.minQual), 'w')
    #Open input fastq files
    fin_r1 = open(opts.paired1, 'r')
    fin_r2 = open(opts.paired2, 'r')
    fin_i1 = open(opts.index1, 'r')
    fin_i2 = open(opts.index2, 'r')
    
#    good_index_count=0
    while 1:
        first_mate=[]
        second_mate=[]
        first_index=[]
        second_index=[]
        #Extract a single sequence from each read file
        for i in range(4):
            first_mate.append(fin_r1.readline())
            second_mate.append(fin_r2.readline())
            first_index.append(fin_i1.readline())
            second_index.append(fin_i2.readline())
            
        if first_mate[0]=='': 
            #Close all open files
            fout_r1.close()
            fout_r2.close()
#            fout_i1.close()
#            fout_i2.close()
            fin_r1.close()
            fin_r2.close()
            fin_i1.close()
            fin_i2.close()
            return '%s_R1_Q%d.fastq' % (opts.out, opts.indexFilt), '%s_R2_Q%d.fastq' % (opts.out, opts.indexFilt)
        else: 
            if good_qual(first_index[3], opts.indexFilt, opts.offset) and good_qual(second_index[3], opts.indexFilt, opts.offset):
#                good_index_count+=1
                fout_r1.writelines(first_mate)
                fout_r2.writelines(second_mate)
#                fout_i1.writelines(first_index)
#                fout_i2.writelines(second_index)
#    print 'Number of good pairs remaining: %d' % (good_index_count)

    return '%s_R1_Q%d.fastq' % (opts.out, opts.indexFilt), '%s_R2_Q%d.fastq' % (opts.out, opts.indexFilt)

def good_qual(index_quals, min_qual, offset):
    index_quals=[ord(x)-offset for x in index_quals.strip()]
    if np.average(index_quals)>=min_qual: return True
    else: return False
    
###---------End of index filter scripts---------------

def quality_filter(bases, quals, opts):
    new_bases=[]
    new_quals=[]
    for index in range(len(bases)):
        if ord(quals[index])-opts.offset >= opts.baseQual:
            new_bases.append(bases[index])
            new_quals.append(quals[index])
    return new_bases, new_quals


def create_consensus(opts):

    #Create temporary working direcory and move to that directory
    os.mkdir(opts.temp)
    os.chdir(opts.temp)

    if opts.Nbeg or opts.Nend:
        Ns_ref=add_Ns(opts.ref, opts.Nbeg, opts.Nend)
        opts.ref=Ns_ref
    else: Ns_ref='fake_name_that_hopefully_wont_match_a_real_file'

    
    #If you are staring with reads that need to be mapped
    if not opts.pile and not opts.bam and not opts.sam:
        
        #Not adding trimmomatic just now
        #Run Trimmomatic, unless NOtrimmomatic flag was thrown
#        if not opts.NOtrimmomatic:
            
        
        #Run Cutadapt, unless NOcutadapt flag was thrown
        if not opts.NOcutadapt:
            #Run cutadapt
            one_cut, one_info_file = run_cutadapt(opts.paired1, opts.i1, opts)
            two_cut, two_info_file = run_cutadapt(opts.paired2, opts.i2, opts)

            #Delete random primers associated with sispa adapters
#            if opts.num:
#                one_cut_trim=trim_fastq(one_cut, one_info_file, opts)
#                two_cut_trim=trim_fastq(two_cut, two_info_file, opts)
#                delete_files(one_cut, two_cut, one_info_file, two_info_file)
#                one_cut=one_cut_trim
#                two_cut=two_cut_trim

            #Make new fastqs with just reads that are still paired. Delete intermediate files
            opts.paired1, opts.paired2 = paired_sort(one_cut, two_cut)
            opts.paired1=os.path.abspath(opts.paired1)
            opts.paired2=os.path.abspath(opts.paired2)
            delete_files(one_cut, two_cut)

        #Hard clip a defined number of bases from the beginning of R1 and the end of R2, if specified
        #This is to make sure to get rid of any random mer incorporated during library prep
        if opts.hc:
            opts.paired1, opts.paired2 = hard_clip(opts)

        #Run prinseq, unless NOprinseq flag was thrown
        if not opts.NOprinseq:
            opts.paired1, opts.paired2 = paired_prinseq(opts)
        
        #Add Ns to the beginning and ends of contigs if specified. Ns_ref variable makes it easy to remove this temporary file at the end
        #First alignment of reads to reference with bowtie2
        sam = bowtie2_index_align(opts.ref, opts)
    
    if not opts.pile and not opts.bam:
        if opts.sam: sam=opts.sam
        
        if opts.NOonlyMapped:
            bam, pileup = make_sort_bam_pileup(sam, opts.ref, opts)
            delete_files(sam)
        else:
            mapped_sam = sam_just_mapped(sam)
            bam, pileup = make_sort_bam_pileup(mapped_sam, opts.ref, opts)
            num_changes, current_consensus = make_consensus(pileup, opts.ref, opts)    
            sam_to_bam(sam)
            delete_files(mapped_sam, sam)
            
    elif not opts.pile and opts.bam:
        bam, pileup = sort_bam_pileup(opts.bam, opts.ref, opts)
        num_changes, current_consensus = make_consensus(pileup, opts.ref, opts)    

    elif opts.pile:
        pileup=opts.pileup
        num_changes, current_consensus = make_consensus(pileup, opts.ref, opts)    
        
    #Create final files, leave working directory and clean it
    new_cons_fasta = '%s/%s.fasta' % (opts.startDir, opts.out)
    os.rename(current_consensus, new_cons_fasta)
    os.rename('change_log.txt', '%s/%s.changelog' % (opts.startDir, opts.out))    
    
    for file in ['ref_index.1.bt2', 'ref_index.2.bt2', 'ref_index.3.bt2', 'ref_index.4.bt2', 'ref_index.rev.1.bt2', 'ref_index.rev.2.bt2',]: 
        if os.path.isfile(file):
            os.remove(file)
    os.chdir(opts.startDir)
#    os.rmdir(opts.temp)

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

def hard_clip(opts):
    #Clip the specified # bases from the beginning of R1
    out_name1 = '%s_trim%dbeg' % ('.'.join(opts.paired1.split('/')[-1].split('.')[:-1]), opts.hc)
    cmd='%s -fastq %s -out_good %s -out_bad null -trim_left %d' % (opts.ps, opts.paired1, out_name1, opts.hc)
    print cmd
    ps=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    ps.wait()
    for line in ps.stderr:
        print line.rstrip()
    
    #Clip the specified # bases from the end of R2
    out_name2 = '%s_trim%dend' % ('.'.join(opts.paired2.split('/')[-1].split('.')[:-1]), opts.hc)
    cmd='%s -fastq %s -out_good %s -out_bad null -trim_right %d' % (opts.ps, opts.paired2, out_name2, opts.hc)
    print cmd
    ps=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    ps.wait()
    for line in ps.stderr:
        print line.rstrip()

    return os.path.abspath("%s.fastq" % out_name1), os.path.abspath("%s.fastq" % out_name2)

###----------Start of Cutadapt and Prinseq functions----------------------------------------

def rev_comp(seq):
    dna = Seq(seq, generic_dna)
    return str(dna.reverse_complement())

#Use this to delete some files, given that those files exist    
def delete_files(*done_files):
    for to_remove in done_files:
        if os.path.isfile(to_remove):
            os.remove(to_remove)

def run_cutadapt(fastq, ill_adapt, opts):
    if fastq.endswith('.gz'): out_name = '%s_cutadapt.fastq' % '.'.join(fastq.split('.')[:-2])
    else: out_name = '%s_cutadapt.fastq' % '.'.join(fastq.split('.')[:-1])
    info_file_name =  '%s_ca_info.txt' % '.'.join(fastq.split('.')[:-1])
    if opts.sispa: cmd='%s -a illumina=%s %s -o %s -m %d --info-file %s --times 3 --match-read-wildcards %s --trim-n' % (opts.ca, ill_adapt, opts.sispa, out_name, opts.minLength, info_file_name , fastq)
    else: cmd='%s -a illumina=%s -o %s -m %d --info-file %s --times 3 --match-read-wildcards %s --trim-n' % (opts.ca, ill_adapt, out_name, opts.minLength, info_file_name , fastq)
    print cmd
    cutit=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    cutit.wait()
    for line in cutit.stdout:
        print line.rstrip()
    return out_name, info_file_name

###-----Functions for removing random primers after running cutadapt with the sispa seqs

def trim_fastq(fq, ifile, opts):

    info_dict=make_info_dict(ifile)
    
    fin=open(fq, 'r')
    fout=open('%s_trimmed' % fq, 'w')
    
    while 1:
        this_read = []
        for i in range(4):
            this_read.append(fin.readline())
        if this_read[0]=='': 
            #Close all open files
            fout.close()
            fin.close()
            return '%s_trimmed' % fq
        else:
            name=this_read[0].strip()[1:]
            if name in info_dict: 
                if 'sispaF' in info_dict[name]:
                    this_read[1]=this_read[1][opts.num:]
                    this_read[3]=this_read[3][opts.num:]
                if 'sispaR' in info_dict[name]:
                    this_read[1]='%s\n' % this_read[1].strip()[-opts.num:]
                    this_read[3]='%s\n' % this_read[3].strip()[-opts.num:]
                fout.writelines(this_read)
            else: fout.writelines(this_read)

    
def make_info_dict(info_file):
    idict={}
    fin=open(info_file, 'r')
    for line in fin:
        cols=line.strip().split('\t')
        if cols[0] not in idict: idict[cols[0]]=[]
        idict[cols[0]].append(cols[7])
    fin.close()
    return idict

###-------------------------------------------------------

def paired_prinseq(opts):
    out_name = '%s_good' % '.'.join(opts.paired1.split('/')[-1].split('.')[:-1])
    cmd='%s -fastq %s -fastq2 %s -out_good %s -out_bad null -min_len %d -lc_method dust -lc_threshold 3 -derep 14 -trim_qual_right %d -trim_qual_type %s -trim_qual_window %d -min_qual_mean %d' % (opts.ps, opts.paired1, opts.paired2, out_name, opts.minLength, opts.trimQual, opts.trimType, opts.trimWin, opts.meanQual)
    print cmd
    ps=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    ps.wait()
    for line in ps.stderr:
        print line.rstrip()
    
    #Deleting singleton files because they are not being used
    delete_files('%s_1_singletons.fastq' % out_name, '%s_2_singletons.fastq' % out_name)
    
    return '%s_1.fastq' % out_name, '%s_2.fastq' % out_name
    
#For sifting through fastqs and just keeping what is still paired

def paired_sort(read1, read2):
    names1 = get_names(read1)
    names2 = get_names(read2)
    paired = set(names1) & set(names2)

    del names1
    del names2

    pair1_file = write_new_file(read1, paired)
    pair2_file = write_new_file(read2, paired)
    
    return pair1_file, pair2_file

def get_names(file):
    fin = open(file, 'r')	
    names=[]
    linenum=0

    for line in fin:
        linenum+=1
        #First name line
        if linenum%4==1:
            names.append(line.strip().split()[0])
    fin.close()
    return names

def write_new_file(fastq, paired_names):

    fin = open(fastq, 'r')
    fout_pair = open('.'.join(fastq.split('.')[:-1]) + '_paired.fastq', 'w')
    linenum=0
    is_paired=0

    for line in fin:
        linenum+=1
        #First name line
        if linenum%4==1:
            name=line.strip().split()[0]
            if name in paired_names:
                is_paired=1
                fout_pair.write(line)	
            else:
                is_paired=0
        #Other lines
        else:
            if is_paired: fout_pair.write(line)
    fin.close()
    fout_pair.close()
    return '.'.join(fastq.split('.')[:-1]) + '_paired.fastq'


###----------Start of samtools functions----------------------------------------

def sam_to_bam(sam):
    sam_prefix='.'.join(sam.split('/')[-1].split('.')[:-1])
    bam_name='%s.bam' % sam_prefix
    cmd='samtools view -bS %s >%s' % (sam, bam_name)
    print cmd
    make_bam=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_bam.wait()
    

def make_sort_bam_pileup(sam, ref, opts):
    #Make new faidx for current reference
    cmd='samtools faidx %s' % ref
    print cmd
    faidx=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    faidx.wait()
    
    #Make bam from sam, and sort the bam
    sam_prefix='.'.join(sam.split('/')[-1].split('.')[:-1])
    bam_name='%s_sorted.bam' % sam_prefix
    cmd='samtools view -bS %s | samtools sort -m 800000000 - %s' % (sam, sam_prefix + '_sorted')
    print cmd
    make_bam=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_bam.wait()
    
#***NEW in _picard version ***    
    #Make a new version of the bam with duplicates removed using the MarkDuplicates function in picard
    if not opts.NOpicard:
        dup_bam_name='%s_sorted_nodups.bam' % sam_prefix
        dup_cmd="java -XX:ParallelGCThreads=%s -jar %s MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true O=%s M=%s_metrics.txt I=%s" % (opts.procs, opts.pc, dup_bam_name, dup_bam_name, bam_name)
        print dup_cmd
        rmdups=Popen(dup_cmd, shell=True, stdout=PIPE, stderr=PIPE)
        rmdups.wait()
        bam_name=dup_bam_name

    
    #Make pileup from bam
    #I added the -B to prevent the altering of the quality values
    pile_name='%s.pile' % sam_prefix
    if opts.useAnomPairs: cmd='samtools mpileup -f %s -q %d -d %d  -BA %s >%s' % (ref, opts.mapQ, opts.maxCov, bam_name, pile_name)
    else: cmd='samtools mpileup -f %s -q %d -d %d  -B %s >%s' % (ref, opts.mapQ, opts.maxCov, bam_name, pile_name)
    print cmd
    make_pileup=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_pileup.wait()
    
    return bam_name, pile_name
#***End of NEW ***


def sort_bam_pileup(bam, ref, opts):
    #Make new faidx for current reference
    cmd='samtools faidx %s' % ref
    print cmd
    faidx=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    faidx.wait()
    
    #Sort the bam
    bam_prefix='.'.join(bam.split('/')[-1].split('.')[:-1])
    bam_name='%s_sorted.bam' % bam_prefix
    cmd='samtools sort -m 800000000 %s %s' % (bam, bam_prefix + '_sorted')
    print cmd
    make_bam=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_bam.wait()
    
#***NEW in _picard version ***    
    #Make a new version of the bam with duplicates removed using the MarkDuplicates function in picard
    if not opts.NOpicard:
        dup_bam_name='%s_sorted_nodups.bam' % bam_prefix
        dup_cmd="java -XX:ParallelGCThreads=%s -jar %s MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true O=%s M=%s_metrics.txt I=%s" % (opts.procs, opts.pc, dup_bam_name, dup_bam_name, bam_name)
        print dup_cmd
        rmdups=Popen(dup_cmd, shell=True, stdout=PIPE, stderr=PIPE)
        rmdups.wait()
        bam_name = dup_bam_name

#*** Also needed to change the bam being used in the pileup formation
    
    #Make pileup from bam
    pile_name='%s.pile' % bam_prefix
    if opts.useAnomPairs: cmd='samtools mpileup -f %s -q %d -d %d -BA %s >%s' % (ref, opts.mapQ, opts.maxCov, bam_name, pile_name)
    else: cmd='samtools mpileup -f %s -q %d -d %d  -B %s >%s' % (ref, opts.mapQ, opts.maxCov, bam_name, pile_name)
    print cmd
    make_pileup=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    make_pileup.wait()
    
    return bam_name, pile_name
#***End of NEW ***


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

#Doesn't actually return anything, just changes zero coverage bases to Ns
def zero_to_Ns(ref_dict, pileup):
    #Make dict of bases in pileup
    pile_dict={}
    fin=open(pileup, 'r')
    for line in fin:
        cols=line.strip().split('\t')
        if cols[0] not in pile_dict: pile_dict[cols[0]]=[]
        pile_dict[cols[0]].append(int(cols[1]))
    fin.close()

    #Find missing bases and replace them with Ns
    for chrom, base_list in ref_dict.iteritems():
        if chrom in pile_dict:
            missing_bases = set(range(1,len(base_list)+1))-set(pile_dict[chrom])
            for each in missing_bases: ref_dict[chrom][each-1]='N'
        #If the whole chrom has no coverage
        else:
            ref_dict[chrom]= ['N']*len(ref_dict[chrom])
    
def make_consensus(pileup, ref, opts):
    #Make a dictionary form the reference fasta
    ref_dict=read_fasta_dict_simple_names(ref)
    #Turn each seq into a list 
    for name, seq in ref_dict.iteritems():
        ref_dict[name]=list(seq)
    
    #Replace all 0 coverage bases (not included in the pileup) to Ns
    zero_to_Ns(ref_dict, pileup)
    
    #Open file to append info about changes that are made
    change_log='change_log.txt'
    fout=open(change_log, 'a')
    num_changes=0					#Will keep track of the number of changes that are made to the consensus
    
    #Will hold positions that have been deleted and therefore the pileup positions for these will be skipped
    to_ignore={}
        
    fin=open(pileup, 'r')
    current_chrom=''					#Need this so that we can reset base_base each time a new contig is encountered
    for line in fin:
        chrom, pos, bases, indels, ends, quals = pileup_info(line)
        
        #Extract only 'good' quality bases
        good_bases, good_quals = quality_filter(bases, quals, opts)

        if chrom != current_chrom:
            base_base=-1				#This value changes with insertions and deletions to keep the position to index relationship correct, it starts at -1 b/c pos is 1-based and indices are 0-based
            current_chrom=chrom
        
        #Make sure this position has not been deleted
        if pos not in to_ignore:
            #Change the base at a given position if it is more than the specified % of bases support an alternative
            #Also change to an N if less that opts.covThresh coverage for any particular base 
            if len(good_bases)<1:alt='N'
            else: alt=majority_alt(good_bases, opts)
            if alt: 
                num_changes+=1
                ref_base=ref_dict[chrom][pos+base_base]
                if len(good_bases)==0: fout.write('%s\t%d\t%s\t%s\t%d\t%.2f\t%.2f\t%d\n' % (chrom, pos, ref_base, alt, len(bases), bases.count(ref_base)/len(bases), bases.count(alt)/len(bases), len(good_bases)))
                else: fout.write('%s\t%d\t%s\t%s\t%d\t%.2f\t%.2f\t%d\t%.2f\t%.2f\n' % (chrom, pos, ref_base, alt, len(bases), (bases.count('.')+bases.count(','))/len(bases), bases.count(alt)/len(bases), len(good_bases), (good_bases.count('.')+good_bases.count(','))/len(good_bases), good_bases.count(alt)/len(good_bases)))
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
                        if indels['insert'].count(each)/(len(bases)-ends)>opts.propThresh and (indels['insert'].count(each)-ends)>=opts.covThresh:
                            num_changes+=1
                            #!!!!!This output has not been updated to match the changes made in the SNP output
                            fout.write('%s\t%d\t-\t%s\t%.2f\t%d\n' % (chrom, pos, each, indels['insert'].count(each)/(len(bases)-ends), len(bases)-ends))
                            ref_dict[chrom].insert((pos+base_base+1), each)
                            base_base+=1							#Need to increment this because all subsequent indices have been shifted forwards
                            break    
                #If some reads show deletions after this base
                elif 'delete' in indels:
                    for each in set(indels['delete']):        
                        if indels['delete'].count(each)/(len(bases)-ends)>opts.propThresh and (indels['delete'].count(each)-ends)>=opts.covThresh:
                            pos_to_ignore=pos
                            num_changes+=1
                            #!!!!!This output has not been updated to match the changes made in the SNP output
                            fout.write('%s\t%d\t%s\t-\t%.2f\t%d\n' % (chrom, pos, each, indels['delete'].count(each)/(len(bases)-ends), len(bases)-ends))
                            for position in each:
                                pos_to_ignore+=1
                                to_ignore[pos_to_ignore]=''
                                del(ref_dict[chrom][pos+base_base+1])
                            base_base-=len(each)							#Need to decrement this because all subsequent indices have been shifted backwards
                            break
    
    #Turn each seq into a string
    for name, seq in ref_dict.iteritems():
        ref_dict[name]=''.join(seq)
    
    #Write new consensus fasta file
    new_cons_name = 'new_consensus.fasta'
    write_fasta(["%s_%s" % (opts.out ,x) for x in ref_dict.keys()], ref_dict.values(), new_cons_name)
    
    return num_changes, new_cons_name
                        
def majority_alt(bases, opts):
    for b in ['A', 'C', 'G', 'T']:
        if bases.count(b)/len(bases)>opts.propThresh and bases.count(b)>=opts.covThresh: return b
    #If none of the alts are above the propThresh
    if bases.count('.') + bases.count(',')<opts.covThresh: return 'N'
    return False
    

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
    indels={}		#This line is needed to ensure an empty dict named indels for strings without indel present
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
#    print bases
    start_index=bases.index(type_char)
    digitcount=0
    for each in bases[start_index+1:]:
        if each.isdigit(): digitcount+=1
        else: break
    indel_bases=int(bases[start_index+1:start_index+1+digitcount])
#    print indel_bases
    indel_str=bases[start_index+1+digitcount:start_index+1+digitcount+indel_bases]
#    print indel_str
        
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

###----------Start of bowtie2 functions----------------------------------------

##!!! Expecting fastq quality scores to be 33-offset
def bowtie2_index_align(ref, opts):
    
    #Make index for reference
    cmd='bowtie2-build %s ref_index' % (ref)
    print cmd
    os.popen(cmd)
    
    sam_name='%s_aligned.sam' % opts.out
    
    #Align unpaired reads
    if opts.unpaired and opts.paired1 and opts.paired2:
        cmd='bowtie2 -p %d --phred33 -N 1 --n-ceil %s -x ref_index -U %s -1 %s -2 %s -S %s' % (opts.procs, opts.nceil, opts.unpaired, opts.paired1, opts.paired2, sam_name)
        print cmd
        unpaired_bowtie=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        unpaired_bowtie.wait()
    elif opts.unpaired:
        cmd='bowtie2 -p %d --phred33 -N 1 --n-ceil %s -x ref_index -U %s -S %s' % (opts.procs, opts.nceil, opts.unpaired, sam_name)
        print cmd
        unpaired_bowtie=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        unpaired_bowtie.wait()
    #Align paired reads    
    elif opts.paired1 and opts.paired2:
        cmd='bowtie2 -p %d --phred33 -N 1 -I %d -X %d --n-ceil %s -x ref_index -1 %s -2 %s -S %s' % (opts.procs, opts.minIns, opts.maxIns, opts.nceil, opts.paired1, opts.paired2, sam_name)
        print cmd
        paired_bowtie=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        paired_bowtie.wait()
    else: 
        print 'Must provide at least 1) an unpaired .fastq file (-u) or 2) two paired end fastq files (-1 and -2)!!!!'
        return
    
    return sam_name

###----------End of bowtie2 functions----------------------------------------

        
        
###------------->>>

if __name__ == "__main__":
    main()
