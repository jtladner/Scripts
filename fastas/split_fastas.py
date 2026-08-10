#!/usr/bin/env python2

import os, optparse
import shutil
import math


def main():
    
    usage = "usage: %prog [options]"
    option_parser = optparse.OptionParser( usage )

    option_parser.add_option(
            '-q', '--query', help = 'Fasta query file, can be multiple sequences'
        )
    option_parser.add_option(
            '-t', '--temp'
        )
    option_parser.add_option(
            '--numProcs', default=4, type=int
        )
        
    options, arguments = option_parser.parse_args()

    input_directory = options.temp
    
    if os.path.exists( input_directory):
        shutil.rmtree( input_directory, ignore_errors=True )

    os.mkdir( input_directory )

    os.chdir( input_directory )

    split_fasta( options )

def split_fasta(opts):
    #Will hold the names of all the files created
    sub_files=[]
    names, seqs = read_fasta_lists(opts.query)
    num_seqs=len(names)
    if num_seqs>=opts.numProcs: sub_size=int(math.ceil(num_seqs/opts.numProcs))                                          #Rounds up so that the first few subsets might have slightly more than the last
    elif num_seqs>0: 
        opts.numProcs=num_seqs
        sub_size=1
    else: return sub_files
    
    for start in range(0, num_seqs, sub_size):
        sub_names=names[start:start+sub_size]
        sub_seqs=seqs[start:start+sub_size]
        new_filename='%d_%d.fasta' % (start+1, start+sub_size)
        sub_files.append(new_filename)
        write_fasta(sub_names, sub_seqs, new_filename)
    return sub_files

def read_fasta_lists(file):
    fin = open( "../" + file, 'r')
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
    
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()



if __name__ == '__main__':
   main() 
