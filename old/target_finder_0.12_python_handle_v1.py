#!/usr/bin/python

# Usage: targetfinder.pl -s <sequence> -d <target database> [OPTIONS]
# REQUIRED:
# -s     Small RNA sequence
# -d     Target database path (FASTA-format)
# OPTIONAL:
# -q     Query sequence name (DEFAULT = query)
# -c     Score cutoff value (DEFAULT = 4)
# -r     Search reverse strand for targets? (BOOLEAN, DEFAULT=FALSE)
# -h     Show this menu
#
# Type perldoc targetfinder.pl for more he

import subprocess
import os

##I/O files
#1. miRNA fasta file, make sure the fasta header is clean: use sftp://kakrana@pikachu.dbi.udel.edu/data2/homes/kakrana/src/clean_fasta_header_v1.py
#2. cDNA or other RNA type fasta file, make sure the fasta header is clean: use sftp://kakrana@pikachu.dbi.udel.edu/data2/homes/kakrana/src/clean_fasta_header_v1.py
#3. File to write target finder results


##1 miRBASE File
fh_in=open('miRinput.fa_new_head.fa', 'r')
##2 cDNA File
cDNA_file = 'TAIR10_cdna_20101214.cdna_new_head.fa'
##3
#fh_out=open('target_finder_results', 'w')

cutoff = str(7)

##Open a directory to store all target files
os.mkdir('./target_finder_results')



##Read file
mir_base = fh_in.read()
mir_blocks= mir_base.split('>')

###Main


for i in mir_blocks[1:]:
    block = i.strip('\n')##Remove newline from end of entry
    ent = block.split('\n')##Use the newline between header and sequence to split
#    print (ent)
    query = ent[0]
    mir = ent[1]
    print(query,mir)
    
    ###To write result for every miRNA to separate file
    filename = ent[0].split()[0]
#    print(filename)
    fh_out=open('./target_finder_results/%s_targets' % (filename), 'w')
    
    pipe = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/TargetFinder_1.6/targetfinder.pl", "-c", cutoff, "-s", mir , "-d", cDNA_file, "-q", query], stdout=subprocess.PIPE, universal_newlines=True)
    
    ##This will work as working on terminal
    #pipe = subprocess.call(["/home/setu/Analysis/TargetFinder_1.6/targetfinder.pl", "-s", mirna , "-d", cDNA_file, "-q", query])
    
    result = pipe.stdout.read()
    
    print (result)
    
    fh_out.write('%s' % (result))
    fh_out.close()





fh_in.close()




'''
Created on Apr 25, 2012

@author: atul
'''

