#!/usr/local/bin/python3

###A new miRNA target prediction with fast mapping, multiprocessing and exhaustive algorithm for whole genome target prediction
###Possible tools - shrimp/segemehl, OPENMPI

import os
import sys
import subprocess
import time
import glob
import shutil
import pyfasta
import mysql.connector as sql
import itertools as it
import operator
import multiprocessing
from multiprocessing import Process, Queue, Pool
from operator import itemgetter


afile = 'REAL_targets'

def tarParse(afile):
    
    print ('File for parsing:%s\n' % (afile))
    fh_in = open(afile,'r')
    #fh_in.readline()##Waste header
    
    ##
    
    #print(afile)
    #out_file = './scoring/%s.tar_parsed' % (afile.split('/')[2].rpartition('.')[0])
    out_file = 'parsed'
    fh_out = open(out_file,'w')
    
    tarlist = fh_in.read().split('\n\n')
    mirna_set = set()
    for i in tarlist:
        #print(i,'\n')
        blocks = i.split('\n')
        #print(blocks)
        head = blocks[0].split('\t')
        query = head[0][1:]##remove '>'
        if query not in mirna_set:
            mirna_set.add(query)
        tar = head[2]
        score = head[4]
        strand = head[5]### Single orientation FASTA file used by CL# script so strand always 1
        #strand =1
        info_block = blocks[2].split('\t')
        bind_start = info_block[1].split()[2]
        bind_end = info_block[1].split()[4]
    
        tar_seq = blocks[2].split('\t')[0]
        alingment = blocks[3].split('\t')[0].translate(str.maketrans("|o",":."))
        mir_seq = blocks[4].split('\t')[0]
        ##print(tar_seq,mir_seq)
        #
        #print(alingment)
        #print("query=%s, target=%s, score=%s, range=%s-%s, strand=%s\n\ntarget 5' %s 3'\n          %s\nquery  3' %s 5'\n\n" % (query,tar,score,bind_start,bind_end,strand,tar_seq,alingment,mir_seq))
        
        
        if strand == '+': ## Feray predicts both strands
            newstrand = 1
            #print ('positive strand')
            fh_out.write("query=%s, target=%s, score=%s, range=%s-%s, strand=%s\n\ntarget 5' %s 3'\n          %s\nquery  3' %s 5'\n\n" % (query,tar,score,bind_start,bind_end,newstrand,tar_seq,alingment,mir_seq))
    
    fh_in.close()
    fh_out.close()
    
    return out_file,mirna_set

def phaseOut(mirna_set,out_file):
    
    fh_in = open(out_file,'r')
    #list = fh_in.read().split('\n\n')
    for i in mirna_set:
        #print(i)
        mirna = '%s' % (i)
        fh_out = open('./temp/%s_targets' % (i), 'w')
        pipe = subprocess.Popen(["grep", "-A", "5", mirna, out_file],stdout=subprocess.PIPE, universal_newlines=True)
        results = pipe.stdout.read()
        #print(results)
        fh_out.write('%s' % (results))
        #for line in results.readlines():
        #    print (line)
        ##    line_strpd = line.strip()
        ##    if line_strpd != '--':
        ##        fh_out.write('%s' % (line))
        fh_out.close()
        
        ###Read and remoce '--' due to grep
        fh_in2 = open('./temp/%s_targets' % (i),'r')
        fh_out2=open('./target_finder_results/%s_targets' % (i), 'w')
        alllines = fh_in2.readlines()
        for i in alllines:
            if i.strip('\n') != '--':
                #print (i)
                fh_out2.write(i)
        
        fh_in2.close()
        fh_out2.close()
        #break

    return mirna_set
  
def main():
    out_file,mirna_set = tarParse(afile)
    
    
    
    shutil.rmtree('./temp', ignore_errors=True)
    os.mkdir('./temp')
    shutil.rmtree('./target_finder_results', ignore_errors=True)
    os.mkdir('./target_finder_results')
    phaseOut(mirna_set,out_file)
    
    ##Remoce temp directory
    shutil.rmtree('./temp', ignore_errors=True)


if __name__ == '__main__':
    
    ####Processors to use####
    #if nproc == 'Y':###Use default 80% of processors
    #    nproc = int(multiprocessing.cpu_count()*0.8)
    #else:##As mannually entered by the user
    #    nproc == int(nproc)
    #    ###############
    
    ##Get time
    start = time.time()###time start
    main()
    end = time.time()
    print ('The script run time is %s' % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting........\n')
    ###Echo time