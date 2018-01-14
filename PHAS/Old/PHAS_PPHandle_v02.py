#!/usr/local/bin/python3

import os
import sys
import mysql.connector as sql
import subprocess
from multiprocessing import Process, Queue, Pool
import multiprocessing

#############

server = 'raichu.dbi.udel.edu'
db = 'RICE_drt_PARE'
deg = 'Y'
PARE = 'RICE_drt_PARE_RUN_MASTER_744_1195_1196_1199_1200_1827_1831'
geno_index  = '/alldata/Genomic/Rice/MSU7/BowtieGenomicIndexes/RICE_MSU7_genome'
nproc ='Y'
geno = 'Y' ### Run on whole genome of transcriptome file
GetLib = 'N' ### If you want to run on specific libs than 'N' and specifiy libs below

#############

def ConnectToDB(server):
    
    ##infile values are '0' when you dont want to pulaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    
    print ('\nTrying to connect to mySQL server on %s' % (server))
    # Try to connect to the database
    try:
        con=sql.connect(host= server, user='kakrana', passwd='livetheday')###local_infile = 1 not supported yet so a table has to be updated on row basis
        print ('Connection Established\n')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con

def GetLibs(con,db):
    ##Function that gets just the list names, required to run script in parts.
    cur = con.cursor()
    cur.execute('select distinct(lib_id) from %s.library' % (db))
    libs = cur.fetchall()
    #print (libs)
    print ('\nTotal number of sRNA libraries found: %s\n' % (len(libs)))
    
    return libs###

def PHASBatch(con,libs,geno,geno_index,deg):
    
    #os.mkdir('./%s' % (lib))
    #output_path = './%s' % (lib)
    
    for lib in libs:
        print (lib)
        cur = con.cursor()
        cur.execute('SELECT processed_path FROM master.library_info where lib_id = %s' % (lib))
        path = cur.fetchall()
        #print(path[0][0])
        
        pro_file = path[0][0].replace('$ALLDATA', '/alldata')###Processed sRNA file
        out_file = '%s.txt' % (lib)
        rl = str(21)
        nproc2 = str(nproc)
        sRNAratio = str(75)
        print (pro_file)
        
        if geno == 'Y':###Uses Whole genome as input
            if deg == 'Y':
                retcode = subprocess.call(["/usr/local/bin/perl", "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call(["/usr/local/bin/perl", "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file,"-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
        
        else: ### Uses FASTA file of genes as input         
            #pipe =subprocess.Popen(["/usr/local/bin/perl5.18", "-v"])
            if deg == 'Y':
                retcode = subprocess.call(["/usr/local/bin/perl", "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call(["/usr/local/bin/perl", "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
                    
        
        if retcode == 0:
            pass
        else:
            print("Problem with Phasing script - Return code not 0")
            sys.exit()
        
    return lib


def main(server):
    con = ConnectToDB(server)
    if GetLib ==  'Y':
        libs = GetLibs(con,db)
    else:
        #libs=[(1627,),(2317,),(1628,),(2318,),(1629,),(2319,),(1630,),(2320,)]###Specify libs here
        libs=[(759,),(760,),(762,),(1782,),(1783,),(1784,)]
        
    PHASBatch(con,libs,geno,geno_index,deg)

 



if __name__ == '__main__':

    ###Processors to use####
    if nproc == 'Y':###Use default 70% of processors
        nproc = int(multiprocessing.cpu_count()*0.9)
    else:##As mannually entered by the user
        nproc == int(nproc)
    ###############
    
    main(server)
    sys.exit()
    
    
###Version 01 -> v02
###Added PARE switch
###Added sRNA ratio option
###Added option to specify libs

