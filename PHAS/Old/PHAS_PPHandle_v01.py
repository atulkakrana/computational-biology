#!/usr/local/bin/python3

import os
import sys
import mysql.connector as sql
import subprocess
from multiprocessing import Process, Queue, Pool
import multiprocessing

#############

server = 'raichu.dbi.udel.edu'
db = 'RICE_sbsQIFA_sRNA'
deg = 'RICE_sbsQIFA_PARE_MASTER'
geno_index  = '/alldata/Genomic/Rice/MSU7/BowtieGenomicIndexes/RICE_MSU7_genome'
nproc ='Y'
geno = 'Y'

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
        print (pro_file)
        
        if geno == 'Y':###Uses Whole genome as input
            retcode = subprocess.call(["/usr/local/bin/perl", "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file, "-q", deg, "-f", "t", "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            
        else: ### Uses FASTA file of genes as input         
            #pipe =subprocess.Popen(["/usr/local/bin/perl5.18", "-v"])
            retcode = subprocess.call(["/usr/local/bin/perl", "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-q", deg, "-f", "t", "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
        
        if retcode == 0:
            pass
        else:
            print("Problem with Phasing script - Return code not 0")
            sys.exit()
        
    return lib


def main(server):
    con = ConnectToDB(server)
    libs = GetLibs(con,db)
    PHASBatch(con,libs,geno,geno_index,deg)

 



if __name__ == '__main__':

    ###Processors to use####
    if nproc == 'Y':###Use default 70% of processors
        nproc = int(multiprocessing.cpu_count()*0.8)
    else:##As mannually entered by the user
        nproc == int(nproc)
    ###############

    main(server)
    sys.exit()

