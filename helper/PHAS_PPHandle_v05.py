#!/usr/local/bin/python3

## This script uses the phasiRNA scripy from SVN So, make sure your svn is updated
## Author: Atul Kakrana kakrana@udel.edu

import os
import sys
import mysql.connector as sql
import subprocess
from multiprocessing import Process, Queue, Pool
import multiprocessing
import os.path

#############

db = 'RICE_sbsQIFA_sRNA' ## sRNA DB
geno = 'N' ### Run on whole genome (Y) or transcriptome file (N) - If (N) then provide index of transcriptome
geno_index = '../index/RiceAssemblyTrans'

deg = 'N' ## Use Degradome validation, IF yes enter PARE db in line below
PARE = 'GuturGu'

fetchLibIDs = 'N' ### (Y): All libs in the DB (N): If you want to run on specific libs than 'N' and specifiy libs below
# userLibs = [(2435,),(2436,),(2437,),(2438,),(2439,),(2440,),(2441,),(2495,),(2496,),(2497,),(2498,),(2499,),(2500,),(2501,)] ## Used only if fetchLib == 'N'
userLibs = [(2438,),(2439,)]

Local = 'N' ## (Y): Get the libs from $ALLDATA with raw reads (N) Get library from srna db with reads filtered on number of hits
phase = 21

## Developer options ####
nproc ='Y'
server = 'raichu.dbi.udel.edu'
perl = "/usr/local/bin/perl_5.18" ## Josh updated the perl on Tarkan and its not ready yet for PHAS script FORK is missing and somemore modules -Check with Pingchuan help
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

## Get lib ids in the sRNA DB
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
        rl = str(phase)
        nproc2 = str(nproc)
        sRNAratio = str(75)
        print (pro_file)
        
        if geno == 'Y':###Uses Whole genome as input
            if deg == 'Y':
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file,"-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
        
        else: ### Uses FASTA file of genes as input         
            #pipe =subprocess.Popen(["perl5.18", "-v"])
            if deg == 'Y':
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
                    
        
        if retcode == 0:
            pass
        else:
            print("Problem with Phasing script - Return code not 0")
            sys.exit()
        
    return lib

### sRNA Libraries are fetched from server
def TagAbundanceFile(con,db,libs):
    
        for alib in libs:##For all the libraries
            
            ## Check if file already exsits in directory - This saves a lot of time downloading the same file
            filePath = '%s' % (alib)
            if os.path.isfile(filePath) == False:
                print ('\nPreparing sRNA reads file for library: %s' % (alib[0]))
                #print (lib[0])
                #print ('Caching tag and count information from server for PARE alib %s' % (alib[0]) )
                cur = con.cursor()
                cur.execute("SELECT tag,norm from %s.run_master where lib_id = %s and hits between 1 and 20" % (db,alib[0]))
                lib_info = cur.fetchall()
                #print('These are the tags:',lib_info[:10])
                
                fh_out = open('%s' % (alib), 'w')##Naming file with lib_ids name
                print ('Library cached, writing abundance file')
                tag_num = 1
                for ent in lib_info:## All the PARE tags in a library
                    #print (ent)
                    fh_out.write('%s\t%s\n' % (ent[0],ent[1]))
                    tag_num += 1
                    
                fh_out.close()
            else:
                print('tag abundance file exists for library: %s' % (alib))
                pass

##PHasing anlysis
def PHASBatch2(con,libs,geno,geno_index,deg):
    for lib in libs:
        pro_file = './%s' % (lib)###Processed sRNA file
        out_file = '%s.txt' % (lib)
        
        rl = str(phase)
        nproc2 = str(nproc)
        sRNAratio = str(75)
        print (pro_file)

        if geno == 'Y':###Uses Whole genome as input
            if deg == 'Y':
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.genome.pl", "-i", pro_file,"-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
        
        else: ### Uses FASTA file of genes as input         
            #pipe =subprocess.Popen(["perl5.18", "-v"])
            if deg == 'Y':
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-q", PARE, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])
            else:
                retcode = subprocess.call([perl, "/data2/homes/kakrana/svn/users/pingchuan/phasiRNA_prediction_pipeline.ver.MUL.pl", "-i", pro_file, "-f", "-t", sRNAratio, "-d", geno_index, "-px", out_file, "-rl", rl, "-cpu", nproc2])

        if retcode == 0:
            pass
        else:
            print("Problem with Phasing script - Return code not 0")
            sys.exit()
        
    return lib

def main(server):
    con = ConnectToDB(server)
    if fetchLibIDs ==  'Y': ## Run on all libraries in the DB
        if Local == 'Y': ## Get sRNA lib data from $ALLDATA
            libs = GetLibs(con,db)
            print('These are the libs: %s' % (libs))
            PHASBatch(con,libs,geno,geno_index,deg)
        else: ## Get sRNA lib data from DB
            libs = GetLibs(con,db)
            print('These are the libs: %s' % (libs))
            TagAbundanceFile(con,db,libs)
            PHASBatch2(con,libs,geno,geno_index,deg)
            
    else: ## Run on specific libraries
        #libs=[(1627,),(2317,),(1628,),(2318,),(1629,),(2319,),(1630,),(2320,)]###Specify libs here
        #userLibs=[(3381,),(3416,),(3417,),(3418,),(3419,),(3420,),(3421,),(3422,),(3423,)]
        if Local == 'Y': 
            print('These are the libs: %s' % (userLibs))
            PHASBatch(con,userLibs,geno,geno_index,deg)
        else:
            #libs = GetLibs(con,db)
            print('These are the libs: %s' % (userLibs))
            TagAbundanceFile(con,db,userLibs)
            PHASBatch2(con,userLibs,geno,geno_index,deg)

if __name__ == '__main__':

    ###Processors to use####
    if nproc == 'Y':###Use default 70% of processors
        nproc = int(multiprocessing.cpu_count()*0.8)
    else:##As mannually entered by the user
        nproc == int(nproc)
    ###############
    
    main(server)
    print ('\n\nPhasing Analysis finished successfully')
    sys.exit()

### Version 01 -> v02
### Added PARE switch
### Added sRNA ratio option
### Added option to specify libs

## v02 -> v03
## Added option to get libs from the server with hits filter
## COrrected bug in main(), repaced libs with userlibs for specific librarues part
## Perl location added as variable

##v03 -> v04
## Changed order of user settings to make them more clear
## Added functionality to check if the abundance file for library already exists in folder - Saves a lot of time

## v04->v05
## Modified the sRNA tags download module with tagg hits between 1 and 20

## Add automatic index resolution
