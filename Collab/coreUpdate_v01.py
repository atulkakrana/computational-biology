#!/usr/local/bin/python3
##Script written by karana@udel.edu
##Script for updation of core tables i.e. isyte, sysCLFT tables

################# LIBRARIES #######################

import mysql.connector as sql
import re
import sys
import os
import subprocess
import csv 
import glob
import time

###################################################
#
#
#
#### CONFIG ######################################

#1A. Input microarray filename
coreFile = 'AF_RESULTS_GLOB_pval_adjust_v2.txt'##The FASTA header will be cleaned automatically

#1B. Number of samples - Will be counted automatically from file and checked from this value to ensure safe update
nSamples = 10

#1C. delimiter in file - Usually tab delimited text file
delim = '\t'

##2A. Server address
server = 'tarkan.dbi.udel.edu'

##2B. Genome database name
DB = 'kakrana2'
table = 'isyte2'

##3.  Drop the Results table before uploading results - Options - 'Y' and 'N' -- Be careful selecting 'Y' will wipe out the table and you will lose all the data
TableUpdate = 'Y'
TableWipe = 'Y'

##################################################
#
#
#
#### FUNCTIONS ###################################

##Module to connect to DB
def ConnectToDB(server, infile):
    
    ##infile values are '0' when you dont want to pulaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    
    print ("\n** Trying to connect to mySQL server on '%s' **" % (server))
    # Try to connect to the database
    try:
        con=sql.connect(host= server, user='kakrana', passwd='livetheday')###local_infile = 1 not supported yet so a table has to be updated on row basis
        print ("** Connection Established **\n")

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con

## Function to make a long format list of core table input file
## Column format - ID, symbol, FC cols, coef cols, t-stat cols, p-val cols,exp cols
def processCoreTab(coreFile):
    
    coreList = [] ## List to hold all entries - Will be passed to table update function in end
    infoCol = 2 ## Python format - Column from where to start a frame for extracting all information from different time points
    frame = 0 ## Increment frame after every timepoint
    
    fh_in = open(coreFile,'r')
    header = fh_in.readline().split(delim)## Get header to extract stage - EM and PN ; and timepoint i.e 10.5, 11.5 etc
    coreRead = fh_in.readlines()
    
    ## For all the stages
    for sample in range(nSamples):
        
        tempList = []## To record entries for every sample/stage i.e E10.5 - Will be appnded to core List in end
        head = header[infoCol+frame].split()[0].split('FC.') ## Example header 'FC.E10.5 - WB'
        #print(head)
        print('\nStage/Sample being processed: %s\n' % (head[1]))
        stage = head[1][0]
        timepoint = head[1][1:]
        #print(stage,timepoint)
        
        ## For entry at every stage
        for entry in coreRead:
            ent = entry.strip('\n').split(delim)
            #print('\nEntry:',entry)
            #print('Info:',stage,timepoint)
            #print (ent[0],ent[1],ent[infoCol+frame],ent[infoCol+nSamples+frame],ent[infoCol+2*nSamples+frame],ent[infoCol+3*nSamples+frame],ent[infoCol+4*nSamples+frame],ent[-1],stage,timepoint,ent[infoCol+5*nSamples+frame])
            tempList.append((ent[0],ent[1],ent[infoCol+frame],ent[infoCol+nSamples+frame],ent[infoCol+2*nSamples+frame],ent[infoCol+3*nSamples+frame],ent[infoCol+4*nSamples+frame],ent[-1],stage,timepoint,ent[infoCol+5*nSamples+frame]))## Tuple is appended -id,symbol,FC,coef,t-stat.p-val,gene expression, WB expression, stage, time point and tissue
       
        ## After one stage/sample - Increment sample column through frame variable as nsample is used as distance parameter for columns i.e distance between columns pertaining to different stage
        coreList.extend(tempList) ## Extend used and not append because 'append' will add the tempList as one item not as list
        frame += 1
        
    fh_in.close()
    print('This is the length of Corelist',len(coreList))
    return coreList

####Module that uploads data to already made table in a DB
def TableUpload(con,coreList,DB,table):###Con is connection and res_upload is file from last module which has genomic coords and need to be upload
    ##Upload scoring_input_extend to table - please make sure that columns are in file are in order of the columns in the table
    
    cur = con.cursor()##Connect to destination server with table
   
    if TableWipe is 'Y': ## If the Global setting is 'Y' than table will be dropped - Be careful   
        print ('\n**Clearing the table before updating.....**')
        cur.execute("TRUNCATE TABLE %s.%s" % (DB,table))## Clear table before writing
        #con2.commit()
        print ('\n**Table cleared successfully, update in process**\n')
        
    ##Current implementation of mysql.connector does not support instant upload by local file - see ConnectToDB() module for implementation
    ##Original query - LOAD DATA LOCAL INFILE "./scoring_input_extend" INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ',';
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    ##So fill table on row by row basis
    ##Check for miRNA table variable if not than Local file used else - public/priv
    print('**Updating table: %s.%s**' %(DB,table))
    add_row = "INSERT INTO kakrana2.isyte2 (ID,symbol,fc,coef,tstat,pval,avexp,wbexp,stage,time,tissue) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
    entcount = 0 ## Count the entries
    for ent in coreList:
        #print(ent)
        #print(ent[0], ent[1], ent[2], ent[3], ent[4], ent[5], ent[6], ent[7],ent[8],ent[9],ent[10])
        res_upload = (ent[0], ent[1], ent[2], ent[3], ent[4], ent[5], ent[6], ent[7],ent[8],ent[9],ent[10])
        cur.execute(add_row,res_upload)        
        con.commit()
        entcount+=1 ## Incremented with every entry- should be equal to entry in input file * nSample
    
    ### Test if uploaded entry count is correct
    with open(coreFile) as f:
        f.readline()##Waste header
        flength = len(f.readlines())
        print('File Length: %s' % (flength))
    
    if int(entcount) == int(flength)*nSamples:
        print('**Entries uploaded matches with total entries in file multiplied by number of samples - Whoooo!!**')## If entries uploaded  is equal to numbe ro fentries in main file minus header * number of samples than job done
    else:
        print('^^^Uploaded entry count does not matches with number of entries in input file * nSamples - Something is fishy - Exiting^^^')
        print('Entries updated: %s AND entries from file: %s' % (entcount,flength))
        sys.exit()
    cur.close()

##################################################
#
#
#
#### MAIN ########################################

def main():
    con = ConnectToDB(server,0) #####The second input '0' is for future use in case of local_infile upload to update table
    coreList = processCoreTab(coreFile)
    
    if TableUpdate == 'Y':
        TableUpload(con,coreList,DB,table)
    else:
        print ('* Table not updated *')


##################################################

if __name__ == '__main__':
    main()
    print ("\nScript finished sucessfully\n")
    sys.exit()    

#### LOG ########################################

##v01 started on Dec-2014
## 1. Make check for every stage so that correct stage is read and loaded
## 2. Table name is hard coded in update string - CRITICAL
## 3. Entry check in Tableupload is giving false alarm


#### ASSUMPTIONS ###############################

###FC column is used to find stage and timepoint - FC.E10.5, here single alphabet 'E' will be captured as stage and later as time point