#!/usr/local/bin/python3

###This makes use of targetfinder program from carringtons lab, runs it in batch mode.
###Cleans the header of both input files and predicts target
###Written by Atul Kakrana: kakrana@udel.edu
###Usage: ./TarFindPP_v2.5.py
####This is parallel processing version

import sys
import re
import os
import subprocess
import shutil
import time
import multiprocessing
from multiprocessing import Process, Queue, Pool
#import numpy

################################# CONFIG ####################################################


#seq_file = 'Tow_Inter.csv.fa' ###cDNA/Transcript/Chromosome file in FASTA format
seq_file = str(sys.argv[1])
mirna_file = 'miRinput.fa' ### miRNA file in FASTA format
#mirna_file = str(sys.argv[2])
cutoff = '7'### Score cutoff 1-7, usually a cut-off of 7 is used in our lab. Please be informed that lesser the score the better it is as its a penalty score
rev = '0' ### Predict targets on reverse strand. 0 = NO 1 =Yes. Use this when you are prediction targets from whole genome i.e chromosome as Seq_file
nproc = 'Y'###Number of processors to use

miRTable = 'mirBASE'
ver = '20'
miRCode = 'osa'
genome = 'Y' ## Is it a genome file or Contig file - Not in use as of now maybe when input file is too large we can use it

##9. Server with PARE data
dataserver = 'raichu.dbi.udel.edu'####Value should be in ''

################################# CONFIG ####################################################


#########################AA##############################TT#############################UU###########################LL#########################


####Module to clean headers#####
def CleanHeader(filename):
    #read file
    fh_in=open(filename, 'r')
    #write file
    out_file = ('%s_new_head.fa' % (filename))
    fh_out =open(out_file, 'w')
    
    print ('\nProcessing "%s" file to clean FASTA headers' % (filename))
    
    acount = 0 ## count the number of entries
    for i in fh_in:
        if re.match('>', i):
            header = i.split()##change here if separater is not white space
            new_head = header[0].split('|')[0]## Change here for number of fields you want strating from 0
            fh_out.write('%s\n' % new_head)
            acount+=1
    #        print(i)
    #        print(new_head)
        else:
            fh_out.write('%s' % i)
        
    fh_in.close()
    fh_out.close()
    print('The fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file, acount))
    
    return out_file

##Module to connect to DB
def ConnectToDB(server, infile):
    
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

##InputmiRNAs
def miRinput(con):
    miRs = [] ## List to hold miRNAs
    #fh = open(miRNA_file)
    try:## If user specified some local miRNA file in config
        fh = open(miRNA_file)##Test if file exists or not
        
        miRNA_file_clean = CleanHeader(miRNA_file)## Not required when info is downloaded from server in case of miRPage script
        fh_miRNA = open(miRNA_file_clean, 'r')
        mir_base = fh_miRNA.read()
        mir_blocks= mir_base.split('>')
        for i in mir_blocks[1:]:
            #print (i)
            block = i.strip('\n')##Remove newline from end of entry
            ent = block.split('\n')##Use the newline between header and sequence to split
            #print (ent)
            #print ('%s,%s,%s,%s' % (ent[0],'None','None',ent[1]))
            miRs.append((ent[0],'None','None',ent[1]))## None at position 2 and 3 because so that name and seq at same index when data collected from table
        fh_miRNA.close()
        mirTable = 'None'##To be filled in final table
        print ('***miRNAs fetched from the specified file***\n')
        
    except NameError:## NO Local file specified by User - Fetch  from table
        print ('\nFetching miRs from the miRNA table')
        fh_out = open('miRinput.fa', 'w')## File with miRNAs for troupleshooting
        cur = con.cursor()
    
    ######################TEST SECTION - COMMENT OUT WHEN NOT TESTING###########################
        #cur.execute("SELECT mirna_name, mirna_precursor_name, score, mature_sequence FROM master.mirna_filter_result where score > 1 group by mirna_name limit 24")
    
    #####################REAL WORLD - UNCOMMENT WHEN NOT-TESTING#############
        #cur.execute('SELECT mirna_name, mirna_precursor_name, score, mature_sequence FROM master.mirna_filter_result where score > 1 group by mirna_name')##For kevins table
        cur.execute("SELECT mir_name, accession, organism, mir_seq FROM mir_master.mirBASE WHERE mir_name like '%s%%' AND version = %s" % (miRCode,ver))
    #####################Uncomment when testing is done#########################################
    
        miRs = cur.fetchall()
        print ('Total number of unique miRs found in miR table: %s' % (len(miRs)))
        
        ##Extract miRname and miR sequence for PARE validation and output in FASTA format as read by cleaveland
        #fh_out = open('miRinput.fa', 'w')
        
        for miR in miRs:
            mirname = str(miR[0])##Converted to string from unicode, a unicode mirname example: u'meyers-miR544-5p'
            mirseq = str(miR[3]) ##Converted to string from unicode, a unicode mirname example: u'TGAAGATGAAGAAGATGAAGAAGA'
            #print (mirname, mirseq)
            fh_out.write('>%s\n%s\n' % (mirname,mirseq))
        fh_out.close()
        mirTable = miRTable
        print('****miR info cached****\n')
        
    #for i in miRs:
    #    print (i)
        
    return miRs, miRTable##miRs holds the list of miRNA name and query where as miRtable holds flag -table name or local

####Module to predict potential targets using targetfinder in parallel mode
#
def TargetFinder(ent):
##    print (ent)
    query = ent[0]
    mir = ent[1]
    print('\n\nPredicting targets of miR: %s with sequence: %s' % (query,mir))
    
    ###To write result for every miRNA to separate file
    filename = ent[0].split()[0]
    ###Opening file for each miRNA
    fh_out=open('./target_finder_results/%s_targets' % (filename), 'w')
    
    ###Must use clean header file for cDNA
    if rev is '0':
        
        pipe = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/TargetFinder_1.6/targetfinder.pl", "-s", mir , "-d", seq_file_clean, "-q", query, "-c", cutoff], stdout=subprocess.PIPE, universal_newlines=True)
    
    else:
        print ('Targets on reverse strand will also be predicted\n')
        pipe = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/TargetFinder_1.6/targetfinder.pl", "-s", mir , "-d", seq_file_clean, "-q", query, "-c", cutoff, "-r", rev], stdout=subprocess.PIPE, universal_newlines=True)
    
    result = pipe.stdout.read()
    
    print (result)
    fh_out.write('%s' % (result))
    fh_out.close()    

########################################## MAIN ####################################### MAIN ################################################


def main ():
    
    con = ConnectToDB(dataserver,0)###The second input '0' is for future use in case of local_infile upload to update table
    miRs,mirTable = miRinput(con)###Fetch list of miRNAs from either local file if specifed or miRNA table
    
    ###Time for target prediction
    print('''\n\t\t**************************Warning*******************************
          Any existing 'target_finder_results' folder will be emptied
          Press CTRL+C in next 10 seconds to stop the process and save existing folder
          ****************************************************************************''')
    time.sleep(10)
    ##Delete any existing ./target_finder_results directory
    shutil.rmtree('./target_finder_results',ignore_errors=True)
    ##Open a directory to store all target files
    os.mkdir('./target_finder_results')
    ##Find targets of individual miR and save its file
    print ('\nTargetFinder will now run and find targets of every miR in the file')
    print ('****Please be Patient - Target prediction takes a lot of time depending upon number of miRNA involved****\n')
    
    ##Read file
    
    mir_list = []##list that will be feed to pool for parallel processing
    fh_in = open(mirna_file_clean, 'r')
    mir_base = fh_in.read()
    mir_blocks= mir_base.split('>')
    print (mir_blocks)
    
    for i in mir_blocks[1:]:
        block = i.strip('\n')##Remove newline from end of entry
        ent = block.split('\n')##Use the newline between header and sequence to split
        print (ent)
        mir_list.append(ent)
    #print(mir_list)
    
    #result = TF(mir_list)
        
    #print (ent)
    start = time.ctime()###time start
    
    npool = Pool(int(nproc))
    npool.map(TargetFinder, mir_list)###Cannot give multiple parameters at this time
    
    print ('Parallel TargetFinder started at: %s\n Ended at: %s\n' % (start, time.ctime()))
    
    
###############################################################################################################################################  

if __name__ == '__main__':
    seq_file_clean = CleanHeader(seq_file)
    mirna_file_clean = CleanHeader(mirna_file)
    ###Processors to use####
    if nproc == 'Y':###Use default 70% of processors
        nproc = int(multiprocessing.cpu_count()*0.8)
    else:##As mannually entered by the user
        nproc == int(nproc)
    ###############
    main()
    sys.exit()
    
    

###v2.5 -> v3
###Added - Automatic download of miRNAs fromm the mirbase table
    
    
    
    
        
        
        


        
        




   
    
    
