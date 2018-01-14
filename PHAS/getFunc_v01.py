#!/usr/local/bin/python3

import mysql.connector as sql
import re
import sys
import os


#####Input########

csvfile = 'sco_inp_ext_geno'
sep = ',' ## \t for tab , for comma
pos = '2' ## In excel format and not programming format - Gene name

##3B. Genome database name
DBtype = 'A' ### If Genome DB than 'G' if personal annotation as in kakrana.Rice GFF than 'A'
GenomeDB = 'MAIZE_AGPv2_genome'
GFFDB = 'kakrana'

##9. Server with PARE data
dataserver = 'erawan.dbi.udel.edu'####Value should be in ''




#####Functions#####################################################################################################################

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

def GetFunc(con,csvfile):
    
    print('\nGetting functions\n')
    
    ##Declare output file
    funcfile = '%s_info.tsv' % (csvfile)
    fh_out = open(funcfile,'w')
    fh_out.write('Gene\tModel\tTitle\n')
    
    cur= con.cursor()
    fh_in = open(csvfile,'r')
    entries = 0
    found = 0
    for ent in fh_in:
        
        ent_splt = ent.strip('\n').split(sep)
        #print (ent_splt)
        gene = ent_splt[int(pos)-1]##Position converted to programming format
        #print(gene)
        cur.execute("SELECT gene,model_cnt,title FROM %s.gene_master where gene = '%s'" % (GenomeDB,gene))
        info = cur.fetchall()
        if info:
            print (info[0][0])
            fh_out.write('%s\t%s\t%s\n' % (info[0][0],info[0][1],info[0][2]))
            found +=1
        else: ## Not found
            print ('Title not found for: %s' % (gene))
            pass
        
        entries += 1
    print ('Total entries: %s | Found entries: %s' % (found,entries))
    fh_in.close()
    fh_out.close()
    
    return funcfile

def main():
    
    con = ConnectToDB(dataserver,0)
    funcfile = GetFunc(con,csvfile)
    
########################################### MAIN #####################################################################

if __name__ == '__main__':
    main()
    sys.exit()
    
    
##v01 -> Oct10, 2013

