#!/usr/local/bin/python3

###The script grabs mature.fa.gz and organisms.gz of current release from mirbase FTP, parses both files, map phylum as per organism code in mirname and populates mirBASE table in mir_master
###Usage ./GetmiRBASE_1a.py
###Written and maintained by kakrana@udel.edu

from ftplib import FTP
import sys
import mysql.connector as sql
import gzip

#################################Config#####################################

##Server to upload data a.k.a destination server with table to hold results of PARE validation

##**Check Phylum index -miRBASE chnages it quite often**

destserver = 'raichu.dbi.udel.edu'#####Value should be in ''
ftpurl = 'mirbase.org'
ver = '19'
TableWipe = 'N' ## Wipe table before update 'Y' and 'N'

#############################################################################

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

def FTPGet(ftpurl):
    
    ftp = FTP(ftpurl)
    ftp.login() ###Anonymous login, if required
    #ftp.retrlines('LIST') ### List root folder
    #ftp.cwd("/pub/mirbase/CURRENT")## Get to the current version
    ftp.cwd("/pub/mirbase/%s" % (ver))## Get to the current version
    #ftp.retrlines('LIST')
    
    #gzfile = open( 'mature.fa.gz','wb') ### wb is for write binary
    #ftp.retrbinary('RETR mature.fa.gz', gzfile.write)
    files = ['mature.fa.gz','organisms.txt.gz']
    
    
    for file in files:
        
        ##miR_file = 'mature.fa.gz'
        #ftp.retrbinary('RETR mature.fa.gz', open(mir_file,'wb').write)
        ftp.retrbinary('RETR %s' % file, open(file,'wb').write)
        
    ##Lets unpack the file
    ftp.quit()
    return files ###Both files in form of list
    
def mirBASEParse(files):
    
    ##Read and Parse the organism file
    org_in = gzip.open(files[1], 'rb')
    org_raw = org_in.read()
    org_in.close()
    org_code = org_raw.decode('utf8')
    #print (org_code)
    
    ###aqu	AQU	Amphimedon queenslandica	Metazoa;Porifera;
    code_dict = {} ###dictionary to hold codes and phylum
    org_blocks = org_code.split('\n')
    for i in org_blocks[:-1]:##exclude item after last newline 
        ent = i.split()
        #print (ent)
        org_code = ent[0]
        phylum = ent[-1].split(';')[0]##miRBASE 18/19 had phylum as last column
        #phylum = ent[-2].split(';')[0]##miRBASE 20 had phlum column shifted to second last from last
        code_dict[org_code] = phylum        
        #print (phylum)
        #break
    
    mir_list = []##list that will hold result
    #print (afile)
    fh_in = gzip.open(files[0],'rb')
    mir_base_raw = fh_in.read()
    fh_in.close()
    mir_base = mir_base_raw.decode('utf8')
    #print(mir_base)
    mir_blocks= mir_base.split('>')
    
    ##>cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
    ##UGAGGUAGUAGGUUGUAUAGUU
    
    for i in mir_blocks[1:]:
        block = i.strip('\n')##Remove newline from end of entry
        ent = block.split('\n')##Use the newline between header and sequence to split
        info = ent[0].split()
        #print ('This is one entry:', info)
        name = info[0]##miRNA name
        org = name.split('-')[0]##Organism code
        phylum = code_dict[org]
        accession = info[1]
        organism = '%s %s' % (info[2],info[3])
        gen_name = info[4]
        mir_seq = ent[1]
        mir_len = len(mir_seq)
        print (name,accession,organism,gen_name,mir_seq,mir_len, phylum)
        
        mir_list.append((name,accession,organism,gen_name,mir_seq,mir_len, phylum))
    #    
    #for i in mir_list:
    #        print(i)
    
    return mir_list
        
def TableUpload(con, mir_list):###Con is connection and res_upload is file from last module which has genomic coords and need to be upload
    ##Upload scoring_input_extend to table - please make sure that columns are in file are in order of the columns in the table
    
    cur = con.cursor()##Connect to destination server with table
    #res_file = 'scoring_input_extend_upload'   
    
    
    if TableWipe == 'Y':
        print ('\nClearing the table before updating.....')
        cur.execute("TRUNCATE TABLE mir_master.mirBASE" % (dest[0],dest[1]))## Clear table before writing
        #con2.commit()
        print ('\nTable cleared successfully, update in process')
        
    else:
        pass
        
    
    ##Current implementation of mysql.connector does not support instant upload by local file - see ConnectToDB() module for implementation
    ##Original query - LOAD DATA LOCAL INFILE "./scoring_input_extend" INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ',';
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    ##So fill table on row by row basis
    
    add_row = 'INSERT INTO mir_master.mirBASE (mir_name,accession,organism,generic_name,mir_seq,mir_len,phylum,version) VALUES (%s, %s, %s, %s, %s, %s, %s, %s)'
    
    #test_data =cel-let-7-5p MIMAT0000001 Caenorhabditis elegans UGAGGUAGUAGGUUGUAUAGUU 22
    
    print ('\n**Updating Table**\n')
    for ent in mir_list:
        #print(ent[0], ent[1], ent[2], ent[3], ent[4], ent[5])
        
        res_upload = (ent[0], ent[1], ent[2], ent[3], ent[4], ent[5], ent[6],ver)###Index in db and table is just to remove '' marks
        cur.execute(add_row,res_upload)        
        con.commit()
    
    print('**Table Updated sucessfully')
        
    cur.close()
    
    
    
def main():    
    gzfile = FTPGet(ftpurl)
    mir_list = mirBASEParse(gzfile)
    con=ConnectToDB(destserver, 0)
    TableUpload(con, mir_list)
    
    


if __name__ == '__main__':
    #url = 'ftp://mirbase.org/pub/mirbase/CURRENT/'
    afile = './mature.fa'
    main()
    sys.exit()
    
    