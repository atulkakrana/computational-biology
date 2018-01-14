#!/usr/local/bin/python3

## Gets genomic co-ordinates for tag including all hits and gene name

import sys,os
import mysql.connector as sql

########Settings#######
db = 'ASPARAGUS_privPHAS1_sRNA' ## Tag position summary table required
tagFile = 'Final_miRNAs_v02.txt'
sep = '\t'
header = 'Y'
tagpos = 6 ## Excel format

server = 'raichu.dbi.udel.edu'
######################


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

def tag2Pos(con,db,tagFile):
    
    fh_in = open(tagFile,'r')
    outfile = "%sAllHits.csv" % tagFile
    fh_out = open(outfile,'w')
    fh_out.write("tag,len,chr_id,strand,start,end,hits,gene\n")
	
    if header == 'Y':
        fh_in.readline()

    fileRead = fh_in.readlines()

    for i in fileRead:
        ent = i.split(sep)
        print("\n",ent)
        tag = ent[tagpos-1]
        cur = con.cursor()
        cur.execute("select chr_id,strand,position,hits,gene from %s.tag_position where tag = '%s'" % (db,tag))
        temp = cur.fetchall()

        print(tag,temp)

        ## Write results
        for i in temp:
            chr_id = i[0]
            strand = i[1]
            start = i[2]
            end = int(i[2])+len(tag)-1
            hits = i[3]
            gene = i[4]

            fh_out.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (tag,len(tag),chr_id,strand,start,end,hits,gene))

    fh_in.close()
    fh_out.close()
    
    return outfile

def main(db,tagFile):
    print(server)

    con = ConnectToDB(server)
    outfile = tag2Pos(con,db,tagFile)


if __name__ == '__main__':
    main(db,tagFile)
    sys.exit()




