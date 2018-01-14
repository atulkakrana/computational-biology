#!/usr/local/bin/python3

import sys
import mysql.connector as sql

### This script accepts coordinates and return tags of requested size

coordFile   = "Asparagus_het-siRNA_pval0.05_Dcall24.txt"
header      = 1                                 ## 0: NO | 1: Yes
name        = 1
chrPos      = 3
startPos    = 4
endPos      = 5
tagLen      = [24]

DB          = "ASPARAGUS_privPHAS1_sRNA"
hits        = 100                                ## Hits cutoff
abun        = 10
sep         = '\t'

server      = 'raichu.dbi.udel.edu'


def getTags(coordFile,con):

    '''Get tags that pass hit and abundace cutoff'''
    cur = con.cursor()
    fh_in = open(coordFile,'r')
    if header == 1:
        fh_in.readline()
    coordRead = fh_in.readlines()
    fh_in.close()

    fh_out = open("%s.tags.fa" % coordFile.rpartition('.')[0],'w')

    mainSizeQuery = '' ## Empty Query
    for i in tagLen:
        if mainSizeQuery == '':
            mainSizeQuery = 'len = %s' % (i)
        else:
            mainSizeQuery = mainSizeQuery + ' or len = %s' % (i)

    for i in coordRead:
        ent = i.strip("\n").split(sep)
        aname   = ent[name-1]
        achr    = ent[chrPos-1] ## Convert to Python format
        astart  = ent[startPos-1] ## Convert to Python format
        aend    = ent[endPos-1] ## Convert to Python format
        print("Entry: Name:%s | achr:%s | astart:%s | aend:%s" % (aname,achr,astart,aend))
        # print("This is size query:%s" % (mainSizeQuery))

        cur.execute("select tag,position,hits,norm_sum from %s.tag_position where chr_id = %s and (position between %s and %s) and norm_sum >= %s and hits <= %s" % (DB,achr,str(astart),str(aend),str(abun),str(hits)))####
        temp = cur.fetchall()
        # print(temp)

        for res in temp:
            tag  = res[0]
            if len(tag) in tagLen:
                tagName = "%s_%s_%s_%s_%s" % (name,achr,res[1],res[2],res[3])
                fh_out.write("%s\n%s\n" % (tagName,tag))
            else:
                # print("Tag len %s not required" % (len(tag)))
                pass


    fh_out.close()

    return None


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

def main():
    con = ConnectToDB(server)
    getTags(coordFile,con)

if __name__ == '__main__':
    main()
    sys.exit()





