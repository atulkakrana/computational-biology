#!/usr/local/bin/python3
#### Script to get sequence from genome DB - written by kakrana@udel.com

import sys
import mysql.connector as sql

######## SEETINGS ################
server = 'raichu.dbi.udel.edu'
genomedb = 'RICE_MSU7_genome'
#coords = 'V2_Abun_Clustered_Union_Phase_21_score_30.csv'###Coord format should be comma separated
coords = str(sys.argv[1])

SpecStrand = 'Y' ### If strand is available and sequence needs to be from specific strand, If 'N' than you dont care about strand and want from both
strandpos = '6' ##Column for strand - excel format - should be 'w' and 'c'

sep = ',' ### Specify seprator used in spreadsheed or input file
chrpos = '3' ##Column for chromosome - excel format
startpos =  '4' #### Column for start postion- excel format
endpos = '5' ## Column for end position - excel format

buff = int(50)
##################################

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

def Coords2FASTA(con,coords,genomedb):
    cur = con.cursor()
    
    file_out = coords.split('.')[0]+'.fa'
    fh_out = open(file_out,'w')
    
    fh_in = open(coords,'r')
    fh_in.readline()
    
    for i in fh_in:
        
        cur = con.cursor()
        coord = i.split(sep)
        #print (coord)
        start = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
        start2 = int(start)-buff
        end = str(coord[int(endpos)-1]).strip()
        length = int(end)-int(start)+ buff*2 ###Extra 1 to include last position, buffer of 100 bp from start and end
        chr_id = str(coord[int(chrpos)-1]).strip()
        strand = str(coord[int(strandpos)-1]).strip()
        #print ('Strand', strand)
        print ('This is chr_id:%s and Strand: %s and Start:%s and End:%s | Length of sequence: %s bp ' % (chr_id,strand,str(start2),int(end)+buff,length))
        cur = con.cursor()
        
        #print (temp[1][0][::-1])
        if SpecStrand == 'Y': ##Strand info available - Like if coordinates are for gene
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start2),length+buff,genomedb,chr_id,strand))####
            temp = cur.fetchall()
            #print (temp)
            #print (temp[1][0][::-1])
            if strand == 'w': ## Distinguish between strands so that you can reverce 'crick' strand to give correct orientation
                fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end,temp[0][0]))
            else:
                fh_out.write('>%s_c_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end,temp[0][0][::-1])) ### Just reverse to give correct orientation

 
        else: ###No strand information both strands were required   
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s" % (str(start2),length+buff,genomedb,chr_id))####
            temp = cur.fetchall()
            #print (temp)### Both positive and negative strand reported
            fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end,temp[0][0]))
            fh_out.write('>%s_c_%s:%s:%s\n%s\n' % (coord[1],chr_id,start2,end,temp[1][0][::-1]))### Reversed to give correct orientation

    fh_out.close()
    
    
def main(genomedb,coords):
    con = ConnectToDB(server)
    Coords2FASTA(con,coords,genomedb)
    

if __name__ == '__main__':
    main(genomedb,coords)
    sys.exit()




####V01 -> v02 11-Oct-13
###Inputs added for strand, seprator, and other column position sin input file
### If coords has strand info than choose 'SpecStrand' mode and specify 'strandpos' variable i.e position in excel or soreadsheet
### If coords info not available than 'SpecStrand' should be 'N' in this cse sequence from both strands reported back
