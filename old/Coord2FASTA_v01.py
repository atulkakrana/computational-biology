#!/usr/local/bin/python3


import sys
import mysql.connector as sql


########Settings#######
server = 'raichu.dbi.udel.edu'
genomedb = 'MAIZE_AGPv2_genome'
#coords = 'V2_Abun_Clustered_Union_Phase_21_score_30.csv'###Coord format should be comma separated
coords = str(sys.argv[1])
buff = int(20)
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
        coord = i.split(',')
        #print (coord)
        start = str(coord[2]).strip() ##Buffer of 100bp
        start2 = int(start)-buff
        end = str(coord[3]).strip()
        length = int(end)-int(start)+ buff*2 ###Extra 1 to include last position, buffer of 100 bp from start and end
        chr_id = str(coord[1]).strip()
        strand = str(coord[3]).strip()
        print ('This is chr_id:%s and Start:%s and End:%s | Length of sequence: %s bp ' % (chr_id,str(start2),int(end)+buff,length))
        cur = con.cursor()
        cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = %s" % (str(start2),length+buff,genomedb,chr_id,strand))####
        temp = cur.fetchall()
        #print (temp)
        #print (temp[1][0][::-1])
        fh_out.write('>%s\n%s\n' % (coord[0]+'_pos',temp[0][0]))
        fh_out.write('>%s\n%s\n' % (coord[0]+'_neg',temp[1][0][::-1]))

    fh_out.close()
    
    
def main(genomedb,coords):
    con = ConnectToDB(server)
    Coords2FASTA(con,coords,genomedb)
    

if __name__ == '__main__':
    main(genomedb,coords)
    sys.exit()



