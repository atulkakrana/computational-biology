#!/usr/local/bin/python3

###Coordinates used to fetch abundance of miRNAs from different library of tag summary table


import sys
import mysql.connector as sql


########Settings#######
server = 'erawan.dbi.udel.edu'
db = 'kakrana'
tablename = 'QIFA_tag_position_summary'
#coords = 'V2_Abun_Clustered_Union_Phase_21_score_30.csv'###Coord format should be comma separated
#coords = str(sys.argv[1])###File name
coords = 'LDMAR.txt'

sep = '.'
header = 'N'## Does first row of file contains header

useStrand = 'N'### Not used for phased loci as sRNA could be from either strand
strandpos = '2'

chrpos = '2' ##excel format
startpos =  '3' #### coordinate strat- excel format
endpos = '4' ## coordinate end - excel format
alen = [21,24] ### Length of sRNAs to calculate abundance for
startbuff = int(25)### This will be used to reduce are back to actual Loci
endbuff = int(0)
###########################

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


def Coords2Abun(con,coords,db):
    cur = con.cursor()
    
    file_out = coords.split('.')[0]+'Abun.csv'
    fh_out = open(file_out,'w')

    ###get column anmes
    columns = []## empty list
    cur.execute("describe kakrana.QIFA_tag_position_summary")
    tablefield = cur.fetchall()
    print (tablefield)
    for i in tablefield:
        col = i[0]
        columns.append(col)

    #lib_col = columns[8:]
    lib_col =",".join(str(x) for x in columns[8:])### Manually mentioned strating lib column - should work on all tag position summary tables
    print(lib_col)
    ##Write file header
    fh_out.write('Loci,%s,%s\n' % (lib_col,lib_col))
    
    fh_in = open(coords,'r')
    if header == 'Y':
        fh_in.readline()

    
    for i in fh_in:
        cur = con.cursor()
        coord = i.strip('\n').split(sep)
        print ('\nRow',coord)
        start = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
        start2 = int(start)+startbuff ### Reverse the effect on co-ordinate made during fetching fasta
        end = str(coord[int(endpos)-1]).strip()
        end2 = int(end)-endbuff
        length = int(end2)-int(start2) ###Extra 1 to include last position, buffer of 100 bp from start and end
        chr_id = str(coord[int(chrpos)-1]).strip()
        strand = str(coord[int(strandpos)-1]).strip()
        print ('This is chr_id:%s and Start:%s and End:%s | Length of loci: %s bp ' % (chr_id,str(start2),int(end2),length))
        cur = con.cursor()
        ###get alen[0] i.e 21nt
        if useStrand == 'Y':
            #cur.execute("select sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589),norm_sum,max_norm from %s.%s where chr_id = %s and strand = '%s' and len = %s and (position between %s and %s)" % (db,tablename,chr_id, strand, alen[0],str(start2),str(end2)))####          
            queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(") ## String manipulation to automativcally convert library names into query format below - used in all queries
            cur.execute("select %s,norm_sum,max_norm from %s.%s where chr_id = %s and strand = '%s' and len = %s and (position between %s and %s)" % (queryLibs,db,tablename,chr_id, strand, alen[0],str(start2),str(end2)))####
            temp = cur.fetchall()
            
            ###Get alen[1] i.e 24,22 or whatever mentioned in settings
            #cur.execute("select sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589),norm_sum,max_norm from %s.%s where chr_id = %s and strand = '%s' and len = %s and (position between %s and %s)" % (db,tablename,chr_id, strand, alen[1],str(start2),str(end2)))####
            cur.execute("select %s,norm_sum,max_norm from %s.%s where chr_id = %s and strand = '%s' and len = %s and (position between %s and %s)" % (queryLibs,db,tablename,chr_id, strand, alen[1],str(start2),str(end2)))####
            temp2 = cur.fetchall()
            lociname = '%s_%s_%s_%s' % (chr_id,strand,start,end)
        else:
            queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
            #cur.execute("select sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589),norm_sum,max_norm from %s.%s where chr_id = %s and len = %s and (position between %s and %s)" % (db,tablename,chr_id, alen[0],str(start2),str(end2)))####
            cur.execute("select %s,norm_sum,max_norm from %s.%s where chr_id = %s and len = %s and (position between %s and %s)" % (queryLibs,db,tablename,chr_id, alen[0],str(start2),str(end2)))####
            temp = cur.fetchall()
            ###Get alen[1] i.e 24,22 or whatever mentioned in settings
            #cur.execute("select sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589),norm_sum,max_norm from %s.%s where chr_id = %s and len = %s and (position between %s and %s)" % (db,tablename,chr_id, alen[1],str(start2),str(end2)))####
            cur.execute("select %s,norm_sum,max_norm from %s.%s where chr_id = %s and len = %s and (position between %s and %s)" % (queryLibs,db,tablename,chr_id, alen[1],str(start2),str(end2)))####
            temp2 = cur.fetchall()
            lociname = '%s_%s_%s' % (chr_id,start,end)
        
        #print (len(temp[0]))
        abundances = ",".join(str(x) for x in temp[0]) ## in comma separated format
        abundances2 = ",".join(str(x) for x in temp2[0]) ## for different size  i.e 22,24 or what ever defined in settings
        fh_out.write('%s,%s,%s\n' % (lociname,abundances,abundances2))
        #print (temp)
    fh_out.close()
    fh_in.close()
    
    return file_out
        
        
    
def main(db,coords):
    con = ConnectToDB(server)
    Coords2Abun(con,coords,db)


if __name__ == '__main__':
    main(db,coords)
    sys.exit()

    
def main(db,coords):
    con = ConnectToDB(server)
    Coords2Abun(con,coords,db)


if __name__ == '__main__':
    main(db,coords)
    print('\n The script finished succesfully')
    sys.exit()
    
    
    
####sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589)
### v01-> v02
## use string manipulation to covert column names used in header to convert to mysql query for ex.. column header = lib_2435,lib_2436,lib_2437 convert to sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438)