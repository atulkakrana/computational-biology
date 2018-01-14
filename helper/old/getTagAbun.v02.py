#!/usr/local/bin/python3

### Script to compare sRNA abundances in two databases using summary two table

###### IMPORTS ######################
import sys
import mysql.connector as sql
#####################################


##### INPUT #########################
mode            = 1                                                 ##  0: Fetch from tag position summary table 
                                                                    ##  1: Fetch from sRNA DB - Required when sRNA libraries are mapped to other genome and has no hits
tagfile         = 'Streptochaeta.txt'                               ##  Text file with tags
DB              = 'STREPTOCHAETA_priv_sRNA'                         ##  Used with mode:1 - If abundance to be fetched directly from sRNA DB, required if no genome is vaialble, and hit = 0

testTable       = 'STREPTOCHAETA_priv_sRNA_TagPosSummNorm'          ##  Used with mode:0

CompareDBs      = 'N'                                               ##  Used with mode:0 AND ComapareDBs:Y - Fetch abundance from two tables each for one DB for comparision
conTable        = 'STREPTOCHAETA_priv_sRNA_TagPosSummNorm'

db              = 'kakrana'                                         ##  Used with mode:0
server          = 'raichu.dbi.udel.edu'

#####################################

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

def tags2Abun(con,tagfile,testTable,conTable):
    
    cur = con.cursor()
    
    outfile = tagfile.split('.')[0]+'Abun_FC.csv'
    fh_out = open(outfile,'w')

    ### Get column anmes
    columns = []## empty list
    cur.execute("describe kakrana.%s" % (testTable))
    testTablefield = cur.fetchall()
    cur.execute("describe kakrana.%s" % (conTable))
    conTablefield = cur.fetchall()
    
    #print (testTablefield)
    columns = []## empty list
    for i in testTablefield:
        col = i[0]
        columns.append(col)
        
    testTablecol =",".join(str(x) for x in columns[5:])### Manually mentioned starting lib column - should work on all tag summary2 tables
    print(testTablecol)
        
    #print(conTablefield)
    columns = []## empty list - emptied again
    for i in conTablefield:
        col = i[0]
        columns.append(col)
    conTablecol = ",".join(str(x) for x in columns[5:])### Manually mentioned starting lib column - should work on all tag summary2 tables
    print(conTablecol)
    
    ##Write file header
    fh_out.write('tag,NormSumRatio,%s,%s\n' % (testTablecol,conTablecol))
    
    ### Open and read tag file
    fh_in = open(tagfile, 'r')
    
    for i in fh_in:
        tag = i.strip('\n').replace('-','').translate(str.maketrans("U","T"))## Manipulation if tag is from PARE valid results - reversed and U -> T
        print("'%s' being queried:" % (tag))
        ### Query test table
        cur.execute("select * from %s.%s where tag = '%s'" % (db,testTable,tag))####
        temp = cur.fetchall()
        print ("\nTemp:",temp)
        
        if temp: ## Tag exists in normalized table
            testNormSum = temp[0][-2] ## There could be multiple entries for different hits
            print("NormSum:",testNormSum)
            testAbundances = ",".join(str(x) for x in temp[0][5:]) ## in comma separated format, from fifth column onwars i.e abundance columns
            print("Test Abundances:",testAbundances)
        else:
            testNormSum = 'NA' 
            testAbundances = 'NA'
        
        if CompareDBs  == 'Y':
            ### Query control table
            cur.execute("select * from %s.%s where tag = '%s'" % (db,conTable,tag))####
            temp = cur.fetchall() ### What if tag is not found
            if temp: ## Tag is found 
                conAbundances = ",".join(str(x) for x in temp[0][5:]) ## in comma separated format, from fifth column onwars i.e abundance columns
                conNormSum = temp[0][4]
                fh_out.write('%s,%s,%s,%s\n' % (tag,round((testNormSum/conNormSum),2),testAbundances,conAbundances))
            else: ## Tag is not found in control - cheers
                conLib = len(conTablecol.split(',')) ## The number of control column i.e number of libraries
                emptyList = [0]*int(conLib) ## Create abundance list with zeros
                conAbundances = ",".join(str(x) for x in emptyList)
                conNormSum = temp[0][4]
                fh_out.write('%s,%s,%s,%s\n' % (tag,round((testNormSum-conNormSum),2),testAbundances,conAbundances)) ## NormAbundance is substracted to avoid division to zero error
        
        else:
            if temp:
                NormSumRatio = 0 ## There is no comparision between two tables so no two NormSum to calculate ratios - just fetching abundance for single table - calculate ratio for required libraries manually
                fh_out.write('%s,%s,%s\n' % (tag,NormSumRatio,testAbundances)) ## NormAbundance is substracted to avoid division to zero error
            else:
                NormSumRatio = 0 ## There is no comparision between two tables so no two NormSum to calculate ratios - just fetching abundance for single table - calculate ratio for required libraries manually
                fh_out.write('%s,%s,%s\n' % (tag,NormSumRatio,testAbundances)) ## NormAbundance is substracted to avoid division to zero error
            

    fh_in.close()
    fh_out.close()
    
    return outfile

def tag2Abun2(con,tagfile):
    '''This fetches tag abundance from sRNA DB'''

    cur = con.cursor()
    outfile = tagfile.split('.')[0]+'Abun_runmaster.csv'

    ## Get libs
    cur.execute("select lib_id from %s.library" % (DB))
    temp = cur.fetchall()
    print(temp)

    libs = []
    for x in temp:
        libs.append(x[0])
    print(libs)
    
    fh_out = open(outfile,'w')
    fh_out.write("tag\thits\t%s\t%s\n" % ('\t'.join(str(x)+"_raw" for x in libs),'\t'.join(str(x)+"_norm" for x in libs))) ## For raw and norm abundances

    fh_in = open(tagfile, 'r')

    for i in fh_in:
        rawAbun  = [] ## For every tag
        normAbun = [] ## FOr every tag
        tag = i.strip('\n').replace('-','').translate(str.maketrans("U","T"))## Manipulation if tag is from PARE valid results - reversed and U -> T
        print("'%s' being queried:" % (tag))
        cur.execute("select hits from %s.run_master where tag = '%s'" % (DB,tag))
        temp = cur.fetchall()
        hits = temp[0][0]
        ### Query test table
        for lib in libs:
            cur.execute("select raw_value,norm from %s.run_master where tag = '%s' and lib_id ='%s'" % (DB,tag,lib))####
            temp = cur.fetchall()
            # print ("\nTemp:",temp)
            if temp:
                rawAbun.append(temp[0][0])
                normAbun.append(temp[0][1])
            else:
                rawAbun.append(0)
                normAbun.append(0)

        print("\nTag:%s | Hits:%s" % (tag,hits))
        print("raw abundances",rawAbun)
        print("norm abundances",normAbun)
        # print("%s\t%s\t%s\t%s\n" % (tag,str(hits),'\t'.join(str(x) for x in rawAbun),'\t'.join(str(x) for x in normAbun)))
        fh_out.write("%s\t%s\t%s\t%s\n" % (tag,str(hits),'\t'.join(str(x) for x in rawAbun),'\t'.join(str(x) for x in normAbun)))

    fh_out.close()

    return None
  
def main(db,tagfile,testTable,conTable):
    con = ConnectToDB(server)
    if mode == 0:
        tags2Abun(con,tagfile,testTable,conTable)
    elif mode == 1:
        tag2Abun2(con,tagfile)

if __name__ == '__main__':
    main(db,tagfile,testTable,conTable)
    print('script finished sucessfully')
    sys.exit()

## v01 -> v02
## Added fetching tag abundace from run_master tabel of database - This is imporatant because tag_pos_summary might not have tags for which there are no hits.
## These are kind of important for species with no genome