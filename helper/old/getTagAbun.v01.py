#!/usr/local/bin/python3

### Script to compare sRNA abundances in two databases using summary two table

###### IMPORTS ######################
import sys
import mysql.connector as sql
#####################################


##### INPUT #########################
testTable       = 'RICE_sbsQIFA_sRNA_summary2'
CompareDBs      = 'N'                               ## Fetch abundance from two tables each for one DB for comparision
conTable        = 'RICE_sbsQIFA_sRNA_summary2'
db              = 'kakrana'
server          = 'taiji.dbi.udel.edu'
tagfile         = 'sRNA.txt'                        ## Text file with tags

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

    ###get column anmes
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
        tag = i.strip('\n').replace('-','')[::-1].translate(str.maketrans("U","T"))## Manipulation if tag is from PARE valid results - reversed and U -> T
        print("'%s' being queried:" % (tag))
        ### Query test table
        cur.execute("select * from %s.%s where tag = '%s'" % (db,testTable,tag))####
        temp = cur.fetchall()
        print (temp)
        testNormSum = temp[0][4]
        testAbundances = ",".join(str(x) for x in temp[0][5:]) ## in comma separated format, from fifth column onwars i.e abundance columns
        
        if CompareDBs  == 'Y':
            ### Query control table
            cur.execute("select * from %s.%s where tag = '%s'" % (db,conTable,tag))####
            temp = cur.fetchall() ### What if tag is not found
            if temp:## Tag is found 
                conAbundances = ",".join(str(x) for x in temp[0][5:]) ## in comma separated format, from fifth column onwars i.e abundance columns
                conNormSum = temp[0][4]
                fh_out.write('%s,%s,%s,%s\n' % (tag,round((testNormSum/conNormSum),2),testAbundances,conAbundances))
            else: ## Tag is not found in control - cheers
                conLib = len(conTablecol.split(',')) ## The number of control column i.e number of libraries
                emptyList = [0]*int(conLib) ## Create abundance list with zeros
                conAbundances = ",".join(str(x) for x in emptyList)
                conNormSum = temp[0][4]
                fh_out.write('%s,%s,%s,%s\n' % (tag,round((testNormSum-conNormSum),2),testAbundances,conAbundances)) ## NormAbundance is substracted to avoid division to zero error
        
        else: ##
            NormSumRatio = 0 ## There is no comparision between two tables so no two NormSum to calculate ratios - just fetching abundance for single table - calculate ratio for required libraries manually
            fh_out.write('%s,%s,%s\n' % (tag,NormSumRatio,testAbundances)) ## NormAbundance is substracted to avoid division to zero error
            

    fh_in.close()
    fh_out.close()
    
    return outfile



def main(db,tagfile,testTable,conTable):
    con = ConnectToDB(server)
    tags2Abun(con,tagfile,testTable,conTable)
    
    
if __name__ == '__main__':
    main(db,tagfile,testTable,conTable)
    print('script finished sucessfully')
    sys.exit()
