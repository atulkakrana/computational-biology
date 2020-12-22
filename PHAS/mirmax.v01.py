#!/usr/local/bin/python3

import os,sys,time,sqlite3,operator
import mysql.connector as sql
from collections import Counter

## This script is written for to fetch two max tag from miRNAs loci

#### USER SETTINGS ######

coordsF = "zm.mod.uniq.txt"
mname   = 4         ## Excel format
mchr    = 1        ## Excel format
mstart  = 2        ## Excel format
mend    = 3        ## Excel format
mstrand = 6        ## Excel format
msize   = 7        ## Excel format
delim   = "\t"      ## Delimiter
msizes  = [21,22,20]   ## miRNAs for these length to be used for analysis
head    = "Y"
excludeLibs = []                                      ## If you wish to exclude a few libs then enter lib_id here, your tag position summary table shoul dhave libs ids and not lib codes

ntags = 2       ## 1: report ratios based on most abundant tag 2: on basis of two topmost abundant tags

### DEVELOPER SETTINGS ####
dataserver      = "tarkan.ddpsc.org"
msDB            = "kakrana"
tagPosTable     = "MAIZE_pub2_sRNATagPosSummNorm"
tagSizeBuffer   = 1 ## Tags of size +/-buffer to miRNA will be used in abundance calculation, because miRNA-3p are usally one nt short 


#### Functions ####

def miRmax(coords):
    '''This function computes ratio of top one and top two sRNA tags
    from miRNA precursor, against all tags in that precursor'''

    con = ConnectToDB(dataserver)
    cur = con.cursor()
    lib_col   = prepareQuery(excludeLibs,cur)

    ## OutFile
    outFile = "mirFetchMax.txt"
    fh_out  = open(outFile,'w')
    fh_out.write("miRNA\tmaxTag\tmaxtagAbun\tmaxtag2Abun\ttagsSum\tmaxratio\tmaxratio2\n")

    ### Compute ratios of max tags in miRNA precursors
    for i in coords:
        aname,achr,astart,aend,astrand,alen = i
        print("\nEnt:",aname,achr,astart,aend,astrand,alen)

        minlen = alen-1 ## If miRNA is 5p then this length corresponds to 3p
        maxlen = alen+1 ## If miRNA is 3p then this length corresponds to 5p
        cur.execute("select tag,%s from %s.%s where chr_id = %s and strand = '%s' and (len between %s and %s) and (position between %s and %s)" % (lib_col,msDB,tagPosTable,str(achr),str(astrand),str(minlen),str(maxlen),str(astart),str(aend)))####
        temp        = cur.fetchall()
        print("These are the tags:",temp)

        ## Check - If miRNA was expressed in provided linraries
        if not temp:
            continue
        else:
            pass

        ## Sum Abundances of a tag from all libs
        precursorTags   = [] ## A temp list that stores merged abundances from all libs
        tagSet          = set() ## To maintain a uniq set, redundant tags found, problem with how tag position summary is computed
        allsum          = 0
        for ent in temp:
            atag = ent[0]
            if atag not in tagSet:
                asum = sum(ent[1:-1])
                print(atag,asum)
                precursorTags.append((atag,int(asum)))
                tagSet.add(atag)
            else:
                print("Duplicate tag from same locus exists")
                # sys.exit()
                pass

        ## Get abundance of top two tags and toal abundance of all tags in regions
        precursorTags_s = sorted(precursorTags,key=operator.itemgetter(1),reverse=True) ## CHR_ID, strand and start 
        max1 = precursorTags_s[0]

        if len(precursorTags_s) > 1: ## In few cases there is no other tag of miRNA size, so max2 set to 0
            max2 = precursorTags_s[1]
        else:
            max2 = ["None",0]
        tagsSum = sum(i[1] for i in precursorTags_s)
        print("Max tag:",max1)
        print("Second max:",max2)
        print("Sum of all tags:",tagsSum)

        ## Compute rations and write to file
        if ntags == 1:
            aratio1 = round(max1[1]/tagsSum,2)
            print(aname,max1[0],max1[1],max2[1],tagsSum,aratio1)
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (aname,max1[0],str(max1[1]),str(max2[1]),str(tagsSum),str(aratio1))) ### Name of miRNA entry, max tag, maxtag abun, maxtag2 abun, sum of all tags and aratio
        elif ntags == 2:
            aratio1 = round(max1[1]/tagsSum,2)
            aratio2 = round((max2[1]+max1[1])/tagsSum,2)
            print(aname,max1[0],max1[1],max2[1],tagsSum,aratio1,aratio2)
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (aname,max1[0],str(max1[1]),str(max2[1]),str(tagsSum),str(aratio1),str(aratio2))) ### Name of miRNA entry, max tag, maxtag abun, maxtag2 abun, sum of all tags and aratio
        else:
            print("Please input correct value for ntags")
            sys.exit()

        ## Write results

    fh_out.close()

    return outFile

def coordsParser(coordsF):
    '''parses input PHAS file'''

    print("\nFunction: coordsParser")

    fh_in   = open(coordsF,'r')
    if head == 'Y':
        fh_in.readline()
    coordsRead = fh_in.readlines()

    coordsList = []
    for i in coordsRead:
        print("Ent:",i)
        ent     = i.split(delim)
        aname   = ent[mname-1]
        achr    = int(ent[mchr-1].replace("chr",""))
        astart  = int(ent[mstart-1])
        aend    = int(ent[mend-1])
        asize   = int(ent[msize-1])
        astrand = ent[mstrand-1].translate(str.maketrans("+-","wc"))
        print("Phas Name %s | chr:%s | start:%s | end:%s | astrand:%s | Len:%s" % (aname,achr,astart,aend,astrand,asize))
        if asize in msizes:
            coordsList.append(((aname,achr,astart,aend,astrand,asize)))

    print("A list of query features prepared with entries:%s" % len(coordsList))
    print("Exiting function - coordsParser\n")
    
    return coordsList

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

def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### get column names
    columns = [] ## empty list
    cur.execute("describe %s.%s" % (msDB,tagPosTable))
    tablefield = cur.fetchall()
    # print("\nTable fields:",tablefield)
    for i in tablefield:
        col = i[0]
        columns.append(col)

    libs    = columns[8:-2] ## norm_sum and max_norm are last two columns
    lastTwo = columns[-2:] ## norm_sum and max_norm
    print("\nLibs:",libs)

    if excludeLibs:
        print("\n\nLibs specifed in excludeLibs %s will be skipped\n\n" % (excludeLibs))
        selectLibs = [] ## Columns excluding the unwanted libraries
        excludeLibs_s = [str(i) for i in excludeLibs] ## Converting all entries in exclude list to string fro matching below
        
        for i in libs:

            ## Check if user made mistake in libType - as that will give results for all entries
            if type(i) is int and libType == 1:
                print("You seem to have input lib_id and chosen wrong libType")
                print("Check libType and excludeLibs match - Script will exit now")
                sys.exit()
            elif type(i) is str and libType == 0:
                print("You seem to have input lib_id and chosen wrong libType")
                print("Check libType and excludeLibs match - Script will exit now")
                sys.exit()
            else:
                print("All seems well")
                pass

            ## Get lib_if if libType == 0
            if libType == 0:
                lib = i.split('_')[1] ## Closed to match lib_codes too - Columns in tag position summary - 'lib_2499' or 'max_norm'
                print(lib)
            else:
                lib = i
                pass

            ### Filter libraries
            if str(lib) not in excludeLibs_s: ## Tested OK
                selectLibs.append(i)
            else:
                print("excluded:",i)
                # sys.exit()
                pass
        
        finalLibs = selectLibs+lastTwo ## Norm_Sum and max_norm are included

    else:
        finalLibs = libs+lastTwo ## Norm_Sum and max_norm are included

    # print("finalLibs:%s" % (finalLibs))
    
    lib_col =",".join(str(x) for x in finalLibs)### Manually mentioned strating lib column - should work on all tag position summary tables
    
    print("Library Columns:",lib_col)
    # queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    # sumLibs = "%s" % lib_col.replace(",","+") ## From phasiMax
    # queryLibs = "%s" % lib_col.replace(",",",") ## From phasiMax
    # print("\nThese are sumLibs:",sumLibs)
    # print("\nThis is query Libs:",queryLibs)
    # sys.exit()

    return lib_col

def main():

    coordsL = coordsParser(coordsF)
    resList = miRmax(coordsL)

    return None

if __name__ == '__main__':
    main()
    print("Script finished sucessfully")
    sys.exit()


### v01


###
### 3p are usually int shorter than miRNAs, so query to get abundance need to be modified so that +/-1 nt of miRNAs length
### should be used in query