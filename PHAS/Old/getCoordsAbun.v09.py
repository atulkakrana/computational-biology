#!/usr/local/bin/python3
## Script by : Atul - kakrana@udel.edu
###Coordinates used to fetch abundance of sRNAs from different library of tag summary table

import sys
import mysql.connector as sql

####### Settings #######################################
########################################################
mode = 0                                                        ## 0: Get summed sRNA abundance for specific size and noise 
                                                                ## 1: Get sRNA abundance for comparision among sizes, sizes from both mainLen and noise Len used
                                                                ## 2: Get sRNA tag counts for comaprision among sizes, size from both mainLen and noise Len used
## DB and Table
genomeDB    = 'RICE_MSU7'                                       ## Just to use for cluster name i.e that can be used directly on our website
tagPosTable = 'RICE_sbsQIFA2_sRNATagPosSummNorm'                ## Tag position summary table required
libType     = 1                                                 ## 0: Lib_ids (4518) | 1: lib_code ('leaf_1')
excludeLibs = ['73','365']                                      ## If you wish to exclude a few libs then enter lib_id here, your tag position summary table shoul dhave libs ids and not lib codes

## Coord file
coords      = 'Table1_24PHAS_AbunAnno.phasi.v4.txt'
sep         = '\t'                                              ## Comma: ',' Tab: '\t'
header      = 'Y'

## Input columns
useStrand   = 'N'                                               ## Not used for phased loci as sRNA could be from either strand
chrpos      = '6'                                                    ## Excel format
strandpos   = '6'                                               ## Excel format
startpos    = '7'                                               ## End coordinate- excel format
endpos      = '8'                                               ## Start coordinate - excel format

## Small RNA lengths 
mainLen     = [24,24]                                           ## Lengths of sRNAs to calculate abundance and used to compute ratio in mode 0
noiseLen    = [20,22,23,21]                                     ## List of different sizes, used as noise to compute main to noise ratio in mode 0. In mode 1 both main and nois elist are combined and abundances reported

## 5' and 3' Buffer in nucleotide
startbuff   = 0                                                 ## This will be used to reduce are back to actual Loci - 0 for PARE
endbuff     = 0

## Dev options #########################################
########################################################
server      = 'tarkan.ddpsc.org'
db          = 'kakrana'

########################################################

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

def Coords2Abun(con,coords,db,excludeLibs):
    cur = con.cursor()
    lib_col,queryLibs = prepareQuery(excludeLibs,cur)
    
    ## Outfile
    file_out    = coords.split('.')[0]+'Abun.csv'
    fh_out      = open(file_out,'w')
    fh_out.write('Loci,Cluster,Ratio,%s,%s\n' % (lib_col,lib_col)) ## Header 
    
    ## Coords file
    fh_in = open(coords,'r')
    if header == 'Y':
        fh_in.readline()
    
    for i in fh_in:
        # cur = con.cursor()
        coord   = i.strip('\n').split(sep)
        print ('\nRow',coord)
        start   = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
        start2  = int(start)+startbuff ### Reverse the effect on co-ordinate made during fetching fasta
        end     = str(coord[int(endpos)-1]).strip()
        end2    = int(end)-endbuff
        length  = int(end2)-int(start2) ###Extra 1 to include last position, buffer of 100 bp from start and end
        chr_id  = str(coord[int(chrpos)-1]).strip().replace('chr','').replace("Chr",'')
        strand  = str(coord[int(strandpos)-1]).strip()
        print ('This is chr_id:%s and Start:%s and End:%s | Length of loci: %s bp ' % (chr_id,str(start2),int(end2),length))

        ## Prepare queries for input size for both main and noise
        mainSizeQuery   = ''
        noiseSizeQuery  = ''
        for i in mainLen:
            if mainSizeQuery == '':
                mainSizeQuery = 'len = %s' % (i)
            else:
                mainSizeQuery = mainSizeQuery + ' or len = %s' % (i)

        for i in noiseLen:
            if noiseSizeQuery == '':
                noiseSizeQuery = 'len = %s' % (i)
            else:
                noiseSizeQuery = noiseSizeQuery + ' or len = %s' % (i)
        print ('%s | %s' % (mainSizeQuery,noiseSizeQuery))
        # sys.exit()

        if useStrand == 'Y':
            ## Get main size
            cur.execute("select %s from %s.%s where chr_id = %s and strand = '%s' and (%s) and (position between %s and %s)" % (queryLibs,db,tagPosTable,chr_id,strand,mainSizeQuery,str(start2),str(end2)))####
            temp = cur.fetchall()
            main_abun = list(map(int, temp[0]))
            # print(main_abun)
            
            ## Get noise size
            cur.execute("select %s from %s.%s where chr_id = %s and strand = '%s' and (%s) and (position between %s and %s)" % (queryLibs,db,tagPosTable,chr_id,strand,noiseSizeQuery,str(start2),str(end2)))####
            temp2 = cur.fetchall()
            noise_abun = list(map(int, temp2[0]))
            # print(noise_abun)
            
            lociname    = '%s_%s_%s_%s' % (chr_id,strand,start,end)
            clusterName = '%s.%s.%s.%s' % (genomeDB,chr_id,start,end)
        
        else:
            ## Get main abun
            cur.execute("select %s from %s.%s where chr_id = %s and (%s) and (position between %s and %s)" % (queryLibs,db,tagPosTable,chr_id,mainSizeQuery,str(start2),str(end2)))####
            temp = cur.fetchall()
            main_abun = list(map(int, temp[0]))
            print(main_abun)
        
            ## Get noise abun
            cur.execute("select %s from %s.%s where chr_id = %s and (%s) and (position between %s and %s)" % (queryLibs,db,tagPosTable,chr_id,noiseSizeQuery,str(start2),str(end2)))####
            temp2 = cur.fetchall()
            noise_abun = list(map(int, temp2[0]))
            print(noise_abun)
            
            lociname    = '%s_%s_%s' % (chr_id,start,end)
            clusterName = '%s.%s.%s.%s' % (genomeDB,chr_id,start,end)
        

        mainNormSum     = main_abun[-2]
        abundances      = ",".join(str(x) for x in main_abun) ## in comma separated format
        mainNormSum2    = noise_abun[-2]
        abundances2     = ",".join(str(x) for x in noise_abun) ## for different size  i.e 22,24 or what ever defined in settings
        print('Norm sum for main:%s | Norm sum for noise:%s' % (mainNormSum,mainNormSum2))
        
        if mainNormSum != None : ## There are sRNAs of desired length in this cluster
            if mainNormSum2 != None:
                ratioNormSum = int(mainNormSum)/(int(mainNormSum+int(mainNormSum2)))
            else:
                ratioNormSum=1
        else:
            ratioNormSum = 0
        fh_out.write('%s,%s,%s,%s,%s\n' % (lociname,clusterName,str(ratioNormSum),abundances,abundances2))
        
        #print (temp)
    fh_out.close()
    fh_in.close()
    
    return file_out

def compCoordAbun(con,coords,db,excludeLibs):

    cur = con.cursor()
    lib_col,queryLibs = prepareQuery(excludeLibs,cur)

    combSize    = mainLen+noiseLen ## All different lengths for which abundance is required
    allLen      = set(sorted(combSize)) ## Sort according to size from low to high, set to remove redundant sizes

    file_out    = ('%s_%s_CompAbun.csv' % (coords.split('.')[0],'_'.join(str(x) for x in allLen))) ## Complete file with lib wise abundances and norm sum of all sizes
    file_out2   = ('%s_%s_CompAbun2.csv'% (coords.split('.')[0],'_'.join(str(x) for x in allLen))) ## File with just the norm sum of varoius sizes
    fh_out      = open(file_out,'w')
    fh_out2     = open(file_out2,'w')

    
    ##Write file header
    header1 = '' ## Blank Header
    header2 = '' ## Blank header for outfile2
    for i in allLen:
        if header1 == '':
            header1 = '%s_%s' % (i,lib_col.replace(",",",%s_" % (i)))
            # print("\n",header1)
        else:
            header1 = header1+','+'%s_%s' % (i,lib_col.replace(",",",%s_" % (i)))
            # print("\n",header1)

        if header2 == '':
            header2 = str(i)+'_norm_sum'
        else:
            header2 = header2+','+str(i)+'_norm_sum'

    fh_out.write('Loci,Cluster,%s\n' % (header1))
    fh_out2.write('Loci,Cluster,%s\n' % (header2))
    
    fh_in = open(coords,'r')
    if header == 'Y':
        fh_in.readline()

    ## For each coordinates
    for i in fh_in:
        # cur = con.cursor()
        coord   = i.strip('\n').split(sep)
        print ('\nRow',coord)
        start   = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
        start2  = int(start)+startbuff ### Reverse the effect on co-ordinate made during fetching fasta
        end     = str(coord[int(endpos)-1]).strip()
        end2    = int(end)-endbuff
        length  = int(end2)-int(start2) ###Extra 1 to include last position, buffer of 100 bp from start and end
        chr_id  = str(coord[int(chrpos)-1]).strip().replace('chr','').replace("Chr",'')
        strand  = str(coord[int(strandpos)-1]).strip()
        print ('This is chr_id:%s and Start:%s and End:%s | Length of loci: %s bp ' % (chr_id,str(start2),int(end2),length))
        queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(") ## String manipulation to automativcally convert library names into query format below - used in all queries

        resList = [] ## To store all results, required for result format 1
        resList2 = [] ## To store just Norm sum for all sizes required for result format 2
        
        if useStrand == 'Y':
            for size in allLen:
                cur.execute("select %s from %s.%s where chr_id = %s and strand = '%s' and (len = %s) and (position between %s and %s)" % (queryLibs,db,tagPosTable,chr_id,strand,size,str(start2),str(end2)))####
                temp        = cur.fetchall()
                aabun = list(map(int, temp2[0]))
                # print(aabun)
                
                abundances  = ",".join(str(x) for x in aabun) ## in comma separated format
                mainNormSum = aabun[-2]

                resList.append(abundances)
                resList2.append(mainNormSum)

                ###Get alen[1] i.e 24,22 or whatever mentioned in settings
                #cur.execute("select sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589),norm_sum,max_norm from %s.%s where chr_id = %s and strand = '%s' and len = %s and (position between %s and %s)" % (db,tagPosTable,chr_id, strand, alen[1],str(start2),str(end2)))####

            lociname        = '%s_%s_%s_%s' % (chr_id,strand,start,end)
            clusterName     = '%s.%s.%s.%s' % (genomeDB,chr_id,start,end)
        
        else:
            for size in allLen:
                # print("Size of sRNAs being collected:%s" % (str(size)))
                #cur.execute("select sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589),norm_sum,max_norm from %s.%s where chr_id = %s and len = %s and (position between %s and %s)" % (db,tagPosTable,chr_id, alen[0],str(start2),str(end2)))####
                cur.execute("select %s from %s.%s where chr_id = %s and (len = %s) and (position between %s and %s)" % (queryLibs,db,tagPosTable,chr_id,size,str(start2),str(end2)))####
                temp = cur.fetchall()
                aabun = list(map(int, temp[0]))
                # print(aabun)

                abundances = ",".join(str(x) for x in aabun) ## in comma separated format
                mainNormSum = aabun[-2]

                resList.append(abundances)
                resList2.append(mainNormSum)
        
            lociname = '%s_%s_%s' % (chr_id,start,end)
            clusterName = '%s.%s.%s.%s' % (genomeDB,chr_id,start,end)

        allabundances   = ",".join(str(x) for x in resList)
        allNormSum      = ",".join(str(x) for x in resList2)
        
        fh_out.write('%s,%s,%s\n' % (lociname,clusterName,allabundances))
        fh_out2.write('%s,%s,%s\n' % (lociname,clusterName,allNormSum))

    fh_out.close()
    fh_out2.close()

    return file_out,file_out2

def coordSizeCounts(con,coords,db,excludeLibs):
    ''' Module reports number of libwise number of tags and sum of tags of same size'''

    cur = con.cursor()
    lib_col,trash = prepareQuery(excludeLibs,cur) ## queryLibs will be prepepared below

    combSize    = mainLen+noiseLen ## All different lengths for which abundance is required
    allLen      = set(sorted(combSize)) ## Sort according to soze from low to high, set to remove redundant sizes

    file_out    = ('%s_%s_CompCounts.csv' % (coords.split('.')[0],'_'.join(str(x) for x in allLen))) ## Complete file with lib wise abundances and norm sum of all sizes
    file_out2   = ('%s_%s_CompCounts2.csv'% (coords.split('.')[0],'_'.join(str(x) for x in allLen))) ## File with just the norm sum of varoius sizes
    fh_out      = open(file_out,'w')
    fh_out2     = open(file_out2,'w')

    ## Write file header #######################################################
    header1 = '' ## Blank Header
    header2 = '' ## Blank header for outfile2
    for i in allLen:
        if header1 == '':
            header1 = '%s_%s' % (i,lib_col.replace(",",",%s_" % (i)))
            # print("\n",header1)
        else:
            header1 = header1+','+'%s_%s' % (i,lib_col.replace(",",",%s_" % (i)))
            # print("\n",header1)

        if header2 == '':
            header2 = str(i)+'_total_counts'
        else:
            header2 = header2+','+str(i)+'_total_counts'

    fh_out.write('Loci,Cluster,%s\n' % (header1))
    fh_out2.write('Loci,Cluster,%s\n' % (header2))

    ############################################################################

    fh_in = open(coords,'r')
    if header == 'Y':
        fh_in.readline()

    ## For each coordinates
    for i in fh_in:
        # cur = con.cursor()
        coord   = i.strip('\n').split(sep)
        print ('\nRow',coord)
        start   = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
        start2  = int(start)+startbuff ### Reverse the effect on co-ordinate made during fetching fasta
        end     = str(coord[int(endpos)-1]).strip()
        end2    = int(end)-endbuff
        length  = int(end2)-int(start2) ###Extra 1 to include last position, buffer of 100 bp from start and end
        chr_id  = str(coord[int(chrpos)-1]).strip().replace('chr','').replace("Chr",'')
        strand  = str(coord[int(strandpos)-1]).strip() ## Used only if "useStrand = Y"
        print ('This is chr_id:%s and Start:%s and End:%s | Length of loci: %s bp ' % (chr_id,str(start2),int(end2),length))
    
        queryLibs = 'COUNT(%s)' % lib_col.replace(",","),COUNTS(") ## String manipulation to automatically convert library names into query format below - used in all queries
        # print ('These are libraries to be queried',queryLibs)

        ### lib_76 >0 or lib_78 > 0
        querycond1 = '%s > 0' % lib_col.replace(","," > 0 OR ")
        # print ("This is querycondition1",querycond1)

        resList = [] ## To store all results, required for result format 1
        resList2 = [] ## To store just Norm sum for all sizes required for result format 2

        if useStrand == 'Y':
            for size in allLen:
                # print("Size of sRNAs being collected:%s" % (str(size)))

                cur.execute("select count(*) from %s.%s where chr_id = %s and strand = '%s' and (len = %s) and (position between %s and %s) and (%s)" % (db,tagPosTable,chr_id,strand,size,str(start2),str(end2),querycon1))
                temp = cur.fetchall()
                sizeCount =  temp[0][0]
                resList.append(sizeCount)

            lociname = '%s_%s_%s_%s' % (chr_id,strand,start,end)
            clusterName = '%s.%s.%s.%s' % (genomeDB,chr_id,start,end)

        else:
            for size in allLen:
                # print("Size of sRNAs being collected:%s" % (str(size)))
                cur.execute("select count(*) from %s.%s where chr_id = %s and (len = %s) and (position between %s and %s)" % (db,tagPosTable,chr_id,size,str(start2),str(end2)))
                temp = cur.fetchall()
                sizeCount =  temp[0][0]
                resList.append(sizeCount)

            lociname = '%s_%s_%s_%s' % (chr_id,strand,start,end)
            clusterName = '%s.%s.%s.%s' % (genomeDB,chr_id,start,end)


        allSizeCounts = ",".join(str(x) for x in resList)
        
        fh_out2.write('%s,%s,%s\n' % (lociname,clusterName,allSizeCounts))

    fh_out.close()
    fh_out2.close()


    return file_out,file_out2

def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### get column names
    columns = [] ## empty list
    cur.execute("describe %s.%s" % (db,tagPosTable))
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
    queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    # sumLibs = "%s" % lib_col.replace(",","+") ## From phasiMax
    # queryLibs = "%s" % lib_col.replace(",",",") ## From phasiMax
    # print("\nThese are sumLibs:",sumLibs)
    print("\nThis is query Libs:",queryLibs)
    # sys.exit()

    return lib_col,queryLibs

def main(db,coords):
    con = ConnectToDB(server)

    if mode == 0:
        resFile1 = Coords2Abun(con,coords,db,excludeLibs)
    elif mode == 1:
        resFile1,resFile2 = compCoordAbun(con,coords,db,excludeLibs)
    elif mode == 2:
        resFile1,resFile2 = coordSizeCounts(con,coords,db,excludeLibs)
    else:
        print ("Please enter correct mode: 0 or 1\n")
        print ("If you wish to get abundances for a specific size and compare with others then use mode = 0")
        print ("Else if you wish to compare abundances of different size tags in a coordinate then use mode = 1")
        sys.exit()

if __name__ == '__main__':
    main(db,coords)
    print('\n**Target sRNA length:%s nt | Size range for noise ratio calculation: %s **' % (str(mainLen),str(noiseLen)))
    print('**Make sure you input the correct length for target sRNA size and range of other sRNA size to calculate Ratio**\n')
    print('**Script Finished Succesfully')
    sys.exit()

######### Developemnt Logs #################

## v01-> v02
## use string manipulation to covert column names used in header to convert to mysql query for ex.. column header = lib_2435,lib_2436,lib_2437 convert to sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438)

##v02 -> v03
## Fixed a bug where DB and table name was hard coded
## Added functionality to provide a range of sRNA sizes that is to be used to calculate noise
## Removed 'norm_sum,max_norm' from select of mysql statement as we are quering tags within a pair of coords so these needs to be summend for all which is already included in querySum
## Added cluster name

## v03 -> v04
## Added functionality to have abundance for tags of multiple lengths - This functionality will allow usage of this script to fetch PARE abundance for cleavge site from mutiple libraries
## Added IF loop in case main abundance or noise abundance is = 0
## Outfile name now contains abundance

## v04 -> v05
## Fix critical bug that was not using all lengths in noiseLen list to compute the abundances and hence gave wrong main to noise ratio in earlier versions
## Added funtionality to exclude libraries
## Added new mode [1], if the abunances of all different sizes are required for loci, uses length is both mainLen and noiseLen list

## v05 -> v065
## modified the excludeLibs part to work on string names
## Fixed the excludeLibs not working bug

## v065 -> v07
## Added functionality work on files woth no Chromosome column i.e. if column has 'chr3/Chr3' like inverted repeat files then it will work

## v7 -> v08
## Added functionality  (mode-2) to compute tags in coordinates

## v08 -> v09 [stable]
## Added prepare query module from fetchmax scripts
## updated all modules to recieve lib_cols and queryLibs from this new module
## Cleaner and slimmer functions

### Extra
####sum(lib_2435),sum(lib_2436),sum(lib_2437),sum(lib_2438),sum(lib_2439),sum(lib_2440),sum(lib_2441),sum(lib_2495),sum(lib_2496),sum(lib_2497),sum(lib_2498),sum(lib_2499),sum(lib_2500),sum(lib_2501),sum(lib_2502),sum(lib_2503),sum(lib_2504),sum(lib_2505),sum(lib_2506),sum(lib_2507),sum(lib_2508),sum(lib_2509),sum(lib_2510),sum(lib_2511),sum(lib_2512),sum(lib_2513),sum(lib_2514),sum(lib_2515),sum(lib_2516),sum(lib_2517),sum(lib_2518),sum(lib_2519),sum(lib_2520),sum(lib_2521),sum(lib_2522),sum(lib_2596),sum(lib_2597),sum(lib_2598),sum(lib_2599),sum(lib_2600),sum(lib_2601),sum(lib_2602),sum(lib_2523),sum(lib_2524),sum(lib_2525),sum(lib_2526),sum(lib_2527),sum(lib_2528),sum(lib_1588),sum(lib_1589)
