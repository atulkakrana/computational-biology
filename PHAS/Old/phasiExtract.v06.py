#!/usr/local/bin/python3

###Script to extract phasiRNA from a cluster file for specific library
###Used phased loci name as in PARE validation results or coords in csv format to extract the phasiRNAs
###Usage

import os,sys
import difflib
import time
import mysql.connector as sql
import operator


##########Settings###################

clustfile       = 'All.24Phas.score_p5e-07_sRNA_24_out.cluster'        ## Cluster file from phasing analysis, coul dbe for one lib or concatanated file for all libs
fetchMax        = "Y"                                               ## Fetch the max abundance phasi from tag position summmary, write in separate file
tagPosTable     = "RICE_sbsQIFA_sRNA_TagPosSummNorm"                ## if fetchMax == Y | USe this DB to get tag with max abundance


#### User specifed coords in input file
coordsfile      = 'Table2_21PHAS_AbunAnno.txt'                      ## File with either phased loci ID 'Phas-100_w_7:6353898:6354260' or coords 'w,7,6353898,6354260'
head            = 'Y'
coordsep        ='\t'                                               ## Seprator used in the file
phasedID        = 'N'                                               ## If given than phas ID is used to extract coords, if not than coords specified by user are used
chrcol          = 5                                                 ## As in excel format
startcol        = 6                                                 ## As in excel format
endcol          = 7                                                 ## As in excel format

phasiLenFilter  = 'Y'                                               ## 'Y' then tags of phase length will be written from the cluster | 'N' - All tags will be written
minAbun         = 1                                                 ## Minimum abundance of tag to be written
matchThres      = 0.95                                              ## Ratio of match required by the cluster to phased loci
phase           = 21
startbuff       = 0                                                 ## While extracting sequence through coords2FASTA a buffer is added to start, 
                                                                    ## start position in phased ID has this buffer added, ## so minus this buffer 
                                                                    ## for better matching with real loci
excludeLibs = [2523,2524,2525,2526,2527,2528,4165,4166,4167,4168,4169]   ##  Used if fetchMax == 'Y' | Libs you wish to exclude, your tag position summary table should have libs ids and not lib codes

dataserver  = 'taiji.dbi.udel.edu'
DB          = "kakrana"

############ FUNCTIONS ##############
#####################################

def reader(coordsfile):
    '''
    Reads coordinates file and return a list of PHASE ID and their coords values
    '''

    fh_in = open(coordsfile)
    if head == 'Y':
        fh_in.readline()

    phasCount = 0       ## Count of entries read
    phasList = []       ## List to store results 
    for ent in fh_in:

        # print('\n\nEntry from input file being matched: %s' % (ent))
        coords = ent.split(coordsep)
        #print('coords:',coords)
        if phasedID == 'Y':                                                         ## Extract coords from phased ID
            loci = coords[0].split('_')                                             ## Phas-10_w_3:27574117:27574772
            # phasID = loci[0]
            fusedcoords     = loci[2].split(':')                                    ## 3:27574117:27574772
            get_chr_id      = fusedcoords[0].replace("chr","").replace("Chr","")
            astart          = fusedcoords[1]
            aend            = fusedcoords[2]
            get_start       = int(astart)+startbuff                         ## It should be added to reduce length of loci
            get_end         = int(aend)+1 
            phasID          = '%s_%s_%s' % (get_chr_id,astart,aend)           ## 1 added because when opening a range using start and end, end number is not included in range
            phasCount       += 1
        
        elif phasedID == 'N': ## User specified coords
            # print("File with user specified columns will be used")
            get_chr_id      = coords[chrcol-1]
            astart          = coords[startcol-1]
            aend            = coords[endcol-1] 
            get_start       = int(int(astart)+startbuff)
            get_end         = int(aend)+1                               ## 1 added because when opening a range using start and end, end number is not included in range
            phasID          = '%s_%s_%s' % (get_chr_id,astart,aend)           ## Name to be used in output file
            print("Phased Loci: %s #######################################################" % (phasID))
            phasCount       += 1

        else:
            print("Please input correct value for 'phasedID' or check your file")

        get_value  = (list(range(int(str(get_start)),int(str(get_end)))))
        phasList.append((phasID,get_chr_id,get_start,get_end,get_value))
        
    print("Entries read:%s | Entries cached:%s" % (phasCount,len(phasList)))

    return phasList

def getClust(clustfile,phasList):
    
    fh_in = open(clustfile,'r')
    clusters = fh_in.read().split('>')
    
    resList   = [] ## Store final results as (phas,[(phasiRNA),(PhasiRNA)],[extra info])
    
    phasCount       = 0                                                          ## Total phased loci in file
    uniqMatchCount  = 0                                                          ## Atleast one cluster present for one phased loci
    allMatchCount   = 0                                                          ## Total number of matched cluster
                                                            
    for ent in phasList:                                                         ## Given an entry in coords file
        phasID,get_chr_id,get_start,get_end,get_value = ent
        # print("This is the PhasId: %s | values:%s" % (phasID,get_value))
        print("\n\nPhaseID being queried:%s ##############" % (phasID))
        phasCount +=1 
        print("%s/%s phasID" % (phasCount,len(phasList)))

        ## Find matching cluster
        matchCount = 0                                                          ## Total maching clusters for a phased loci - if same cluster in multiple libraries
        finalMatchList = []                                                     ## Holds best cluster, from multiple libraries
        
        
        for aclust in clusters[1:]:
            tempMatchList   = [] ## To hold results of current matching cluster
            aclust_splt     = aclust.split('\n')
            header          = aclust_splt[0].split()
            clust_id        = header[2]
            chr_id          = header[6].replace("chr","").replace("Chr","")
            start           = header[10]
            end             = int(header[12])+1 ##1 added because when opening a range using start and end, end number is not included in range
            value           = (list(range(int(str(start)),int(str(end)))))
            #print ('Cluster:', (value))
            
            if int(get_chr_id) == int(chr_id):
                sm=difflib.SequenceMatcher(None,get_value,value)
                
                if round(sm.ratio(),2) >= matchThres:
                    ### Matched - phasiRNA from this cluster
                    print ('\nMatching cluster found:%s' % ''.join(header))
                    matchCount +=1
                    
                    phasiCyc = 0 ## Stores phasing cycles
                    phasiSig = 0 ## Stores total abundance of phased sRNAs
                    
                    for i in aclust_splt[1:-1]:## Because header was the first entry of block and not required here, Last entry is always empty
                        # print ("Matched Cluster:\n",i)
                        phasient    = i.split('\t')
                        phasiname   = phasient[4].replace("|","_")
                        phasiseq    = phasient[5]
                        phasilen    = int(phasient[6])
                        phasiabun   = int(phasient[7])
                        phasihits   = int(phasient[10].split("=")[1])
                        phasipos    = int(phasient[3])

                        # print(phasiname,phasiabun,phasiseq,phasilen)
                        tempMatchList.append((phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos))

                        if int(phasilen) == phase:
                            phasiCyc +=1
                            phasiSig += phasiabun
                
                    print("Current Cycles:%s | Current sig. strength:%s" % (phasiCyc,phasiSig))
                    tempMatchList.append((phasiCyc,phasiSig,phasID,clust_id))

                    ## Decide the best and remove other from list ##############################
                    ############################################################################
                    if finalMatchList:
                        ## There exists a previosly matched cluster
                        exist_phasiCyc = finalMatchList[-1][0]
                        exist_phasiSig = finalMatchList[-1][1]
                        print("Existing Cycles:%s | Existing sig. strength:%s" % (exist_phasiCyc,exist_phasiSig))

                        if phasiCyc > exist_phasiCyc: ## New cluster has more cycles
                            del finalMatchList[0:]
                            finalMatchList = list(tempMatchList)
                            print("--- New cluster selected ---")

                        elif phasiCyc == exist_phasiCyc: ## Both have same cycles
                            if phasiSig > exist_phasiSig: ## New one has more total abundance of phased siRNAs
                                del finalMatchList[0:]
                                finalMatchList = list(tempMatchList)
                                print("--- New cluster selected ---")
                        
                        else: ## Existing/old one was long i.e. had more cycles
                            print("Earlier recorded cluster is retained")
                            pass

                    else: ## This is the first cluster
                        finalMatchList = list(tempMatchList)

                    # print("\nFinal Match List:",finalMatchList)
            else:
                # print("No Match with this cluster")
                pass

        resList.append((phasID,finalMatchList)) ## Add best matched entry to the final list, there has to be one best matched entry per PHAS
        
        allMatchCount += matchCount ## Add the matched cluster for each entry
        if matchCount > 0 :
            uniqMatchCount+=1
        
    print("\nTotal phas loci: %s | Matched: %s" % (len(phasList),len(resList)))
    print ("SUMMARY: Phased loci in input file:%s | Loci match threshold: %s |Uniq matched cluster: %s | Total matched clusters found:%s" % (phasCount,matchThres,uniqMatchCount,allMatchCount))
    print("NOTE: If matched clusters more then phased loci that means same cluster was present in different libs\n")
    # print("NOTE: Don't forget to uniq the miRNAs")
    fh_in.close()

    return resList

def writer(resList,con):
    '''
    write the results
    '''

    print("\n\nFunction:Writer\n")

    outfile = clustfile+'thres_%s_phasi.csv' % matchThres                       ## 2437.txt.score_p1e-07_sRNA_21_out.cluster
    fh_out = open(outfile,'w')


    if fetchMax == "Y": ### Fetch max phasi for each loci 
        outfile2 = "fetchMax.txt"
        fh_out2 = open(outfile2,'w')
        fh_out2.write("loci\tmaxtag\tmaxTagAbun\ttotalPhasAbun\n")
        cur= con.cursor()
        queryLibs = prepareQuery(excludeLibs,cur)
        # print(queryLibs) ## Libs whose abindance will be summed to give final abundance of tags
        abunList =   [] ### List to capture tags and abundance for each phased loci

    for ent in resList: ## entry corresponds to one phased loci
        print("\nEntry",ent)
        phasID      = ent[0]
        phasCycles  = ent[1][-1][0]
        phasSig     = ent[1][-1][1]
        clustID     = ent[1][-1][2]
        phasiList   = ent[1][0:-1]
        # print("%s | phasCycles:%s | phasSig:%s" % (phasID,phasCycles,phasSig))

        for i in phasiList:
            print("-Phasi",i)  ## Phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos
            tag = i[2]
            
            ## Write phasiRNAs 
            if phasiLenFilter == 'Y': ## If tags filter is ON
                if int(i[3]) == int(phase) and int(i[1]) >= minAbun: ### Size specified in settings
                    print(phasID,clustID,i[0],i[1],i[2])
                    fh_out.write('>%s_Clust%s_%s,%s,%s\n%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],tag,i[1],i[4]) )  ## phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun

            else:
                if int(i[1]) > minAbun:
                    print(phasID,clustID,i[0],i[1],i[2])
                    fh_out.write('>%s_Clust%s_%s,%s,%sn%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],tag,i[1],i[4]) )  ## phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun
                    pass

            ## Get max abundance phasiRNA for specified phase
            if fetchMax == "Y":
                if len(tag) == int(phase): ### Size specified in settings
                    atag,abun_sum = getAbundance(cur,tag,phase,queryLibs)
                    abunList.append((atag,abun_sum))

        ## fetch max results fro this phased loci
        if fetchMax == "Y":
            abunList_sort = sorted(abunList, key=operator.itemgetter(1),reverse=True) ## Sort list on abundances to get max abundant phasiRNA
            print("\nExample sorted values:%s" % (abunList_sort[0:10]))
            maxTag      = abunList_sort[0][0]   ## max abundant tag
            maxAbun     = abunList_sort[0][1]   ## Abundance of max abundant tag
            totalAbun   = sum(int(i[1]) for i in abunList_sort) ## Phased abundance

            fh_out2.write("%s\t%s\t%s\t%s\n" % (phasID,maxTag,maxAbun,totalAbun))
            abunList =   [] ## Empty before next phased loci


    fh_out.close()
    fh_out2.close()

    return outfile,outfile2

def getAbundance(cur,tag,phase,queryLibs):
    '''Input is tag for each loci and out put is tag with maximmum abumdance and sum of phasiRNAs'''
    
    cur.execute("SELECT tag,%s FROM %s.%s where tag = '%s'" % (queryLibs,DB,tagPosTable,tag))### Convert intergenic to gene name so as to get strand
    info = cur.fetchall()
    print("From tagpos table",info)
    atag,abun_sum = info[0]

    return tag,abun_sum
            
def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### get column names
    columns = [] ## empty list
    cur.execute("describe %s.%s" % (DB,tagPosTable))
    tablefield = cur.fetchall()
    # print (tablefield)
    for i in tablefield:
        col = i[0]
        columns.append(col)

    libs = columns[8:-2] ## norm_sum and max_norm are last two columns
    lastTwo = columns[-2:] ## norm_sum and max_norm

    if excludeLibs:
        print("\n\nLibs specifed in excludeLibs %s will be skipped\n\n" % (excludeLibs))
        selectLibs = [] ## Columns excluding the unwanted libraries
        excludeLibs_s = [str(i) for i in excludeLibs] ## Converting all entries in exclude list to string fro matching below
        
        for i in libs:
            lib = i.split('_')[1] ## Columns in tag position summary - 'lib_2499' or 'max_norm'
            print(lib)
            if str(lib) not in excludeLibs_s: ## Tested OK
                selectLibs.append(i)
            else:
                print("excluded:",i)
                # sys.exit()
                pass
        
        finalLibs = selectLibs ## Norm_Sum and max_norm are not included

    else:
        finalLibs = libs ## ## Norm_Sum and max_norm are not included

    # print("finalLibs:%s" % (finalLibs))
    
    lib_col =",".join(str(x) for x in finalLibs)### Manually mentioned strating lib column - should work on all tag position summary tables
    
    print("Library Columns:",lib_col)
    # queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    queryLibs = "%s" % lib_col.replace(",","+")

    return queryLibs

def ConnectToDB(server, infile):
    
    ##infile values are '0' when you dont want to pulaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE bioinfo_data.mir_page_results FIELDS TERMINATED BY ','")
    
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

######## MAIN #######################
#####################################
def main():
    phasList = reader(coordsfile)
    time.sleep(1)

    resList = getClust(clustfile,phasList)

    con = ConnectToDB(dataserver,0)
    results,results2 = writer(resList,con)


if __name__ == '__main__':
    main()
    print('\n**script finished sucessfully**\n')
    sys.exit()

#### v0.1 -> v0.3
## updated phas matching with chromosome not part of value instead matched before any ratio computation
## If clusters for all libs are concatanated and same phased locus matches clusters from different libs then the longest phase cycle is slected and is phase cycle are same most abundant is selcted

### v0.3 -> 0.4 [major][stable]
### Organized script into functions
### Added functionality to fetch max abundance phasiRNA for each locus

## v0.4 -> v05 
### Fixed bug fetching abundance of most abundant tag, norm sum was being fetched unstead of max norm for tag
## Added fetching total abundance of phased tags - Useful for plotting

## v05 -> v06
## Added functionality to exclude libraries for which abundnace is not required
## Made design more modular for functionality