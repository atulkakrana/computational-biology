#!/usr/local/bin/python3

###Script to extract phasiRNA from a cluster file for specific library
###Used phased loci name as in PARE validation results or coords in csv format to extract the phasiRNAs
###Usage

import os,sys
import difflib
import time
import mysql.connector as sql
import operator
from operator import add


##########Settings###################

clustfile       = 'All.24Phas.score_p5e-07_sRNA_24_out.cluster'        ## Cluster file from phasing analysis, coul dbe for one lib or concatanated file for all libs
fetchMax        = 1                                                ## Fetch the max abundance phasi from tag position summmary, write in separate file
fetchLibAbun    = 1                                                ## Fetch lib-wise abudnaces using all the tags
tagPosTable     = "DAYLILY_priv_sRNA_TagPosSummNorm"                ## if fetchMax == Y | Use this DB to get tag with max abundance


#### User specifed coords in input file
coordsfile      = "Final_PHAS_Loci_5e-07_ALL.csv"                   ## File with either phased loci ID 'Phas-100_w_7:6353898:6354260' or coords 'w,7,6353898,6354260'
seqType         = 1                                                 ## 0: Genomic coords and normal seq files 1: PacBio/Trinity - with lots of stuff in header
phase           = 24
head            = "Y"
coordsep        ='\t'                                               ## Seprator used in the file
phasedID        = 'N'                                               ## If given than phas ID is used to extract coords, if not than coords specified by user are used
chrcol          = 3                                                 ## As in excel format
startcol        = 4                                                 ## As in excel format
endcol          = 5                                                 ## As in excel format

phasiLenFilter  = 'Y'                                               ## 'Y' then tags of phase length will be written from the cluster | 'N' - All tags will be written
minAbun         = 1                                                 ## Minimum abundance of tag to be written
matchThres      = 0.95                                              ## Ratio of match required by the cluster to phased loci
startbuff       = 0                                                 ## While extracting sequence through coords2FASTA a buffer is added to start, 
                                                                    ## start position in phased ID has this buffer added, ## so minus this buffer 
                                                                    ## for better matching with real loci

libType         = 1                                                 ## 0: Lib_ids (4518) | 1: lib_code ('leaf_1')
excludeLibs = ['Daylily_leaf_1','Daylily_wbud','Daylily_01mm','Daylily_02mm','Daylily_03mm','Daylily_ant_0_r1','Daylily_ant_1_r1','Daylily_ant_2_r1','Daylily_ant_3_r1','Daylily_ant_4_r1','Daylily_fem_0_r1','Daylily_fem_1_r1','Daylily_fem_2_r1','Daylily_fem_3_r1','Daylily_fem_4_r1']
# excludeLibs = [4518,5006,5003,5004,5005,4519,4521,4523,4525,4527,4520,4522,4524,4526,4528]   ##  Used if fetchMax == 'Y' | Libs you wish to exclude, your tag position summary table should have libs ids and not lib codes

dataserver  = 'tarkan.dbi.udel.edu'
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

            if seqType == 0: 
                ## Normal genomic coordinates with integer chr_id
                get_chr_id = int(get_chr_id)
                chr_id = int(chr_id)
            else:
                ## Chromosomes are transcript names  i.e. strings
                pass

            
            if get_chr_id == chr_id:
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


    if fetchMax == 1: ### Fetch max phasi for each loci 
        cur= con.cursor()
        queryLibs,sumLibs = prepareQuery(excludeLibs,cur)
        
        outfile2 = "fetchMax.txt"
        fh_out2 = open(outfile2,'w')
        libsHead = queryLibs.split(",")
        fh_out2.write("loci\tmaxtag\tmaxTagAbun\ttotalPhasAbun\t#phaseTags\t%s\n" % ('\t'.join(x for x in libsHead)))
        # print(queryLibs) ## Libs whose abindance will be summed to give final abundance of tags
        abunList =   [] ### List to capture tags and abundance for each phased loci
        libAbunList = [] ## Lib-wise abundances of tag

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
                    # print(phasID,clustID,i[0],i[1],i[2])
                    fh_out.write('>%s_Clust%s_%s,%s,%s\n%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],tag,i[1],i[4]) )  ## phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun

            else:
                if int(i[1]) > minAbun:
                    # print(phasID,clustID,i[0],i[1],i[2])
                    fh_out.write('>%s_Clust%s_%s,%s,%sn%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],tag,i[1],i[4]) )  ## phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun
                    pass

            ## Get max abundance phasiRNA for specified phase
            if fetchMax == 1:
                atag,abun_sum,lib_abun = getAbundance(cur,tag,phase,queryLibs,sumLibs)
                libAbunList.append((lib_abun))

                if len(tag) == int(phase): ### Size specified in settings
                    abunList.append((atag,abun_sum))

        print("Elements in phasiList:%s" % (len(phasiList)))
        print("Elements in libAbunList:%s" % (len(libAbunList)))


        ## Process Fetch max results for this phased loci
        #################
        if fetchMax == 1:
            abunList_sort = sorted(abunList, key=operator.itemgetter(1),reverse=True) ## Sort list on abundances to get max abundant phasiRNA
            print("\nExample sorted values:%s" % (abunList_sort[0:10]))
            maxTag      = abunList_sort[0][0]   ## max abundant tag
            maxAbun     = abunList_sort[0][1]   ## Abundance of max abundant tag
            totalAbun   = sum(int(i[1]) for i in abunList_sort) ## Phased abundance

            ## Sum Lib-wise abundances for all tags
            libAbunSum = [0]*len(libAbunList[0]) ## This will hold sum of all tags, intialized for number of libraries
            for tag in libAbunList:
                # print(tag)
                libAbunSum  = [sum(x) for x in zip(libAbunSum,tag)]
            
            print("Tag:%s | maxPhasTag:%s | totalPhasAbun:%s" % (maxTag,maxAbun,totalAbun))
            print("Libwise Abundances",libAbunSum)

            ## Write - Loci, most abundant tag, most abundant tag abun,total phased abun, number of phased tags, and lib-wise abundances
            fh_out2.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (phasID,maxTag,maxAbun,totalAbun,str(len(abunList_sort)),'\t'.join(str(x) for x in libAbunSum)))
            abunList    = [] ## Empty before next phased loci
            libAbunList = [] ## Empty lib-iwse abundances of tag list before next entry
            # sys.exit()

    fh_out.close()
    fh_out2.close()

    return outfile,outfile2

def getAbundance(cur,tag,phase,queryLibs,sumLibs):
    '''Input is tag for each loci and out put is tag with maximmum abumdance and sum of phasiRNAs'''
    
    cur.execute("SELECT tag,%s FROM %s.%s where tag = '%s'" % (sumLibs,DB,tagPosTable,tag))### Convert intergenic to gene name so as to get strand
    info = cur.fetchall() ## Entries are redundant, one for every hit
    # print("--Tag summed abundance",info[0])
    atag,abun_sum = info[0]

    cur.execute("SELECT %s FROM %s.%s where tag = '%s'" % (queryLibs,DB,tagPosTable,tag)) ## Tag count from all libraries
    info2 = cur.fetchall() ## Entries are redundant, one for every hit
    print("--Lib-wise tag abundances:",info2[0])
    lib_abun = list(map(int, info2[0]))

    # sys.exit()

    return tag,abun_sum,lib_abun
            
def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### get column names
    columns = [] ## empty list
    cur.execute("describe %s.%s" % (DB,tagPosTable))
    tablefield = cur.fetchall()
    # print("\nTable fields:",tablefield)
    for i in tablefield:
        col = i[0]
        columns.append(col)

    libs = columns[8:-2] ## norm_sum and max_norm are last two columns
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
        
        finalLibs = selectLibs ## Norm_Sum and max_norm are not included

    else:
        finalLibs = libs ## ## Norm_Sum and max_norm are not included

    # print("finalLibs:%s" % (finalLibs))
    
    lib_col =",".join(str(x) for x in finalLibs)### Manually mentioned strating lib column - should work on all tag position summary tables
    
    print("Library Columns:",lib_col)
    # queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    sumLibs = "%s" % lib_col.replace(",","+")
    queryLibs = "%s" % lib_col.replace(",",",")
    print("\nThese are sumLibs:",sumLibs)
    print("\nThis is query Libs:",queryLibs)
    # sys.exit()

    return queryLibs,sumLibs

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

## v06 -> v07
## Added functionality to work on transcipts -but not tested on PacBio

## v07 -> v08
## Added functionality to get lib-wise abundances of all tags (not filtered on size like maxTag, nTags, and totalAbun)