#!/usr/local/bin/python3

###Script to extract phasiRNA from a cluster file for specific library
### Used phased loci name as in PARE validation results or coords in csv format to extract the phasiRNAs
### Requires 1) combined cluster fle matching p-value used to uniue phased loci 2) Non-redundant phased loci list 3)sRNA DB

import os,sys,difflib,time,operator
import mysql.connector as sql
from operator import itemgetter


##########Settings###################

clustfile       = "all.txt.score_p5e-07_sRNA_21_out.cluster"        ## Cluster file from phasing analysis, coul dbe for one lib or concatanated file for all libs
fetchMax        = 1                                                 ## Fetch the max abundance phasi from tag position summmary, write in separate file
fetchLibAbun    = 0                                                 ## 0: Fetch libwise abundances for tags with 'phase' len required fpr heatmaps and other comparision between loci | 1: All tags from phaster, required to see abundance od phasiRNA tag and tags against all tags in these loci. The latter would be sum of abundances from all libs computed manualy in excel
sRNADB          = "DAYLILY_priv_sRNA"                               ## if fetchMax == Y | Use this DB to get tag with max abundance

#### User specifed coords in input file
coordsfile      = "Final_PHAS_Loci_5e-07_ALL.csv"                   ## File with either phased loci ID 'Phas-100_w_7:6353898:6354260' or coords 'w,7,6353898,6354260'
seqType         = 1                                                 ## 0: Genomic coords and normal seq files 1: PacBio/Trinity - with lots of stuff in header
phase           = 21
head            = "Y"
coordsep        ='\t'                                               ## Seprator used in the file
phasedID        = 'N'                                               ## If given than phas ID is used to extract coords, if not than coords specified by user are used
chrcol          = 3                                                 ## As in excel format
startcol        = 4                                                 ## As in excel format
endcol          = 5                                                 ## As in excel format

phasiLenFilter  = 'N'                                               ## 'Y' then tags of phase length will be written from the cluster | 'N' - All tags will be written
minAbun         = 1                                                 ## Minimum abundance of tag to be written
matchThres      = 0.25                                              ## Ratio of match required by the cluster to phased loci | For transcripts, to capture biggest cluster like for mapping to IR, use a lower ratio so that longer transcript with smaller match ratio can be included
startbuff       = 0                                                 ## While extracting sequence through coords2FASTA a buffer is added to start, 
                                                                    ## start position in phased ID has this buffer added, ## so minus this buffer 
                                                                    ## for better matching with real loci

libType         = 0                                                 ## 0: Lib_ids (4518) | 1: lib_code ('leaf_1')
excludeLibs     = [4518,5006,5003,5004,5005,4519,4521,4523,4525,4527,4520,4522,4524,4526,4528]   ##  Used if fetchMax == 'Y' | Libs you wish to exclude, your tag position summary table should have libs ids and not lib codes

dataserver      = 'raichu.ddpsc.org'
DB              = "kakrana"

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
    resList2  = [] ## Store phasiRNAs from all clusters as (phas,[(phasiRNA),(PhasiRNA)],[extra info])
    
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
        finalMatchList  = []       ## Holds best cluster, from multiple libraries
        
        tempAllList    = [] ## Hold phasiRNAs from all matching clusters, of use for phased transcripts to capture allphasiRNAs, must be used with low matchThres
        for aclust in clusters[1:]:
            tempMatchList   = [] ## To hold results of current matching cluster
            aclust_splt     = aclust.split('\n')
            header          = aclust_splt[0].split()
            clust_id        = header[2]
            chr_id          = header[6].replace("chr","").replace("Chr","")
            start           = header[10]
            end             = int(header[12])+1 ##1 added because when opening a range using start and end, end number is not included in range
            value           = (list(range(int(str(start)),int(str(end)))))
            # print ('Cluster:', (value))

            if seqType == 0: 
                ## Normal genomic coordinates with integer chr_id
                get_chr_id = int(get_chr_id)
                chr_id = int(chr_id)
            else:
                ## Chromosomes are transcript names  i.e. strings
                pass

            
            if get_chr_id == chr_id:
                sm=difflib.SequenceMatcher(None,get_value,value) ## The rratio corresponds to larger loci i.e. match/length of longer nucleotides
                # print("Get Value:",get_value)
                # print("Current value",value)
                aratio  = round(sm.ratio(),2)
                # print("Ratio:%s" % (aratio))
                
                if round(sm.ratio(),2) >= matchThres:
                    ### Matched - phasiRNA from this cluster
                    print ('\nMatching cluster found:%s' % ''.join(header))
                    print("Allowed Ratio:%s | Current Ratio:%s" % (matchThres,aratio))
                    matchCount +=1
                    
                    phasiCyc = 0    ## Stores phasing cycles
                    phasiSig = 0    ## Stores total abundance of phase size sRNAs
                    otherSig = 0    ## Stores total abundance of other size sRNAs
                    kvalsL   = []   ## List to store kvalues for each phasiRNAs
                    
                    for i in aclust_splt[1:-1]:## Because header was the first entry of block and not required here, Last entry is always empty
                        # print ("Matched Cluster:\n",i)
                        phasient    = i.split('\t')
                        phasiname   = phasient[4].replace("|","_")
                        phasiseq    = phasient[5]
                        phasilen    = int(phasient[6])
                        phasiabun   = int(phasient[7])
                        phasihits   = int(phasient[10].split("=")[1])
                        phasipos    = int(phasient[3])
                        phasistrand = phasient[2].translate(str.maketrans("+-","wc"))
                        phasipval   = phasient[12]
                        phasikval   = int(phasient[9].split('=')[1])
                        # print(phasipos)
                        # sys.exit()

                        # print(phasiname,phasiabun,phasiseq,phasilen)
                        kvalsL.append(phasikval)
                        tempMatchList.append((phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos,phasistrand,phasipval))
                        tempAllList.append((phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos,phasistrand,phasipval)) ## Records all phasiRNAs from all clusters

                        if int(phasilen) == phase:
                            phasiCyc +=1
                            phasiSig += phasiabun
                        else:
                            otherSig += phasiabun
                    sizeRatio   = round(phasiSig/(phasiSig+otherSig),2) 
                    bestkval    = max(kvalsL) ## Best k-value achieved by this cluster
                    # print("Current Cycles:%s | Current sig. strength:%s" % (phasiCyc,phasiSig))
                    print("Current Cycles:%s | Current sig. strength:%s" % (bestkval,phasiSig))
                    tempMatchList.append((bestkval,phasiSig,phasID,clust_id,sizeRatio))

                    ## Decide the best and remove other from list ##############################
                    ############################################################################
                    if finalMatchList:
                        ## There exists a previosly matched cluster
                        exist_bestkval = finalMatchList[-1][0]
                        exist_phasiSig = finalMatchList[-1][1]
                        print("Existing Cycles:%s | Existing sig. strength:%s" % (exist_bestkval,exist_phasiSig))

                        if bestkval > exist_bestkval: ## New cluster has more cycles
                            del finalMatchList[0:]
                            finalMatchList = list(tempMatchList)
                            print("--- New cluster selected ---")

                        elif bestkval == exist_bestkval: ## Both have same cycles
                            if phasiSig > exist_phasiSig: ## New one has more total abundance of phased siRNAs
                                del finalMatchList[0:]
                                finalMatchList = list(tempMatchList)
                                print("--- New cluster selected ---")
                        
                        else: ## Existing/old one was long i.e. had more cycles
                            print("Earlier recorded cluster is retained")
                            pass

                    else: ## This is the first cluster
                        finalMatchList  = list(tempMatchList)
                        allphasiList    = list(tempMatchList) 

                    # print("\nFinal Match List:",finalMatchList)
            else:
                # print("No Match with this cluster")
                pass

        tempAllList.append((bestkval,phasiSig,phasID,clust_id,sizeRatio)) ## This list has phasiRNAs from all cluters but to keep the structure same as original resList, helpful while writing results, this info is added

        resList2.append((phasID,tempAllList))
        resList.append((phasID,finalMatchList)) ## Add best matched entry to the final list, there has to be one best matched entry per PHAS
        
        allMatchCount += matchCount ## Add the matched cluster for each entry
        if matchCount > 0 :
            uniqMatchCount+=1
        
    print("\nTotal phas loci: %s | Matched: %s" % (len(phasList),len(resList)))
    print ("SUMMARY: Phased loci in input file:%s | Loci match threshold: %s |Uniq matched cluster: %s | Total matched clusters found:%s" % (phasCount,matchThres,uniqMatchCount,allMatchCount))
    print("NOTE: If matched clusters more then phased loci that means same cluster was present in different libs\n")
    # print("NOTE: Don't forget to uniq the miRNAs")
    fh_in.close()

    return resList,resList2

def allphasiWriter(resList):
    '''
    This takes all the phasiRNAs on a transcript and return unique phasiRNAs file from all clusters - You need thos file for study of precursors
    '''

    print("\n\nFunction:allphasiWriter\n")

    outfile2     = outfile = clustfile+'thres_%s_allphasi.csv' % matchThres                       ## 2437.txt.score_p1e-07_sRNA_21_out.cluster
    fh_out      = open(outfile2,'w')

    ## Find entries with unique seq
    for ent in resList: ## entry corresponds to one phased loci
        # print("\nEntry",ent)
        phasID      = ent[0]
        phasCycles  = ent[1][-1][0]
        phasSig     = ent[1][-1][1]
        clustID     = ent[1][-1][2]
        sizeRatio   = ent[1][-1][4]
        phasiList   = ent[1][0:-1]
        phasiList_s = sorted(phasiList,key=itemgetter(5)) ## Sort on position so that results from multiple clusters are at-least sorted in postions
        # print("Sorted all phasiRNAs:", phasiList_s)
        # sys.exit()
        # print("%s | phasCycles:%s | phasSig:%s" % (phasID,phasCycles,phasSig))

        tagset = set() ## To record only the uniq tags, in case of lilium many a times same tag is documented twice
        for i in phasiList:
            print("-Phasi",i)  ## Phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos
            tag = i[2]
            
            if tag not in tagset: ## Avoid recordinf same tag agains as observed in casne on lilium
                tagset.add(tag)
                
                ## Write phasiRNAs 
                if phasiLenFilter == 'Y': ## If tags filter is ON
                    if (int(i[3]) == int(phase)) and (int(i[1]) >= minAbun): ### Size specified in settings
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)

                else:
                    if int(i[1]) > minAbun:
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)
                        pass


    fh_out.close()
    print("\n\nExiting :allphasiWriter\n")


    return outfile2

def writer(resList,con):
    '''
    write the results
    '''

    print("\n\nFunction:Writer\n")

    outfile = clustfile+'thres_%s_phasi.csv' % matchThres                       ## 2437.txt.score_p1e-07_sRNA_21_out.cluster
    fh_out  = open(outfile,'w')


    if fetchMax == 1: ### Fetch max phasi for each loci 
        cur= con.cursor()
        queryLibs,sumLibs,finalLibs = prepareQuery(excludeLibs,cur)
        
        outfile2    = "phassummary.txt"
        fh_out2     = open(outfile2,'w')
        libsHead    = queryLibs.split(",")
        fh_out2.write("loci\tbest-K-value\tsizeRatio\t%s\ttotalPhasAbun\tmaxtag\tmaxTagAbun\tmaxtag2\tmaxTagAbun2\n" % ('\t'.join(x for x in libsHead)))
        # print(queryLibs) ## Libs whose abindance will be summed to give final abundance of tags
        abunList    = [] ### List to capture tags and abundance for each phased loci
        libAbunList = [] ## Lib-wise abundances of tag

    for ent in resList: ## entry corresponds to one phased loci
        print("\nEntry",ent)
        phasID      = ent[0]
        phasCycles  = ent[1][-1][0]
        phasSig     = ent[1][-1][1]
        clustID     = ent[1][-1][2]
        sizeRatio   = ent[1][-1][4]
        phasiList   = ent[1][0:-1]
        # print("%s | phasCycles:%s | phasSig:%s" % (phasID,phasCycles,phasSig))

        tagset = set() ## To record onlu uniq tags, in case of lilium many a times same tag is documented twice
        for i in phasiList:
            print("-Phasi",i)  ## Phasiname,phasiabun,phasiseq,phasilen,phasihits,phasipos
            tag = i[2]
            
            if tag not in tagset: ## Avoid recording same tag agains as observed in casne on lilium
                tagset.add(tag)
                
                ## Write phasiRNAs 
                if phasiLenFilter == 'Y': ## If tags filter is ON
                    if (int(i[3]) == int(phase)) and (int(i[1]) >= minAbun): ### Size specified in settings
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s,%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)

                else:
                    if int(i[1]) > minAbun:
                        # print(phasID,clustID,i[0],i[1],i[2])
                        fh_out.write('>%s_Clust%s_%s,%s,%s,%s,%s,%s,%s\n%s%s,%s,%s,%s,%s,%s\n' % (phasID,clustID,i[0],i[1],i[4],i[5],i[3],i[6],i[7],tag,i[1],i[4],i[5],i[3],i[6],i[7]) )  ## phasID,clust_id,phasiname,phasiabun,phasihits,phasiseq,phasiabun,phasihits,phasipos,phasistrand,phasipval
                        tagset.add(tag)
                        pass


                ## Get max abundance phasiRNA for specified phase
                if fetchMax == 1:

                    ## Get lib-wise abundaces mode
                    if fetchLibAbun == 0:
                        if len(tag) == int(phase):
                            atag,abun_sum,lib_abun = getAbundance(cur,tag,finalLibs)
                            libAbunList.append((lib_abun))

                    elif fetchLibAbun == 1: ## All the tags
                        atag,abun_sum,lib_abun = getAbundance(cur,tag,finalLibs)
                        libAbunList.append((lib_abun))
                    
                    else:
                        print("Libwise abundances won't be fetched")
                        pass

                    ## Tag specific abundances for fetching most abundant tag
                    if len(tag) == int(phase): ### Size specified in settings
                            abunList.append((atag,abun_sum))
            else:
                print("Tag recorded once already#####################################\n")
                # sys.exit()
                pass


        print("Elements in phasiList:%s" % (len(phasiList)))
        print("Elements in libAbunList:%s" % (len(libAbunList)))


        ## Process Fetch max results for this phased loci
        #################
        if fetchMax == 1:
            abunList_sort = sorted(abunList, key=operator.itemgetter(1),reverse=True) ## Sort list on abundances to get max abundant phasiRNA
            print("\nExample sorted values:%s" % (abunList_sort[0:10]))
            maxTag      = abunList_sort[0][0]   ## Max abundant phasiRNA sequence
            maxAbun     = abunList_sort[0][1]   ## Abundance of max abundant phasiRNA
            if len(abunList_sort) > 1: ## In one case no second tag was found
                maxTag2     = abunList_sort[1][0]   ## Max abundant phasiRNA sequence
                maxAbun2    = abunList_sort[1][1]   ## Abundance of max abundant phasiRNA
            else:
                maxTag2     = "na"
                maxAbun2    = "0"
                print("The predeiction is of really low quality - Just one tag of phase size found")
                time.sleep(2)
                # sys.exit()
            totalAbun   = sum(int(i[1]) for i in abunList_sort) ## Phased abundance


            ## Sum Lib-wise abundances for all tags
            libAbunSum = [0]*len(libAbunList[0]) ## This will hold sum of all tags, intialized for number of libraries
            for tag in libAbunList:
                # print(tag)
                libAbunSum  = [sum(x) for x in zip(libAbunSum,tag)]
            
            print("Tag:%s | maxPhasTag:%s | totalPhasAbun:%s" % (maxTag,maxAbun,totalAbun))
            print("Libwise Abundances",libAbunSum)

            ## Write - Loci, most abundant tag, most abundant tag abun,total phased abun, number of phased tags, and lib-wise abundances
            fh_out2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (phasID,str(phasCycles),sizeRatio,'\t'.join(str(x) for x in libAbunSum),totalAbun,maxTag,maxAbun,maxTag2,maxAbun2))
            abunList    = [] ## Empty before next phased loci
            libAbunList = [] ## Empty lib-iwse abundances of tag list before next entry
            # sys.exit()

    fh_out.close()
    fh_out2.close()

    return outfile,outfile2

def getAbundance(cur,tag,finalLibs):
    '''Input is tag for each loci and out put is tag with maximmum abumdance and sum of phasiRNAs - 
    rewritten in v1.0 for fetching tags from run master'''

    lib_abun = [] ## list to hold lib-wise abudnances
    
    for alib in finalLibs:
        # print("Lib:",alib)
    
        cur.execute("SELECT tag,norm FROM %s.run_master where tag = '%s' and lib_id = %s" % (sRNADB,tag,alib))### Convert intergenic to gene name so as to get strand
        info = cur.fetchall() ## Entries are redundant, one for every hit
        # print("Query fetched", info)

        if info:
            atag,norm_abun = info[0]
            lib_abun.append(norm_abun)
            # print("--Tag abundance:%s for lib:%s"% (tag,norm_abun))
        else:
            norm_abun = 0
            lib_abun.append(norm_abun)
            # print("--Tag abundance:%s for lib:%s"% (tag,norm_abun))


    abun_sum = sum(lib_abun)
    print("--Lib-wise abundances",lib_abun)
    print("--Sum of abundances:%s\n" % (abun_sum))

    # sys.exit()

    return tag,abun_sum,lib_abun
            
def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### Get lib names
    columns = [] ## empty list
    cur.execute("SELECT DISTINCT(lib_id) FROM %s.run_master" % (sRNADB))
    info = cur.fetchall()
    libs = [x[0] for x in info]

    print("\nLibs:",libs)

    if excludeLibs:
        print("\n\nLibs specifed in excludeLibs %s will be skipped\n\n" % (excludeLibs))
        selectLibs      = [] ## Columns excluding the unwanted libraries
        excludeLibs_s   = [str(i) for i in excludeLibs] ## Converting all entries in exclude list to string for matching below
        
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

            ### Filter libraries
            if str(i) not in excludeLibs_s: ## Tested OK
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
    
    print("\nLibrary Columns:",lib_col)
    # queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    sumLibs = "%s" % lib_col.replace(",","+")
    queryLibs = "%s" % lib_col.replace(",",",")
    print("\nThese are sumLibs:",sumLibs)
    print("\nThis is query Libs:",queryLibs)
    # sys.exit()

    return queryLibs,sumLibs,finalLibs

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
####################################
def main():
    phasList = reader(coordsfile)
    time.sleep(1)

    resList,resList2 = getClust(clustfile,phasList)
    if seqType == 1:
        allphasiFile = allphasiWriter(resList2)
    else:
        print("File with all phas will not be generated for this 'seqtype'")
        pass


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
## Fixed library exclusion - This needs to be propoagated to getCoords ABundances script

## v08 -> v09
## Added "fetchLibAbun" switch that fetches abudnances for all libs, either for tags that have len equal to phase or all tags irrespective of size

## v09 -> v1.0
## Elimnate need of tag position summary table requirement which is kind of unrealiable in case of non-genome species

## v1.0 -> v1.01
## Added functionality to output second max tag seq and abundance, so that it can be used to compute ratio of top two tags 
## against PHAS abundance or all sRNA abundance in fetchLibAbun = 1

## v1.01-> v1.02
## Added functionality to use only unique tags for a loci, as in case of lilium same tag map to two different arms and recorded twice

## v1.02 -> v1.03
## In case of really low p-value predictions, sometimes there is just one tag of phassize, in such cases fetchMax does not have second best tag and gave error. Now if second max tag is not available results will have it as "na" and 0 abundance

## v1.03 -> v1.04
## Added extra columns to phasiRNA fasta/csv file - position, p-val and strand

## v1.04 -> v1.05
## Added an additional module allPhasiFile - that generates an addition phasiRNA file that includes all phasiRNAs for any matching clsuter. FOr this low matchThres should be used ~ 0.25-0.50. This is useful to capture all phasiRNAs for transcripts to get seconadry structure

## v1.05 -> v1.06
## Fixed bug added in 1.05 that has a comma before phasiRNA seq in fasta file. This comma caused problem in automatic secondary structure detection as phasiRNA would show empty

## v1.06 -> v1.07
## Added best kvalue for every cluster in phas summary
## Added a size-speciifc ratio of abudnances for each cluster i.e. abundance of all 21-/24-phasiRNAs comapred to total
## Organized the phassummary file for clarity 