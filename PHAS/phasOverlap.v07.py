#!/usr/local/bin/python3

## This script takes a PHAS file (bed format) and a GTF file of assembly/genome featires and find overlapping transcripts to PHAS, plus computes degree of overlap with reference to PHAS, also if provided a de novo PHAS file (from PHASredunndant) with precursors/transcipts then it adds a validation flag "p" suggesting the that transcript overlapping a PHAS was indeed detected as PHAS  - You can use this identify precursors from assembly 

## Written by Atul Kakrana - kakrana@udel.edu

import os,sys,time,sqlite3,operator
import mysql.connector as sql
from collections import Counter

### User inputs #####
aStep           = 1                                                             ## 1: Generate an overlapping file   
phasFile        = "PNAS_24PHAS.bed"                                             ## Phased file which needs to be checke for overlap
delim           = "\t"                                                          ## Delimiter for summary file
head            = 'N'                                                           ## Header is summary file: 'Y' else: 'N'
phasName        = 4
chrpos          = 1
startpos        = 2
endpos          = 3
# strandpos   = 5

exonOverlapPerc     = 0.2                                                       ## Minimum percentage of PHAS covered by exon. Keep 0 for finding overlapping genes and 0.20 or above for finding precursors
exonsOverlapCutoff  = 100                                                       ## Values of minimum overlap in nts. This is required if PHAS overlap is checked for precursors and not overlap with any regio of genes. FOr precursors one will expect an overlap of 4 phases or more
featureFile     = "Trinity.collapsed.gtf"                                       ## Feature file in GTF format, with whome overlap needs to be checked
featureFile2    = "v3.genome.merged.gtf"                                        ## Second feature file if overlap needs to be checked with two lists, but independently
mode            = 3                                                             ## 1: PASA feature (GTF) file | 2: Trinity | 3: Rocket feature (GTF) file or generic gtf file with missing "transcript" entry | 4: Both rocket and trinity feature files | 5: PacBio GTF | 6: Comprehensive transciptome (from PASA) - NOTE: New mode should be registered in overlapchecker function

phasedFeature   = "Final_PHAS_Loci_5e-07_ALL.txt"                               ## [Optional, give an empty file] Phasredundant file of transcripts from feature file, to label overapping transcripts as true phased. Provide concataned file if both rocket and trinity GTF file used in mode3

strandbiasStep  = 0                                                             ## 0: No | 1: yes
mainLen         = [24,24]                                                       ## For stand Ratio
noiseLen        = [24,23,22]                                                    ## For strand ratio
includeLibs     = ["lib_1739"]        ## Libraries to compute strand ratio, pre-meiotic:  | Meiotic:"ant1_07mm_r1","Asp_1_ant_bud"

### Developer inputs ###
gtfScore        = 100
dataServer      = "tarkan.ddpsc.org"
msDB            = "kakrana"
tagPosTable     = "MAIZE_pub_sRNA_WGSTagPosSummNorm"
startbuff       = 0                                                         ## This will be used to reduce are back to actual Loci - 0 for PARE
endbuff         = 0

##### Functions #############
#############################
def overlapChecker(phasList,gtfList,phasedFeature,includeLibs):
    '''Checks for overlap between genomic PHAS and transcipts frpm GTF file'''
    
    print("\nFunction: overlapChecker")
    
    ## Prepare output file ####################################
    ###########################################################
    outFile = "overlapResults.txt"
    fh_out  = open(outFile,'w')
    fh_out.write("PHAS\tphasChr\tphasStart\tphasEnd\tphasStrand\tOverlapTrans\toverlapNucl\toverlapPerc\tstrand\ttransLen\tdataFlag\tphasStatus\tphasiRatio\tstrandRatio\tpStrand\n")
    
    ## Prepare DB and featureTables ##########################
    ##########################################################
    DB = 'tempdb'
    try:
        os.remove(DB)
    except OSError:
        pass
    ## Prepare features table - This might include entries from two GTF files (mode-2), identified by the flags
    conn            = sqlite3.connect(DB)
    featureTable    = tableMaker(gtfList,conn) ## TURN ON ****
    
    # # Test Query
    # cur = conn.cursor()
    # cur.execute("SELECT * FROM %s where chr = 1 AND strand = 'w' limit 10" % (featureTable))
    # test = cur.fetchall()
    # print("Test results:",test)
    # # sys.exit()

    ## Flags i.e. data to query - This allows combining results from Trinity and rocket to one sheet
    flags     = [] ## Datatype to query
    if mode     == 1 or mode == 6:
        flags.append(('P'))
    elif mode   == 2:
        flags.append(('T'))
    elif mode   == 3:
        flags.append(('R'))
    elif mode   == 4:
        flags.append(('R'))
        flags.append(('T'))
    elif mode == 5:
        flags.append(('PB'))
    else:
        print("Please check the mode selected for analysis")
        pass

    ## Find transcript overlapping PHAS ########################
    ############################################################
    print("\nAnalysis to identify overlapping transcripts initialized")
    featurePhasList,featurePhasName = featurePhasReader(phasedFeature) ## PHAS transcripts identified from rocket, trinity, pasa transcripts - These will be used to label the overlapping transcripts as phased
    
    nonOverlapList  = [] ## Temp list so that these non-overlapping ones are just added once,as there is no uniq transcipt attached to them
    overlapNameList = [] ## Just the naemes of phased entries that has overlapping results
    for aflag in flags: ##FOr provided datatypes
        print("\nFlag:",aflag)
        
        for ent in phasList: ## FOr every PHAS
            print("\n###################Entry:",ent)
            aname,achr,astart,aend,astrand = ent
            transList = overlapTrans(ent,conn,featureTable,aflag) ## TURN ON ****

            ## Check overlap with exons of every overlapping transcript
            for trans in transList: ## FOr every overlapping trasc
                atrans,tstrand,tlen = trans
                print("\n-Computing overlap for transcript:%s" % (atrans))
                exonsOverlap = overlapExons(ent,conn,featureTable,atrans,aflag)
                aperc = round((exonsOverlap/(aend-astart)),2)

                ## Get PHAS info of overlapping trans ###
                #########################################
                if atrans in featurePhasName:
                    print("%s overlapping transcript is phased" % (atrans))
                    pStatus = 'p'
                    # sys.exit() ## Works 
                else:
                    pStatus = 'no'

                ## Get strand based on sRNA ############
                #######################################
                if strandbiasStep == 1:
                    mscon   = ConnectToDB(dataServer) ## MySQL Connection
                    mscur   = mscon.cursor()          ## MySQL Cursor
                    phasiRatio,strandRatio,pstrand = srnaStrand(ent,msDB,tagPosTable,includeLibs,mscur)
                else:
                    phasiRatio  = 0
                    strandRatio = 0
                    pstrand     = "n"


                ## Filters and output results############
                ########################################
                if aperc >= exonOverlapPerc and exonsOverlap >= exonsOverlapCutoff:
                    overlapNameList.append(aname)
                    fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('\t'.join(str(x) for x in ent),atrans,str(exonsOverlap),str(aperc),astrand,str(tlen),aflag,pStatus,str(phasiRatio),str(strandRatio),pstrand))
                else:
                    nonOverlapList.append((ent,phasiRatio,strandRatio,pstrand)) ## This in current form wont catch those entries for which there is no overlap with transcripts i.e. PHAS with no overlap to any gene. The filter is designed to include those cases for which there is an overlap to transcripts but not with exons
                    pass

    #### Write results for which no overlapping transcript was found
    writtenList = [] ## Temp list to check the result for entry has been written to avoid duplicate for PHAS which has no results
    for i in nonOverlapList:
        ent,phasiRatio,strandRatio,pstrand = i
        phasName = ent[0]
        if phasName not in overlapNameList and phasName not in writtenList:
            fh_out.write("%s\tna\tna\tna\tna\tna\tna\tna\t%s\t%s\t%s\n" % ('\t'.join(str(x) for x in ent),str(phasiRatio),str(strandRatio),pstrand))
            writtenList.append(phasName)
     

    # sys.exit()
    fh_out.close()

    return outFile

def phasParser(phasFile):
    '''parses input PHAS file'''

    print("\nFunction: phasParser")

    fh_in = open(phasFile,'r')
    if head == 'Y':
        fh_in.readline()
    phasRead = fh_in.readlines()

    phasList = []
    for i in phasRead:
        ent = i.split(delim)
        aname   = ent[phasName-1]
        achr    = int(ent[chrpos-1].replace('chr','').replace("omosome",""))
        astart  = int(ent[startpos-1])
        aend    = int(ent[endpos-1])
        # astrand = ent[strandpos-1]
        astrand = "na"
        print("Phas Name %s | chr:%s | start:%s | end:%s | astrand:%s" % (aname,achr,astart,aend,astrand))
        phasList.append(((aname,achr,astart,aend,astrand)))

    print("A list of query features prepared with entries:%s" % len(phasList))

    print ("Exiting function - phasParser\n")
    return phasList

def gtfParser(afeatureFile):

    '''Parses gtf file into featurename,chr,start,end,strand,feature'''

    print("\nFunction: gtfParser")
    
    with open(afeatureFile) as fh_in:
        lines = (line.rstrip() for line in fh_in) 
        gtfRead = list(line for line in lines if line) # Non-blank lines in a list
    fh_in.close()

    gtfList = [] ## List to hold parsed GTF entries

    for i in gtfRead:
        # print(i)
        ent = i.split("\t")
        # print("\nEnt",ent)
        if ent[2] == "transcript" or ent[2] == "exon":
            # print(ent)
            gchr    = ent[0].replace("chr","")
            gtype   = ent[2]
            gstart  = ent[3]
            gend    = ent[4]
            gstrand = ent[6].translate(str.maketrans("+-","wc"))
            aflag   = "P" ## PASA 
            info    = ent[8].strip("\n").split(";")
            # print(info,len(info))
            if len(info) == 3: ## With last one empty 
                ## Protein coding gene with a version number
                gid     = info[0].split()[1].replace('"','') ## Gene ID
                tid     = info[1].split()[1].replace('"','') ## Transcript ID
                # print(gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))
            # if mode == 6 and len(info) == 4: ## With last one empty 
            #     ## Protein coding gene with a version number
            #     gid     = info[0].split()[1].replace('"','') ## Gene ID
            #     tid     = info[1].split()[1].replace('"','') ## Transcript ID
            #     # print(gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
            #     gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))
            else:
                print("This entry has more than expected info")
                print(info)
                print("Debug!!")
                sys.exit()

        else:
            print("We don't need this info for current script")
            pass

    if len(gtfList) == 0:
        print("Check if the feature used for extraction i.e. gene, transcript, exon is correct")
        print("Debug!!")
        sys.exit()

    print("First 10 entries of gtfList list: %s" % (gtfList[:10]))
    print("Total entries fetched from GTF file:%s" % (str(len(gtfList))))

    print ("Exiting function - gtfParser\n")
    time.sleep(1)
    
    return gtfList

def gtfParser2(afeatureFile):

    '''This function parses Trinity and Rocket GTF file
    to give a trascript entry and entry for all exons - Basically this is the parser for gffread geenrated files'''

    print("\nFunction: gtfParser")
    
    with open(afeatureFile) as fh_in:
        lines = (line.rstrip() for line in fh_in) 
        gtfRead = list(line for line in lines if line) # Non-blank lines in a list
    fh_in.close()

    gtfList     = [] ## List to hold parsed GTF entries
    tempName    = [] ## Stores current trans name
    tempCoords  = [] ## Temp coords
    for i in gtfRead:
        # print(i)
        ent = i.split("\t")
        # print("\nEnt",ent)
        gScore  = ent[5] ## Score provided for mapping accuracy in Trinity GTF from GMAP 0 to 100
        gtype   = ent[2]
        if gtype == "exon" and (gScore == '.' or float(gScore) > 90.0):
            # print("\nExon:",ent)
            gchr    = ent[0].replace("chr","").replace("omosome","")
            gstart  = int(ent[3])
            gend    = int(ent[4])
            gstrand = ent[6].translate(str.maketrans("+-","wc"))
            info    = ent[8].strip("\n").split(";")
            # print(info)
            
            ## Parse the info and add the exon entries  #################
            #############################################################

            if mode == 2 or (mode == 4 and len(info)==4): ## This is  trinity GTF info
                ## Protein coding gene with a version number
                tid     = info[0].split()[1].replace('"','').split(".")[0] ## Transcript ID
                gid     = info[2].split()[1].replace('"','').rpartition("_")[0] ## Gene ID
                aflag   = 'T' ## Trinity
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))

            elif mode == 3 or (mode ==4 and len(info) >= 7): ## This is rocket GTF info, with and w/o p_id "P1" in the end
                tid     = info[0].split()[1].replace('"','') ## Transcript ID
                gid     = info[1].split()[1].replace('"','') ## Gene ID
                aflag   = 'R' ## Rocket
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))
            
            elif mode == 5:
                tid     = info[0].split()[1].replace('"','').split(".")[0] ## Transcript ID
                gid     = info[1].split()[1].replace('"','').split(".")[0] ## Gene ID
                aflag   = 'PB' ## PacBio
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))

            elif mode ==6:
                gid     = info[0].split()[1].replace('"','') ## Gene ID
                tid     = info[1].split()[1].replace('"','') ## Transcript ID
                aflag   = 'P' ## Pasa
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))


            else:
                print("-This entry has more than expected info")
                print("-Info block",info)
                print("-Info block len:%s" % len(info))
                print("-Debug!!")
                sys.exit()

            ## Check trascript change #########
            ###################################

            ## Check if transcript entry needs to be added i.e. all the exons from previous one has been captured
            if not tempName or tempName[0] == tid:
                tempName.append(tid)
                tempCoords.append(gstart)
                tempCoords.append(gend)
                tempStrand  = gstrand   ## strand of last transcript
                tempgid     = gid       ## gene name of last transcript
                tempChr     = gchr      ## chr of last transcript
                # print(tempCoords)
            elif tempName[0] != tid: ## Check if the new transcript is being read. if yes then add transcript entry for old exons using tempName and tempCoords
                # print("-New transcript read - summarizing transcript entry for earlier one")
                tstart      = min(tempCoords)
                tend        = max(tempCoords)
                ttid         = tempName[0]
                ttype       = 'transcript'
                # print('-',tempgid,ttid,tempChr,tstart,tend,tempStrand,ttype,aflag)
                gtfList.append((tempgid,ttid,tempChr,tstart,tend,tempStrand,ttype,aflag))
                # sys.exit()
                
                ## Empty lists and fill with current exon from new transcript
                tempName    = [] ## Empty trans name
                tempName.append(tid)
                tempCoords  = []
                tempCoords.append(gstart)
                tempCoords.append(gend) ## Empty trans coords
                tempStrand  = gstrand   ## strand of last transcript
                tempgid     = gid       ## gene name of last transcript
                tempChr     = gchr      ## chr of last transcript
                # sys.exit()
            else:
                print("-Unforseen scenario encountered")
                sys.exit()

        else:
            # print("We don't need this info for current script") ## CDS or GENE
            pass

    if len(gtfList) == 0:
        print("Check if the feature used for extraction i.e. gene, transcript, exon is correct")
        print("Debug!!")
        sys.exit()
    else:
        print("First 10 entries of gtfList list: %s" % (gtfList[:10]))

    print("Total entries fetched from GTF file:%s" % (str(len(gtfList))))

    print ("Exiting function - gtfParser2\n")
    # sys.exit()
    time.sleep(1)
    
    return gtfList

def tableMaker(alist,conn):
    '''makes SQLlite table for query'''

    print("-Preparing table with features")

    cur = conn.cursor()

    ## Make table
    featuretable = "tempTable"
    cur.execute('''DROP TABLE IF EXISTS %s''' % (featuretable)) ### Drop Old table - while testing    
    conn.commit()
    
    try:
        cur.execute('''CREATE TABLE %s (gene varchar(255),trans varchar(255),chr integer, start integer, end integer, strand varchar(10), type varchar(255),flag varchar(255))''' % (featuretable))
        conn.commit()
        for i in alist:
            # print(i)
            gid,tid,gchr,gstart,gend,gstrand,gtype,aflag = i
            # print ('Gene:',gid,'Trans:',tid,'| Chr:',gchr,'| Start:',gstart,'| End:',gend, '| Strand:',gstrand,'| Type:',gtype,'| Flag:',aflag)
            cur.execute("INSERT INTO %s VALUES ('%s','%s',%d,%d,%d,'%s','%s','%s')" % (featuretable,str(gid),str(tid),int(gchr),int(gstart),int(gend),str(gstrand),str(gtype),str(aflag)))
            # featureList.append((str(aname),int(achr),int(astart),int(aend),str(astrand)))
            #conn.commit()
    
    except sqlite3.Error:
        print('ERROR:',Error)
        sys.exit()

    print("\n-Feature table made sucessfully")

    return featuretable

def overlapTrans(ent,conn,featureTable,aflag):
    '''This function returns a list of transcipts that overlaps with a PHAS'''
    transList = [] ## Store the transcipts that overlap

    aname,achr,astart,aend,astrand = ent
    cur = conn.cursor()

    print("-FLag being queried:%s" % (aflag))
    print("-Query entry:",ent)
    
    ## trans flanking phas or enclaved in phas
    cur.execute("SELECT * FROM %s where chr = %s AND flag = '%s' AND ((end between %s and %s) or (start between %s and %s)) AND type = 'transcript'" % (featureTable,achr,aflag,astart,aend,astart,aend))
    flankTrans = cur.fetchall()

    ## PHAS enclaved in trans
    cur.execute("SELECT * FROM %s where chr = %s AND flag = '%s' AND (%s between start and end) AND (%s between start and end) AND type = 'transcript'" % (featureTable,achr,aflag,astart,aend))
    bigTrans = cur.fetchall()

    print("flanking Trans:",flankTrans)
    print("Long Trans:",bigTrans)

    ## Combine both lists and report overlapping trans
    allTrans = flankTrans + bigTrans
    for i in allTrans:
        # print(i)
        atrans  = i[1]
        astrand = i[5] 
        alen    = i[4]-i[3]
        transList.append((atrans,astrand,alen))

    print("-Number of overlapping trans:%s" % (len(transList)))
    print("-transList",transList)

    return transList

def overlapExons(ent,conn,featureTable,atrans,aflag):
    '''This function will compute overlaps wit exons of overlapping transcripts'''

    aname,achr,astart,aend,astrand = ent
    
    cur = conn.cursor()
    cur.execute("SELECT * FROM %s where flag = '%s' AND trans = '%s' AND type = 'exon'" % (featureTable,aflag,atrans))
    exons = cur.fetchall()
    print("-These are the exons",exons)

    ## compute overlap with exons of this transcript
    exonsOverlap = 0
    for aexon in exons:
        print("-Checking exons:",aexon)
        xoverlap = 0
        xstart  = aexon[3]
        xend    = aexon[4]

        ## Exon is enclaved
        if astart <= xstart and aend >= xend:
            print("-exon enclaved")
            xoverlap = xend - xstart + 1

        ## PHAS is enclaved
        elif xstart <= astart and xend >= aend:
            print("-PHAS enclaved")
            xoverlap = aend - astart + 1
        
        ## Overlap at 5' of PHAS, use phas start as reference
        elif xstart <= astart and xend >= astart:
            print("- 5' flank")
            xoverlap = xend - astart + 1
        
        ## Overlap at 3' of PHAS, use PHAS end as reference
        elif xstart <= aend and xend >= aend:
            print("- 3' flank")
            xoverlap = aend - xstart + 1

        else:
            print("-Unexpected overlap encountered - Transcript Overlaps but not Exons?")
            pass
            # sys.exit()

        exonsOverlap+=xoverlap

    print("-Overlap found:%s" % exonsOverlap)
    return exonsOverlap

def featurePhasReader(phasedFeature):
    '''Reads de novo phased file and reports name, p-val,start,end coords as a list, plus names in a separate list - This is done to validate the overlap'''

    print("\nFunction: featurePhasReader")

    fh_in = open(phasedFeature,'r')
    fh_in.readline() ## PhasRedundant file has header
    phasRead = fh_in.readlines()
    fh_in.close()

    phasList = [] ## Full PHAS details
    phasTransList = [] ## Just the names of phased transcripts - This transcripts are from featureFiles i.e. gtf
    for i in phasRead:
        print(i.strip("\n").split('\t'))
        aname,apval,atrans,astart,aend,trash,alib= i.strip("\n").split('\t')
        
        if atrans.startswith("XLOC"): ## The GTF extractor v02 and above merges genename and transcript name, plus rc if no strand info available - These are cached and fixed to just have transcript name so that it matched with GTF file
            btrans = atrans.replace("_rc","").split("_",2)[-1]
            # sys.exit()
        else:
            btrans = atrans
            pass

        # print("phased Feture Trans:%s" % (btrans))
        
        phasList.append((aname,apval,btrans,astart,aend,alib))
        phasTransList.append(atrans)

    print("-Snippet of feature PHAS full list",phasList[1:10])
    print("-Snippet of feature PHAS names list",phasTransList[1:10])
    print("-Entries in full info list:%s | Names list:%s" % (len(phasList),len(phasTransList)))

    print ("Exiting function - phasReader\n")
    # sys.exit()

    return phasList,phasTransList

def resReader(resFile):
    '''Reads results file generated from step one, used as a checkpoint to avoid redoing
    the overlap analysis again'''


    print("\nFunction: resReader")
    fh_in = open(resFile,'r')
    resRead = fh_in.readlines()
    fh_in.close()

    overlapResList = []
    for i in resRead:
        aphas,pChr,pStart,pEnd.pStrand,tOverlap,tOverNucl,tOverRatio,tStrand,tlen,flag,tPhasStatus,phasiRatio,strandRatio,pstrand = i.strip("\n").split("\t")
        overlapResList.append((aphas,pChr,pStart,pEnd.pStrand,tOverlap,tOverNucl,tOverRatio,tStrand,tlen,flag,tPhasStatus,phasiRatio,strandRatio,pstrand))

    print("-Results cached from 'overlapResults' file")
    print("-Entries in full info list:%s" % (len(overlapResList)))

    print ("Exiting function - resReader\n")

    ## PHAS\tphasChr\tphasStart\tphasEnd\tphasStrand\tOverlapTrans\toverlapNucl\toverlapPerc\tstrand\ttransLen\tdataFlag\tphasStatus\n"
    
    return overlapResList

def srnaStrand(ent,DB,tagPosTable,includeLibs,cur):
    ''' Computes strand ratio of sRNA abundnaces and predictes a strand from which sRNAs must have been produced'''

    name,chr_id,start,end,strand = ent
    # print ('\nEnt',ent)
    start   = ent[2] ##Buffer of 100bp or as specified
    start2  = int(start)+startbuff ### Reverse the effect on co-ordinate made during fetching fasta
    end     = ent[3]
    end2    = int(end)-endbuff
    length  = int(end2)-int(start2) ###Extra 1 to include last position, buffer of 100 bp from start and end
    chr_id  = ent[1]
    strand  = ent[4]
    # print ('This is chr_id:%s and Start:%s and End:%s | Length of loci: %s bp ' % (chr_id,str(start2),int(end2),length))

    ## Prepare queries for both main size and noise ###
    ###################################################
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

    ### Prepare and perform query ############
    ##########################################
    lib_col =",".join(str(x) for x in includeLibs)### Manually mentioned strating lib column - should work on all tag position summary tables
    # print("Library Columns:",lib_col)
    queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    # print("Query:",queryLibs)

    ## Get main size
    cur.execute("select %s from %s.%s where chr_id = %s and strand = 'w' and (%s) and (position between %s and %s)" % (queryLibs,DB,tagPosTable,chr_id,mainSizeQuery,str(start2),str(end2)))####
    temp = cur.fetchall()
    main_abun_w = list(map(int, temp[0]))

    cur.execute("select %s from %s.%s where chr_id = %s and strand = 'c' and (%s) and (position between %s and %s)" % (queryLibs,DB,tagPosTable,chr_id,mainSizeQuery,str(start2),str(end2)))####
    temp = cur.fetchall()
    main_abun_c = list(map(int, temp[0]))

    print(main_abun_w)
    print(main_abun_c)
    
    ## Get noise size
    cur.execute("select %s from %s.%s where chr_id = %s and strand = 'w' and (%s) and (position between %s and %s)" % (queryLibs,DB,tagPosTable,chr_id,noiseSizeQuery,str(start2),str(end2)))####
    temp2 = cur.fetchall()
    noise_abun_w = list(map(int, temp2[0]))
    
    cur.execute("select %s from %s.%s where chr_id = %s and strand = 'c' and (%s) and (position between %s and %s)" % (queryLibs,DB,tagPosTable,chr_id,noiseSizeQuery,str(start2),str(end2)))####
    temp2 = cur.fetchall()
    noise_abun_c = list(map(int, temp2[0]))
    
    print(noise_abun_w)
    print(noise_abun_c)

    mainSum_w      = sum(main_abun_w) ## in comma separated format
    mainSum_c      = sum(main_abun_c) ## in comma separated format
    
    noiseSum_w     = sum(noise_abun_w) ## for different size  i.e 22,24 or what ever defined in settings
    noiseSum_c     = sum(noise_abun_c) ## for different size  i.e 22,24 or what ever defined in settings
    print("Sum for main size on 'w':%s | Sum for main size on 'c':%s" % (mainSum_w,mainSum_c))
    print("Sum for noise size on 'w':%s | Sum for noise size on 'c':%s" % (noiseSum_w,noiseSum_c))

    ## Compute rations and infer strand
    if mainSum_w == 0 and mainSum_c == 0:
        print("PHAS loci doesn't have phasiRNAs in specied coords - Please debug!!")
        sys.exit()
    else:
        phasiRatio = round((mainSum_w+mainSum_c)/(noiseSum_c+noiseSum_w+mainSum_w+mainSum_c),2)

        if mainSum_w > mainSum_c:
            print("sRNA predicted PHAS strand is 'w'")
            pstrand = 'w'
            strandRatio = round(mainSum_w/(mainSum_w+mainSum_c),3)
        elif mainSum_w < mainSum_c:
            print("sRNA predicted PHAS strand is 'c'")
            pstrand = 'c'
            strandRatio = round(mainSum_c/(mainSum_w+mainSum_c),3)
        else:
            print("Both strands produce same abundace of tags")
            pStrand = 'no'
            strandRatio = round(mainSum_c/(mainSum_w+mainSum_c),3) ## This will be 0.5 since both have eqeual abundance

    print("phasiRatio:%s | strandRatio:%s | predicted strand:%s" % (phasiRatio,strandRatio,pstrand))

    return phasiRatio,strandRatio,pstrand

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

#### Main ##################
############################
def main():

    if aStep == 1:
        ### 1. Parse featureFile
        if mode == 1: ## PASA GTF for overlapping transcripts
            gtfList = gtfParser(featureFile)
        elif mode == 2 or mode == 3 or mode == 5 or mode == 6: ## Trinity or Rocket GTF for overlapping transcipts
            gtfList = gtfParser2(featureFile)
        elif mode == 4: ## Trinity and Rocket both for overlapping transcipts
            gtfList1    = gtfParser2(featureFile)  ## Both lists can be identified with aflag feature
            gtfList2    = gtfParser2(featureFile2) ## Both lists can be identified with aflag feature
            gtfList     = gtfList1 + gtfList2
        else:
            print("Please input correct mode in user settings - script will exit now")
            sys.exit()

        print("Total features cached:%s" % len(gtfList))
        time.sleep(1)
        
        phasList    = phasParser(phasFile)
        resFile     = overlapChecker(phasList,gtfList,phasedFeature,includeLibs)
    

    ## Part-II
    elif aStep == 2:
        resFile = "overlapResults.txt"
        overlapResList = resReader(resFile)
        pass

if __name__ == "__main__":
    main()
    print("Script finished sucessfully")
    sys.exit()

#### Log
## v01 -> v02 [major]
## Added functionality to use Trinity and Rocket results
## Added functionality to use both of above together in single mode = 3
## Added collation of phas status of overlapping transcripts, this required phasing analysis results from transcripts of feature
#### file, like those from rocket and trinity
## Fixed bugs with flags being added  as one tuple, instead of two
## Added functionality to modify phased transcript names from PHAS analysis of rocket transcrits to match gtf file
## Results with no overlapping are added to final results

## v03 -> vo4
## Added functionality to parse pacbio gtf generated from gmap; and work on generic GTF files that have no transcripts entry i.e. the coordinates are computed in gtfParser2 (mode = 3)

## v04 - > v05
## Added funcationality to [rase gtf from comprehensive transcriptome from pasa, mode = 6

## v04 -> v06 [stable]
## Fixed a bug with transcript coords that was mising the start and end coords with each other i.e. a small coord for 'w' was considered as end coord. Foxed by converting bith to 
## integer and then get max and min 

## v06 -> v07[stable]
## Fixed problem with "chr" and "chromosome", now these both are cleaned off from GTF and PHASfile