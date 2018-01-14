#!/usr/local/bin/python3

## Written by Atul Kakrana - Kakrana@udel.edu

## Script takes a (transcript) summary file with start,end,chr and strand info and finds overlap with GTF file 
## Returns file with overlapping genes, type of overlap and orientation of overlap
## Type of overlap - 5' Annotated genes overlaps at 5' of transcript 3' overlaps at 3' of transcript 8' completly enclaves transcript 
## 0' transcript overlaps annotated gene

## Orientation - F: Overlapping gene is at same strand as transcript R: both are on different strands

## USER SETTINGS #######################################

gtf         = "Zea_mays.AGPv3.27.gtf"        ## GTF file matching genome version used for mapping and assembling reads
summary     = "summary.agp2.v3.csv"               ## Summery file with chr,start,stop and strand
delim       = "\t"                           ## Delimiter for summary file
head        = 'Y'                            ## Header is summary file: 'Y' else: 'N'
geneName    = 1
transName   = 2
chromo      = 6
startpos    = 3
endpos      = 4
strand      = 5
geneTypePos = 29
matScorePos = 19                            ## Final annotScore, which comes from matchAnnot parser, be careful this is different from matcAnnot score
makeDB      = 1                             ## If DB for GTF file is not present in present directory then make : 1 else: 0

## IMPORTS ############################################

import os,sys,operator,sqlite3
from operator import itemgetter


## FUNCTIONS ##########################################

def overlapCheck(summary,conn,annotable):

    print("Function: overlapCheck")

    ## Test DB
    cur = conn.cursor()

    ## Test
    # cur.execute("PRAGMA table_info(%s)" % (annotable))
    # desc = cur.fetchall()
    # print(desc)
    # cur.execute("SELECT geneName FROM %s LIMIT 10" % (annotable))
    # test = cur.fetchall()
    # print(test)

    outFile = "%s.overlap.txt" % summary.rpartition(".")[0]
    fh_out = open(outFile,'w')
    fh_out.write("Trans\toverlapGenes\toverlapFLags\toverlapConf\toverlapStatus\n")

    fh_in = open(summary,'r')
    if head == "Y":
        fh_in.readline()
    sumRead = fh_in.readlines()

    for i in sumRead:

        geneList = []       ## Store overlap genes
        flagList = []       ## Store diffrent overlap flags - 5',3', 0' [if gene is enclaved within our transcript] and 8' [if gene extends our transcript at both ends]
        confList = []       ## Store overlap configuration
        statusList = []     ## Store overlap status
        resList = []        ## Store a single merged list fo results

        ent = i.split(delim)
        agene = ent[geneName-1]
        trans = ent[transName-1]
        achr = ent[chromo-1]
        astart = ent[startpos-1]
        aend = ent[endpos-1]
        astrand = ent[strand-1]
        print("\n***Entry:",trans,achr,astart,aend,astrand)

        ## Gene end overlaps
        cur.execute("SELECT * FROM %s WHERE chr = '%s' AND end between %s and %s ORDER BY start asc" % (annotable,achr,astart,aend))
        prime5 = cur.fetchall()
        # print("5Prime\n%s" % prime5)
        tempP5 = [] ## Temp store gene names for checking later
        if prime5:
            for i in prime5:
                tempP5.append(i[4])

        ## Gene start is overlaps
        cur.execute("SELECT * FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY start asc" % (annotable,achr,astart,aend))
        prime3 = cur.fetchall()
        # print("3Prime\n%s" % (prime3))
        tempP3 = []
        if prime3:
            for i in prime3:
                tempP3.append(i[4])

        ## Gene that completly enclaves our transcript   <------ ----trans---- -------->
        cur.execute("SELECT * FROM %s WHERE chr = '%s' AND start < %s AND end > %s ORDER BY start asc" % (annotable,achr,astart,aend))
        prime8 = cur.fetchall()
        # print("8Prime%s\n" % prime8)

        if prime5:
            for i in prime5:
                if i[4] in tempP3:
                    # print("Transcript enclaves the annotated gene'")
                    flag = 0
                    if i[4] != agene:
                        geneList.append(i[4])
                        statusList.append(i[7])
                        flagList.append(flag)
                        if i[3] == astrand:
                            confList.append("F")
                        else:
                            confList.append("R")
                    # print(i)
                    # sys.exit()

                else:
                    # print("Gene overlaps only at one end")
                    if astrand == "+":
                        flag = 5
                    elif astrand == "-":
                        flag = 3
                    else:
                        print("Wrong Strand indicator - Please check file for strand column")
                    # print("Appending prime5:%s" % (i[4]))
                    
                    if i[4] != agene:
                        geneList.append(i[4])
                        statusList.append(i[7])
                        flagList.append(flag)
                        if i[3] == astrand:
                            confList.append("F")
                        else:
                            confList.append("R")
                    # print(i)
        if prime3:
            for i in prime3:
                if i[4] not in tempP5:
                    # print("Gene Overlaps only at one end")
                    
                    if astrand == "+":
                        flag = 3
                    elif astrand == "-":
                        flag = 5
                    else:
                        print("Wrong Strand indicator - Please check file for strand column")
                    
                    if i[4] != agene:                    
                        # print("Appending prime3:%s" % (i[4]))
                        geneList.append(i[4])
                        statusList.append(i[7])
                        flagList.append(flag)
                        if i[3] == astrand:
                            confList.append("F")
                        else:
                            confList.append("R")
                        # print(i)
        if prime8:
            for i in prime8:
                # print("Annotated gene enclaves our transcript")
                # print("Prime8",i)
                flag = 8

                if i[4] != agene: 
                    geneList.append(i[4])
                    statusList.append(i[7])
                    flagList.append(flag)
                    if i[3] == astrand:
                        confList.append("F")
                    else:
                        confList.append("R")

        # print("geneList",geneList,"flagList",flagList,"confList",confList)

        resList = list(zip(geneList,flagList,confList,statusList))
        # print("Final Results",resList)

        # print("FinalRes:%s\t%s\t%s\t%s\n" % (trans,','.join( str(x) for x in geneList), ','.join(str(x) for x in flagList), ','.join(str(x) for x in confList) ))
        if geneList:
           fh_out.write("%s\t%s\t%s\t%s\t%s\n" % (trans,','.join( str(x) for x in geneList), ','.join(str(x) for x in flagList), ','.join(str(x) for x in confList),','.join(str(x) for x in statusList) ))
        else:
            ## There are no overlaps
            fh_out.write("%s\tNA\tNA\tNA\tNA\n" % (trans))


    fh_out.close()
    fh_in.close()

    print("Exiting function - overlapCheck\n")

    return outFile

def gtfParser(gtf):

    print("\nFunction: gtfParser")

    ## file I/O
    fh_in = open(gtf,'r')
    gtfRead = fh_in.readlines()

    parsedGTF = [] ## List to hold parsed GTF entries

    for i in gtfRead:
        if i[0].isdigit():
            ent = i.split("\t")
            if ent[2] == "gene":
                # print(ent)
                achr = ent[0]
                gStart = ent[3]
                gEnd = ent[4]
                gStrand = ent[6]
                info = ent[8].strip("\n").split(";")
                # print(info,len(info))
                if len(info) == 5: 
                    ## Protein coding gene with a version number
                    gid = info[0].split()[1].replace('"','')
                    gVer = info[1].split()[1].replace('"','')
                    gSource = info[2].split()[1].replace('"','')
                    gType = info[3].split()[1].replace('"','')
                    # print(achr,gStart,gEnd,gStrand,gid,gVer,gSource,gType)
                    parsedGTF.append((achr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                elif len(info) == 4:
                    ## miRNA with no version number
                    gid = info[0].split()[1].replace('"','')
                    gVer = "1"
                    gSource = info[1].split()[1].replace('"','')
                    gType = info[2].split()[1].replace('"','')
                    parsedGTF.append((achr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))


            else:
                pass

    print("First 10 entries of parsedGTF list: %s" % (parsedGTF[:10]))
    print("Total entries fetched from GTF file:%s" % (str(len(parsedGTF))))

    print ("Exiting function\n")
    return parsedGTF

def summParserGene(summary):

    print("\nFunction: summParser")

    '''Parses the summary to get gene coordinates, checks if this gene is annotated i.e. matScore > 2, if not
    then adds to list'''

    fh_in = open(summary,'r')
    
    ## Remove header and read the file 
    if head =='Y':
        fh_in.readline()
    summRead = fh_in.readlines()

    ### Sort on Gene Name, and transcript to have all isoforms together
    entList = []        ## To save sorted summary entries
    for i in summRead:
        ent = i.strip('\n').split(delim)
        entList.append(ent)
    summReadSorted = sorted(entList,key=operator.itemgetter(1)) ## Sort on Gene and Trans

    # for i in summReadSorted:
    #     print (i[0],i[1])
    # sys.exit()

    parsedSumm      = []
    tempCodTrans    = []    ## Temporarily store trans for coding transcripts of a gene
    tempCodScore    = []    ## Temporarily store score for coding transcripts of a gene
    tempCodStart    = []
    tempCodEnd      = []

    tempNCodTrans   = []    ## Temporarily store trans for non-coding transcripts of a gene
    tempNCodScore   = []    ## Temporarily store score for non-coding transcripts of a gene
    tempNCodStart   = []
    tempNCodEnd     = []

    acount = 0              ## Count number of pairs read, to check if we are reaching end 
    alen = len(summReadSorted[1:])
    for g1,g2 in zip(summReadSorted,summReadSorted[1:]):
        # print("\n")
        print("G1:",g1,"\n")
        print("G2:",g2,"\n")
        acount+=1
        # print(acount,alen)

        g1Name      = g1[geneName-1]
        g1Trans     = g1[transName-1]
        g2Name      = g2[geneName-1]
        astart      = int(g1[startpos-1])
        aend        = int(g1[endpos-1])
        g1Chr       = g1[chromo-1]
        g2Chr       = g2[chromo-1]
        g1status    = g1[geneTypePos-1]
        g2status    = g2[geneTypePos-1]
        
        ## ## Those that did not match any known transcripts, make their score -1
        scoreInfo  = g1[matScorePos-1]
        if scoreInfo == 'NA':               
            ascore = -1
        else:
            ascore = scoreInfo

        if acount < alen: ## Entries processed are less then total, then this is not last entry      
            if g1Name == g2Name and g1Chr==g2Chr:
                # print("Same genes")
                ## Keep adding info
                
                ## Classify into coding and non-coding transcripts before selectiong a gene locus and unqiue transcripts
                if g1status == "protein_coding" or g1status == "coding" or g1status == "Coding" or g1status == "Cod" or g1status == "NA" or g1status == "co":
                    tempCodTrans.append(g1Trans)
                    tempCodStart.append(astart)
                    tempCodEnd.append(aend)
                    tempCodScore.append(int(ascore))
                elif g1status == "non-coding" or g1status == "Non-coding" or g1status =="Non-Cod" or g1status =="non-cod" or g1status == "nc" or g1status == "up":
                    tempNCodTrans.append(g1Trans)
                    tempNCodStart.append(astart)
                    tempNCodEnd.append(aend)
                    tempNCodScore.append(int(ascore))
                else:
                    print("Wrong Type:%s" % (g1status))
                    sys.exit()


            elif g1Name != g2Name or g1Chr==g2Chr: ## Either a new genes is encountered or end of list is reached
                ## This is the last transcript for this gene
                # print("This is the last isoform under %s gene" % (g1Name))

                ## Classify into coding and non-coding transcripts before selectiong a gene locus and unqiue transcripts
                if g1status == "protein_coding" or g1status == "coding" or g1status == "Coding" or g1status == "Cod" or g1status == "NA" or g1status == "co":
                    tempCodTrans.append(g1Trans)
                    tempCodStart.append(astart)
                    tempCodEnd.append(aend)
                    tempCodScore.append(int(ascore))
                elif g1status == "non-coding" or g1status == "Non-coding" or g1status =="Non-Cod" or g1status =="non-cod" or g1status == "nc" or g1status == "up":
                    tempNCodTrans.append(g1Trans)
                    tempNCodStart.append(astart)
                    tempNCodEnd.append(aend)
                    tempNCodScore.append(int(ascore))
                else:
                    print("Wrong Type:%s" % (g1status))
                    sys.exit()

                ## Find overlappig and unique transcripts, for both coding and non-coding isoforms
                if tempCodTrans:
                    ## If a single gene is coding and non-coding - NOTE: This approach will give error in case there are two small codRNAs, with same gene name and no overlap in coordinates, then start or one and end of one will be recorded making it a long transcript
                    gChr    = g1[chromo-1]
                    gStrand = g1[strand-1]
                    gVer    = "1"                  ## There is nothing like gene version in PacBio Transcripts
                    gSource = "PB"
                    gType   = "co"   ## Summary file should have all coding and non-coding flags

                    ## Collapse isoforms foro a gene, to merged or uniq, for both coding and non-coding RNA
                    mergedList,uniqList = isoformOverlap(tempCodTrans,tempCodStart,tempCodEnd,tempCodScore)
                    
                    if mergedList:
                        mergedEnt = mergedList[0]
                        gid     = mergedEnt[0]
                        gStart  = mergedEnt[1]
                        gEnd    = mergedEnt[2]
                        gScore  = mergedEnt[3]
                        # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                        ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                        if gScore <= 1:
                            parsedSumm.append((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))

                    if uniqList:
                        # print(uniqList)
                        for ent in uniqList:
                            tid     = ent[0]
                            gStart  = ent[1]
                            gEnd    = ent[2]
                            tScore  = ent[3]
                            # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                            ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                            if tScore <= 1:
                                parsedSumm.append((gChr,gStart,gEnd,gStrand,tid,gVer,gSource,gType))


                if tempNCodTrans:
                    ## If a single gene is coding and non-coding - NOTE: This approach will give error in case there are two small ncRNAs, with same gene name and no overlap in coordinates, then start or one and end of one will be recorded making it a long transcript - Check for overlap, if none then record both 
                    gChr = g1[chromo-1]
                    gStrand = g1[strand-1]
                    gVer = "1"                  ## There is nothing like gene version in PacBio Transcripts
                    gSource = "PB"
                    gType = "nc"                ## Summary file should have all coding and non-coding flags

                    ## Collapse isoforms foro a gene, to merged or uniq, for both coding and non-coding RNA
                    mergedList,uniqList = isoformOverlap(tempNCodTrans,tempNCodStart,tempNCodEnd,tempNCodScore) ## Merged list = [name,start,end] and uniqList [(name,start,end),(name,start,end)]
                    
                    if mergedList:
                        mergedEnt = mergedList[0]
                        gid     = mergedEnt[0]
                        gStart  = mergedEnt[1]
                        gEnd    = mergedEnt[2]
                        gScore  = mergedEnt[3]
                        # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                        ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                        if gScore <= 1:
                            parsedSumm.append((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))

                    if uniqList:
                        for ent in uniqList:
                            tid     = ent[0]
                            gStart  = ent[1]
                            gEnd    = ent[2]
                            tScore  = ent[3]
                            # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                            ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                            if tScore <= 1:
                                parsedSumm.append((gChr,gStart,gEnd,gStrand,tid,gVer,gSource,gType))

                ## Empty temp entries for next gene
                tempCodScore    = []
                tempCodStart    = []
                tempCodEnd      = []
                tempCodTrans    = []
                tempNCodScore   = []
                tempNCodStart   = []
                tempNCodEnd     = []
                tempNCodTrans   = []

        if acount == alen:
            ## This is the last set, get details from second gene
            # print("This was the final set")

            bstart = int(g2[startpos-1])
            bend = int(g2[endpos-1])
            bscoreInfo  = g2[matScorePos-1]
            
            ## Those that did not match any known transcripts
            if bscoreInfo == 'NA':               
                bscore = -1
            else:
                bscore = bscoreInfo

            if g1Name == g2Name and g1Chr==g2Chr:
                # print("Same genes")
                
                if g1status == "protein_coding" or g1status == "coding" or g1status == "Coding" or g1status == "Cod" or g1status == "NA" or g1status == "co":
                    tempCodStart.append(astart)
                    tempCodEnd.append(aend)
                    tempCodScore.append(int(ascore))
                elif g1status == "non-coding" or g1status == "Non-coding" or g1status =="Non-Cod" or g1status =="non-cod" or g1status == "nc" or g1status == "up":
                    tempNCodStart.append(astart)
                    tempNCodEnd.append(aend)
                    tempNCodScore.append(int(ascore))
                else:
                    print("Wrong Type:%s" % (g1status))
                    sys.exit()


                if g2status == "protein_coding" or g2status == "coding" or g2status == "Coding" or g2status == "Cod" or g1status == "NA" or g1status == "co":
                    tempCodStart.append(bstart)
                    tempCodEnd.append(bend)
                    tempCodScore.append(int(bscore))

                elif g2status == "non-coding" or g2status == "Non-coding" or g2status =="Non-Cod" or g2status =="non-cod" or g1status == "nc" or g1status == "up":
                    tempNCodStart.append(astart)
                    tempNCodEnd.append(aend)
                    tempNCodScore.append(int(ascore))
                else:
                    print("Wrong Type:%s" % (g1status))
                    sys.exit()

                if tempCodScore:
                    ## Get gene coordinates for this one
                    gChr    = g1[chromo-1]
                    gStrand = g1[strand-1]
                    gVer    = "1"                  ## There is nothing like gene version in PacBio Transcripts
                    gSource = "PB"
                    gType   = "co"                 ## Summary file should have all coding and non-coding flags
                    
                    ## Collapse isoforms foro a gene, to merged or uniq, for both coding and non-coding RNA
                    mergedList,uniqList = isoformOverlap(tempCodTrans,tempCodStart,tempCodEnd,tempCodScore)
                    
                    if mergedList:
                        mergedEnt = mergedList[0]
                        gid     = mergedEnt[0]
                        gStart  = mergedEnt[1]
                        gEnd    = mergedEnt[2]
                        gScore  = mergedEnt[3]
                        # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                        ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                        if gScore <= 1:
                            parsedSumm.append((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))

                    if uniqList:
                        for ent in uniqList:
                            tid     = ent[0]
                            gStart  = ent[1]
                            gEnd    = ent[2]
                            tScore  = ent[3]
                            # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                            ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                            if tScore <= 1:
                                parsedSumm.append((gChr,gStart,gEnd,gStrand,tid,gVer,gSource,gType))
                
                if tempNCodScore:
                    ## Get gene coordinates for this one
                    gChr    = g2[chromo-1]
                    gStrand = g2[strand-1]
                    gVer    = "1"                  ## There is nothing like gene version in PacBio Transcripts
                    gSource = "PB"
                    gType   = "nc"                 ## Summary file should have all coding and non-coding flags
                    ## Collapse isoforms foro a gene, to merged or uniq, for both coding and non-coding RNA
                    mergedList,uniqList = isoformOverlap(tempNCodTrans,tempNCodStart,tempNCodEnd,tempNCodScore) ## Merged list = [name,start,end] and uniqList [(name,start,end),(name,start,end)]
                    
                    if mergedList:
                        mergedEnt = mergedList[0]
                        gid     = mergedList[0]
                        gStart  = mergedList[1]
                        gEnd    = mergedList[2]
                        gScore  = mergedList[3]
                        # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                        ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                        if gScore <= 1:
                            parsedSumm.append((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))

                    if uniqList:
                        for ent in uniqList:
                            tid     = ent[0]
                            gStart  = ent[1]
                            gEnd    = ment[2]
                            tScore  = ent[3]
                            # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                            ## If transcript matchAnnot score <= 1 i.e. novel, not matching genome annotations
                            if tScore <= 1:
                                parsedSumm.append((gChr,gStart,gEnd,gStrand,tid,gVer,gSource,gType))

            ## Last entry of different lonely gene
            elif g1Name != g2Name:
                gStart  = bstart
                gEnd    = bend
                gChr    = g2[chromo-1]
                gStrand = g2[strand-1]
                gVer    = "1"                  ## There is nothing like gene version in PacBio Transcripts
                gSource = "PB"

                if g2[geneTypePos-1] == "protein_coding" or g2[geneTypePos-1] == "coding" or g2[geneTypePos-1] == "Coding" or g2[geneTypePos-1] == "Cod" or g2[geneTypePos-1] == "NA" or g2[geneTypePos-1] == "co":
                    gtype = "co"
                elif g2[geneTypePos-1] == "non-coding" or g2[geneTypePos-1] == "Non-coding" or g2[geneTypePos-1] =="Non-Cod" or g2[geneTypePos-1] =="non-cod" or g2[geneTypePos-1] == "nc" or g2[geneTypePos-1] == "up":
                    gtype = "nc"

                gid = g2Name
                # print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                
                if bScore <= 1:
                    # print("Added. Score:",bScore)
                    parsedSumm.append((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))

    print("Total entries retained from PacBio:%s" % (str(len(parsedSumm))))
    print ("Exiting function\n")
    
    return parsedSumm

def tableMaker(parsedGTF,parsedSumm):
    '''
    Make a track specific sqlite DB for probe ID and coords that will be used
    to query probes on 20MB interval. Each probe entry has following info:
    probe_id,FC,pval,chr,start,end
    '''


    print("\nFunction: tableMaker")

    mergedList = parsedSumm+parsedGTF
    print("Total entries in final merged list:%s" % (str(len(mergedList))))

    annoDB = '%s.db' % (gtf.rpartition('.')[0])
    annotable = "geneMaster"
    conn = sqlite3.connect(annoDB)

    cur = conn.cursor()
    cur.execute('''DROP TABLE IF EXISTS %s''' % (annotable)) ### Drop Old table - while testing    
    conn.commit()

    try:
        cur.execute('''CREATE TABLE %s (chr integer, start integer, end integer, strand varchar(10), geneName varchar(255), geneVersion varchar(255), geneSource varchar(255), geneType varchar(255) )''' % (annotable))
        conn.commit()

        ### Fill the table
        acount = 0 ## COunt of number of gene entries
        for ent in mergedList:
            gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType = ent
            # print("Inserting to sql table:",gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType)
            cur.execute("INSERT INTO %s VALUES (%d,%d,%d,'%s','%s',%d,'%s','%s')" % (annotable, int(gChr), int(gStart), int(gEnd), str(gStrand), str(gid), int(gVer), str(gSource), str(gType) ))
            acount +=1
        
    except sqlite3.Error:
        print('ERROR:',Error)
        sys.exit()

    cur.execute("SELECT * FROM %s LIMIT 10" % (annotable))
    test = cur.fetchall()
    print("First 10 entries:\n%s" % (test))

    cur.execute("SELECT COUNT(*) FROM %s" % (annotable))
    totalEnt = cur.fetchall()
    print("\nTotal entries in table:%s | Total entries in file: %s" % (totalEnt[0][0],acount) )  

    conn.commit() ## Imporatnt to save table


    print("Exiting function\n")
    return annoDB,annotable,conn

def isoformOverlap(transList,startList,stopList,scoreList):
    ''' This function checks for overlap between isoforms to give merged coordinates in case of overlapping geneSource
    and unique transcripts'''

    overlapName     = set() ## To store all uniq trans name for overlapping genes
    overlapStart    = set() ## To store all start sites for overlapping genes
    overlapEnd      = set() ## To store all end sites for overlapping genes
    overlapScore    = set() ## To store all uniq score for overlapping genes, basically for overlapping genes, only best score will be used

    mergedList  = []    ## List to hold gene name, merged start and merged end
    uniqList    = []    ## List to hold (Trans, start and end) of non-overlapping transripts - multiple tuples
    
    # print("\nChecking for isoform overlaps")
    for trans,start,end,score in zip(transList,startList,stopList,scoreList):
        # print('Isoform:',trans,start,end,score,"\n")
        overlapFlag = 0
        ## Check for ovelap with others and seggraget in overlapping or unique
        for x,y,z,a in zip (transList,startList,stopList,scoreList):
            # print('Target:',x,y,z,a)
            if x == trans:
                ## Same transcript
                # print("Same transcript\n")
                pass
            else:
                prime5 = max(int(start),int(y))
                prime3 = min(int(end),int(z))
                
                if (prime3-prime5) > 0:
                    # print("overlap:%s\n" % (prime3-prime5))
                    overlapName.add(trans)
                    overlapStart.add(start)
                    overlapEnd.add(end)
                    overlapScore.add(score)
                    overlapFlag = 1
                    continue

                else:
                    # overlapFlag = 0
                    # print("Do not overlap, trying next isoform")
                    pass

        ## If trans does not overlaps add to Uniq list
        if overlapFlag == 0:
            uniqList.append((trans,start,end,score))


    ## Prepare merged List
    if overlapName:
        ## Some overlaps found
        # print("Overlapping:",overlapName,overlapStart,overlapEnd,overlapScore)
        mergedStart     = min(overlapStart)
        mergedEnd       = max(overlapEnd)
        bestScore       = max(overlapScore)
        name            = trans.rpartition('.')[0] ## Gene name
        mergedList.append((name,mergedStart,mergedEnd,bestScore))

    # print("Uniq List:",uniqList)
    # print("Merged List:",mergedList)
    return mergedList, uniqList

def main():

    if makeDB == 1:
        parsedGTF = gtfParser(gtf)
        parsedSumm = summParserGene(summary)
        annoDB,annotable,conn= tableMaker(parsedGTF,parsedSumm)

    elif makeDB == 0:
        print("Existing annotation DB will be used, to make a new DB please turn on makeDB from settings")
        parsedSumm = summParserGene(summary)
        annoDB = '%s.db' % (gtf.rpartition('.')[0])
        annotable = "geneMaster"
        conn = sqlite3.connect(annoDB)

    else:
        print("'makeDB' variable takes boolean values")

    resFile = overlapCheck(summary,conn,annotable)
    print("Overlap Check complete - see '%s' for results" % (resFile))

if __name__ == '__main__':
    main()

## v01 <-13th July -15
## v01 -> v02
## Added functionality to check overlap with PAcBio transcripts too
## Added fuctionality to reduce PacBio transcripts to gene before checking for overlap
## Only those PacBio genes that have matScore 0,1, or NA are considered non-overlapping to genome annotations and used to check for overlap to 
#### avoid redundancy

## v02 -> v03
## Corrected bug where a PacBio transcript shows overlap to it's own gene - Assumption made that a gene will contain transcripts only from same strand
## Assumed fix for cases when selecting gene coords, had isoforms with no overlap in co-ordinates, it was an error in summary file

### FIX
## Prevent self overlap if PB1.1 overlaps with PB1.1
## NOTE: This approach will give error in case there are two small ncRNAs, with same gene name and no overlap in coordinates, then start or one and end of one will be recorded making it a long transcript