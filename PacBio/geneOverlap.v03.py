#!/usr/local/bin/python3

## Written by Atul Kakrana - Kakrana@udel.edu


## Script takes a (transcript) summary file with start,end,chr and strand info and finds overlap with GTF file 
## Returns file with overlapping genes, type of overlap and orientation of overlap
## Type of overlap - 5' Annotated genes overlaps at 5' of transcript 3' overlaps at 3' of transcript 8' completly enclaves transcript 
## 0' transcript overlaps annotated gene

## Orientation - F: Overlapping gene is at same strand as transcript R: both are on different strands


## USER SETTINGS #######################################

gtf = "Zea_mays.AGPv3.27.gtf"   
summary = "Summary.txt"                 ## Summery file with chr,start,stop and strand
delim = "\t"                            ## Delimiter for summary file
head = 'Y'                              ## Header is summary file: 'Y' else: 'N'
geneName = 1
name = 2
chromo = 6
startpos = 3
endpos = 4
strand = 5
geneTypePos = 19
matScorePos = 17
makeDB = 1                             ## If DB for GTF file is not present in present directory then make : 1 else: 0

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
    fh_out.write("Trans\toverlapGenes\toverlapFLags\toverlapConf\n")

    fh_in = open(summary,'r')
    if head == "Y":
        fh_in.readline()
    sumRead = fh_in.readlines()

    for i in sumRead:

        geneList = []       ## Store overlap genes
        flagList = []       ## Store diffrent overlap flags - 5',3', 0' [if gene is enclaved within our transcript] and 8' [if gene extends our transcript at both ends]
        confList = []       ## Store overlap configuration
        resList = []        ## Store a single merged list fo results

        ent = i.split(delim)
        agene = ent[geneName-1]
        trans = ent[name-1]
        achr = ent[chromo-1]
        astart = ent[startpos-1]
        aend = ent[endpos-1]
        astrand = ent[strand-1]
        # print("\n***Entry:",trans,achr,astart,aend,astrand)

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
        # print("8Prime\n%s\n" % prime8)

        if prime5:
            for i in prime5:
                if i[4] in tempP3:
                    # print("Transcript enclaves the annotated gene'")
                    flag = 0
                    if i[4] != agene:
                        geneList.append(i[4])
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
                        flagList.append(flag)
                        if i[3] == astrand:
                            confList.append("F")
                        else:
                            confList.append("R")
                        # print(i)
        if prime8:
            for i in prime8:
                # print("Annotated gene enclaves our transcript")
                # print(i)
                flag = 8

                if i[4] != agene: 
                    geneList.append(i[4])
                    flagList.append(flag)
                    if i[3] == astrand:
                        confList.append("F")
                    else:
                        confList.append("R")

        # print("geneList",geneList,"flagList",flagList,"confList",confList)

        resList = list(zip(geneList,flagList,confList))
        # print("Final Results",resList)

        # print("FinalRes:%s\t%s\t%s\t%s\n" % (trans,','.join( str(x) for x in geneList), ','.join(str(x) for x in flagList), ','.join(str(x) for x in confList) ))
        if geneList:
           fh_out.write("%s\t%s\t%s\t%s\n" % (trans,','.join( str(x) for x in geneList), ','.join(str(x) for x in flagList), ','.join(str(x) for x in confList) ))
        else:
            ## There are no overlaps
            fh_out.write("%s\tNA\tNA\tNA\n" % (trans))


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

def summParser(summary):

    print("\nFunction: summParser")

    '''Parses the summary to get gene coordinates, checks if this gene is annotated i.e. matScore > 2, if not
    then adds to list'''

    fh_in = open(summary,'r')
    if head =='Y':
        fh_in.readline()
    summRead = fh_in.readlines()

    entList = []
    for i in summRead:
        ent = i.strip('\n').split(delim)
        entList.append(ent)
    summReadSorted = sorted(entList,key=operator.itemgetter(0))

    # for i in summReadSorted:
    #     print (i[0])

    parsedSumm = []
    tempScore = []
    tempStart = []
    tempEnd = []

    # for i in range(0,int(len(summReadSorted))+1):
    alen = len(summReadSorted[1:])

    acount = 0
    for g1,g2 in zip(summReadSorted,summReadSorted[1:]):
        # print("\n")
        # print(g1)
        # print(g2)
        acount+=1
        # print(acount,alen)

        g1Name = g1[geneName-1]
        g2Name = g2[geneName-1]
        astart = int(g1[startpos-1])
        aend = int(g1[endpos-1])
        
        scoreInfo  = g1[matScorePos-1]
        if scoreInfo == 'NA':               ## Those that did not match any known transcripts
            ascore = -1
        else:
            ascore = scoreInfo

        if acount < alen:
                ## This is  not last entry        
            if g1Name == g2Name:
                # print("Same genes")
                ## Keep adding info
                
                tempStart.append(astart)
                tempEnd.append(aend)
                tempScore.append(int(ascore))

            elif g1Name != g2Name: ## Either a new genes is encountered or end of list is reached
                # print("Different Genes")
                ## This is the last transcript for this gene
                tempStart.append(astart)
                tempEnd.append(aend)
                tempScore.append(int(ascore))

                ## Get gene coordinates for this one
                gStart = min(tempStart)
                gEnd = max(tempEnd)
                gChr = g1[chromo-1]
                gStrand = g1[strand-1]
                gVer = "1"                  ## There is nothing like gene version in PacBio Transcripts
                gSource = "PB"
                gType = g1[geneTypePos-1]   ## Summary file should have all coding and non-coding flags
                gid = g1Name
                
                print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                print(tempScore,max(tempScore))

                if max(tempScore) <= 1:
                    print("Added")
                    ## Gene doesn't matches with known annotation
                    ## should be added to the list to check for overlaps
                    parsedSumm.append((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))

                    ## Empty temp entries for next gene
                    tempScore = []
                    tempStart = []
                    tempEnd = []

                else:
                    ## This entry matches with known annotation with some confidencem skip
                    tempScore = []
                    tempStart = []
                    tempEnd = []


        if acount == alen:
            ## This is the last set, get details from second gene
            # print("This was the final set")

            bstart = int(g2[startpos-1])
            bend = int(g2[endpos-1])
            bscoreInfo  = g2[matScorePos-1]
            
            if bscoreInfo == 'NA':               ## Those that did not match any known transcripts
                bscore = -1
            else:
                bscore = bscoreInfo

            if g1Name == g2Name:
                # print("Same genes")
                ## Keep adding info
                
                tempStart.append(astart)
                tempEnd.append(aend)
                tempScore.append(int(ascore))

                tempStart.append(bstart)
                tempEnd.append(bend)
                tempScore.append(int(bscore))

                ## Get gene coordinates for this one
                gStart = min(tempStart)
                gEnd = max(tempEnd)
                gChr = g1[chromo-1]
                gStrand = g1[strand-1]
                gVer = "1"                  ## There is nothing like gene version in PacBio Transcripts
                gSource = "PB"
                gType = g1[geneTypePos-1]   ## Summary file should have all coding and non-coding flags
                gid = g1Name
                print((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))
                print(tempScore,max(tempScore))

                if max(tempScore) <= 1:
                    print("Added")
                    parsedSumm.append((gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType))

            elif g1Name != g2Name:
                ## Last entry of different lonely gene
                gStart = bstart
                gEnd = bend
                gChr = g2[chromo-1]
                gStrand = g2[strand-1]
                gVer = "1"                  ## There is nothing like gene version in PacBio Transcripts
                gSource = "PB"
                gType = g2[geneTypePos-1]   ## Summary file should have all coding and non-coding flags
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

def main():

    if makeDB == 1:
        parsedGTF = gtfParser(gtf)
        parsedSumm = summParser(summary)
        annoDB,annotable,conn= tableMaker(parsedGTF,parsedSumm)

    elif makeDB == 0:
        print("Existing annotation DB will be used, to make a new DB please turn on makeDB from settings")
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
## Corrected bug where a PacBio transcript shows overlap to it's own gene - Assumption made that a gene will contain transcripts obly from same strand.




