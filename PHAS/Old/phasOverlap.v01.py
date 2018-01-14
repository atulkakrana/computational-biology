#!/usr/local/bin/python3

## This script takes a PHAS file and a GTF file of assembly/genome featires and find overlapping transcripts to PHAS,
## plus computes degree of overlap with reference to PHAS - You can use this identify precursors from assembly 

## Written by Atul Kakrana - kakrana@udel.edu

import os,sys,time,sqlite3,operator
import mysql.connector as sql
from collections import Counter

### User inputs #####
phasFile    = "Final_24PHAS_Loci_ALL_v8.bed"    ## Phased file which needs to be checke for overlap
delim       = "\t"                           ## Delimiter for summary file
head        = 'Y'                            ## Header is summary file: 'Y' else: 'N'
phasName    = 4
chrpos      = 1
startpos    = 2
endpos      = 3
# strandpos   = 5

featureFile = "kakrana_aspa_v23.pasa_assemblies.gtf"    ## Feature file in GTF format, with whome overlap needs to be checked

### Developer inputs ###


#### Functions ########
def tableMaker(alist,conn):
    '''makes SQLlite table for query'''

    print("-Preparing table with features")

    cur = conn.cursor()

    ## Make table
    featuretable = "tempTable"
    cur.execute('''DROP TABLE IF EXISTS %s''' % (featuretable)) ### Drop Old table - while testing    
    conn.commit()
    
    try:
        cur.execute('''CREATE TABLE %s (gene varchar(255),trans varchar(255),chr integer, start integer, end integer, strand varchar(10), type varchar(255))''' % (featuretable))
        conn.commit()
        for i in alist:
            # print(i)
            gid,tid,gchr,gstart,gend,gstrand,gtype = i
            # print ('Gene:',gid,'Trans:',tid,'| Chr:',gchr,'| Start:',gstart,'| End:',gend, '| Strand:',gstrand,'| Type:',gtype)
            cur.execute("INSERT INTO %s VALUES ('%s','%s',%d,%d,%d,'%s','%s')" % (featuretable,str(gid),str(tid),int(gchr),int(gstart),int(gend),str(gstrand),str(gtype)))
            # featureList.append((str(aname),int(achr),int(astart),int(aend),str(astrand)))
            #conn.commit()
    
    except sqlite3.Error:
        print('ERROR:',Error)
        sys.exit()

    print("\n-Feature table made sucessfully")

    return featuretable

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
        achr    = int(ent[chrpos-1].replace('chr',''))
        astart  = int(ent[startpos-1])
        aend    = int(ent[endpos-1])
        # astrand = ent[strandpos-1]
        astrand = "na"
        print("Phas Name %s | chr:%s | start:%s | end:%s | astrand:%s" % (aname,achr,astart,aend,astrand))
        phasList.append(((aname,achr,astart,aend,astrand)))

    print("A list of query features prepared with entries:%s" % len(phasList))

    print ("Exiting function - phasParser\n")
    return phasList

def gtfParser(featureFile):
    '''Parses gtf file into featurename,chr,start,end,strand,feature'''

    print("\nFunction: gtfParser")
    
    with open(featureFile) as fh_in:
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
            info    = ent[8].strip("\n").split(";")
            # print(info,len(info))
            if len(info) == 3: ## With last one empty 
                ## Protein coding gene with a version number
                gid     = info[0].split()[1].replace('"','') ## Gene ID
                tid     = info[1].split()[1].replace('"','') ## Transcript ID
                # print(gid,tid,gchr,gstart,gend,gstrand,gtype)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype))
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

def overlapChecker(phasList,gtfList):
    '''Checks for overlap between genomic PHAS and transcipts frpm GTF file'''
    
    print("\nFunction: overlapChecker")
    
    ## Prpeare outfut file
    outFile = "overlapResults.txt"
    fh_out = open(outFile,'w')
    fh_out.write("PHAS\tphasChr\tphasStart\tphasEnd\tphasStrand\tOverlapTrans\toverlapNucl\toverlapPerc\tstrand\ttransLen\n")
    
    ## Prepare a DB to store feature Table
    DB = 'tempdb'
    try:
        os.remove(DB)
    except OSError:
        pass
    ## Prepare features table 
    conn = sqlite3.connect(DB)
    featureTable = tableMaker(gtfList,conn)

    # featureTable = "tempTable" ## FOr testing when table has been made once
    
    ## Test Query
    # cur = conn.cursor()
    # cur.execute("SELECT * FROM %s where chr = 1 AND strand = 'w' limit 10" % (featureTable))
    # test = cur.fetchall()
    # print("Test results:",test)

    ## Find transcript overlapping PHAS 
    for ent in phasList:
        print("\nEntry:",ent)
        aname,achr,astart,aend,astrand = ent
        transList = overlapTrans(ent,conn,featureTable)

        ## Check overlap with exons of every overlapping transcript
        for trans in transList:
            atrans,tstrand,tlen = trans
            print("\n-Computing overlap for transcript:%s" % (atrans))
            exonsOverlap = overlapExons(ent,conn,featureTable,atrans)
            aperc = round((exonsOverlap/(aend-astart)),2)
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ('\t'.join(str(x) for x in ent),atrans,str(exonsOverlap),str(aperc),astrand,str(tlen)))

    sys.exit()
    fh_out.close()

    return outFile

def overlapTrans(ent,conn,featureTable):
    '''This function returns a list of transcipts that overlaps with a PHAS'''
    transList = [] ## Store the transcipts that overlap

    aname,achr,astart,aend,astrand = ent
    cur = conn.cursor()
    
    ## trans flanking phas or enclaved in phas
    cur.execute("SELECT * FROM %s where chr = %s AND ((end between %s and %s) or (start between %s and %s)) AND type = 'transcript'" % (featureTable,achr,astart,aend,astart,aend))
    flankTrans = cur.fetchall()

    ## PHAS enclaved in trans
    cur.execute("SELECT * FROM %s where chr = %s AND (%s between start and end) AND (%s between start and end) AND type = 'transcript'" % (featureTable,achr,astart,aend))
    bigTrans = cur.fetchall()

    # print("flanking Trans:",flankTrans)
    # print("Long Trans:",bigTrans)

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

def overlapExons(ent,conn,featureTable,atrans):
    '''This function will compute overlaps wit exons of overlapping transcripts'''

    aname,achr,astart,aend,astrand = ent
    
    cur = conn.cursor()
    cur.execute("SELECT * FROM %s where trans = '%s' AND type = 'exon'" % (featureTable,atrans))
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
            print("-Unexpected overlap encountered - Non-overlap?")
            pass
            # sys.exit()

        exonsOverlap+=xoverlap

    print("-Overlap found:%s" % exonsOverlap)
    return exonsOverlap

def main():
    gtfList = gtfParser(featureFile)
    phasList  = phasParser(phasFile)
    resFile = overlapChecker(phasList,gtfList)

if __name__ == "__main__":
    main()
    print("Script finished sucessfully")
    sys.exit()

#### Log
## v01




