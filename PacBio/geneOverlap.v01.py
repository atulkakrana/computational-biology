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
name = 2
chromo = 6
start = 3
end = 4
strand = 5
geneType = 20                            ## Column for coding or non-coding, if info not available then mention 99
makeDB = 0                             ## If DB for GTF file is not present in present directory then make : 1 else: 0

## IMPORTS ############################################

import os,sys
import sqlite3

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
        trans = ent[name-1]
        achr = ent[chromo-1]
        astart = ent[start-1]
        aend = ent[end-1]
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
                    # print("Appending prime5:%s" % (i[4]))
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

    print("Function: gtfParser")

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

    print ("Exiting function\n")
    return parsedGTF

def summParser(summary):

    '''Create a a list of summary file'''
    
    print("\nFunction: summParser")
    
    fh_in = open(summary,'r')
    if head == 'Y':
        fh_in.readline()
    summRead = fh_in.readlines()

    geneSet = set()
    # for i in summRead:
    #     ent = i.split("\t")
    #     agene = ent[gene-1]
    #     geneSet.add(agene)

    parsedSumm = []
    acount = 0 ## To count the entries
    for i in summRead:
        # print(i)
        ent = i.split("\t")  ## It has to be tab seprated file always
        agene = ent[gene-1]

        if agene not in geneSet:
            ## New entry add





    print("Total entries scanned: %s | Length of GTF Dictionary %s" % (acount,len(summDict)))
    print("Exiting function - summaryDict\n")
    return parsedSumm

def tableMaker(parsedGTF,parsedSumm):
    '''
    Make a track specific sqlite DB for probe ID and coords that will be used
    to query probes on 20MB interval. Each probe entry has following info:
    probe_id,FC,pval,chr,start,end
    '''

    mergedInfo = parsedGTF + parsedSumm

    print("Function: tableMaker")
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
        for ent in mergedInfo:
            gChr,gStart,gEnd,gStrand,gid,gVer,gSource,gType = ent
            # print(achr,gStart,gEnd,gStrand,gid,gVer,gSource,gType)
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
        annoDB,annotable,conn= tableMaker(parsedGTF)

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




