#!/usr/local/bin/python3

import os,sys,time,sqlite3,operator
import mysql.connector as sql
from collections import Counter

teSubstract     = 1
teFile          = "Aspa.UGA2.csv" ## Tab separated with Name,chr,start,end,strand - No header

mirSubstract    = 1
mirFile         = "Aspa.miR.v2.csv" ## Tab separated with Name,chr,start,end - No header

phasiSubstract  = 1
phasiFile       = "Aspa.PHAS.v2.csv" ## Tab separated with Name,chr,start,end - No header

dataserver      = "tarkan.ddpsc.org"
genomeDB        = "ASPARAGUS_UGA2_genome"

############################
def getDB(afile,conn,aflag):
    '''
    Prepare SQLlite DB for all files
    '''
    cur = conn.cursor()
    featureList = [] ## This is used later

    if aflag == 1: ## Flag comes from 
        ## Read file
        fh_in = open(teFile,'r')
        teRead = fh_in.readlines()

        ## Make table
        featuretable = "teTable"
        cur.execute('''DROP TABLE IF EXISTS %s''' % (featuretable)) ### Drop Old table - while testing    
        conn.commit()
        
        try:
            cur.execute('''CREATE TABLE %s (name varchar(255),chr integer, start integer, end integer, strand varchar(10))''' % (featuretable))
            conn.commit()
            for i in teRead:
                print(i)
                aname,achr,astart,aend,astrand = i.strip("\n").split("\t")
                print ('Name:',aname,'| Chr:',achr,'| Start:',astart,'| End:',aend, '| Strand:',astrand)
                cur.execute("INSERT INTO %s VALUES ('%s',%d,%d,%d,'%s')" % (featuretable,str(aname),int(achr),int(astart),int(aend),str(astrand)))
                featureList.append((str(aname),int(achr),int(astart),int(aend),str(astrand)))
                #conn.commit()
        
        except sqlite3.Error:
            print('ERROR:',Error)
            sys.exit()

        print("-TE table made sucessfuly")

    if aflag == 2:
        ## Read file
        fh_in = open(mirFile,'r')
        mirRead = fh_in.readlines()

        ## Make table
        featuretable = "mirTable"
        cur.execute('''DROP TABLE IF EXISTS %s''' % (featuretable)) ### Drop Old table - while testing    
        conn.commit()
        
        try:
            cur.execute('''CREATE TABLE %s (name varchar(255),chr integer, start integer, end integer, strand varchar(10))''' % (featuretable))
            conn.commit()
            for i in mirRead:
                aname,achr,astart,aend,astrand = i.strip("\n").split("\t")
                print ('Name:',aname,'| Chr:',achr,'| Start:',astart,'| End:',aend, '| Strand:',astrand)
                cur.execute("INSERT INTO %s VALUES ('%s',%d,%d,%d,'%s')" % (featuretable,str(aname),int(achr),int(astart),int(aend),str(astrand)))
                featureList.append((str(aname),int(achr),int(astart),int(aend),str(astrand)))
                #conn.commit()
        
        except sqlite3.Error as er:
            print ('ERROR:', er.message)
            sys.exit()

        print("-miR table made sucessfuly")

    if aflag == 3:
        ## Read file
        fh_in = open(phasiFile,'r')
        phasiRead = fh_in.readlines()

        ## Make table
        featuretable = "phasiTable"
        cur.execute('''DROP TABLE IF EXISTS %s''' % (featuretable)) ### Drop Old table - while testing    
        conn.commit()
        
        try:
            cur.execute('''CREATE TABLE %s (name varchar(255),chr integer, start integer, end integer)''' % (featuretable))
            conn.commit()
            for i in phasiRead:
                print(i)
                aname,achr,astart,aend = i.strip("\n").split("\t")
                print ('Name:',aname,'| Chr:',achr,'| Start:',astart,'| End:',aend)
                cur.execute("INSERT INTO %s VALUES ('%s',%d,%d,%d)" % (featuretable,str(aname),int(achr),int(astart),int(aend)))
                featureList.append((str(aname),int(achr),int(astart),int(aend)))
                #conn.commit()
        
        except sqlite3.Error:
            print('ERROR:',Error)
            sys.exit()

        print("-phasi table made sucessfuly")

    return featuretable,featureList

    ''' This module merges annotations for protein coding, miRNA, phasiRNA and TE's to give 
    a strand specific list of pure intergenic seqeunces'''

    ## Get gene coords #####
    ########################
    cur= con.cursor()
    cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (genomeDB))###extra % escapes the % in query
    genome_info = cur.fetchall() ## List with information of all the genes in genome ## (1, 'c','AT1G01020', 5928, 8737, protein_coding)
    geneCoords = [] ## List to store gene coords
    for i in genome_info:
        aname   = str(i[2])
        achr    = int(i[0])
        astart  = int(i[3])
        aend    = int(i[4])
        astrand = str(i[1])
        geneCoords.append((aname,achr,astart,aend,astrand))

    ## Prepare SQLite DB - Remove old and make new
    ##############################################
    DB = 'tempdb'
    try:
        os.remove(DB)
    except OSError:
        pass
    conn = sqlite3.connect(DB)

    ## Merge TEs
    if teSubstract == 1:
        print("Merging genes and TE's to geenrate intergenic coords")
        teTable,teList = getDB(teFile,conn,1) ## 1 specifies that atable is of TEs is required
        # teTable = "teTable"
        mergedCoords = merger(geneCoords,teTable,teList,conn,'TE') ## merged coords is list with Name,chr,start,end,strand

    if mirSubstract == 1:
        if teSubstract == 1:
            geneCoords = mergedCoords
        print("Merging miRNAs and TE's to geenrate intergenic coords")
        mirTable,mirList = getDB(mirFile,conn,2) ## 2 specifies that atable is of miRs is required
        mergedCoords = merger(geneCoords,mirTable,mirList,conn,'Mi') ## merged coords is list with Name,chr,start,end,strand

    if phasiSubstract == 1:
        if teSubstract == 1 or mirSubstract == 1:
            geneCoords = mergedCoords
        phasiTable,phasiList = getDB(phasiFile,conn,3) ### 3 specifies that atable is of PHAS is required
        mergedCoords = mergerNoStrand(geneCoords,phasiTable,phasiList,conn,'Ph') ## merged coords is list with Name,chr,start,end,strand

    ## Reformat merged list to be compatible with sPARTA functions

    return mergedCoords

def mergeAnnotations(con):
    ''' This module merges annotations for protein coding, miRNA, phasiRNA and TE's to give 
    a strand specific list of pure intergenic seqeunces'''

    ## Get gene coords #####
    ########################
    cur= con.cursor()
    cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (genomeDB))###extra % escapes the % in query
    genome_info = cur.fetchall() ## List with information of all the genes in genome ## (1, 'c','AT1G01020', 5928, 8737, protein_coding)
    geneCoords = [] ## List to store gene coords
    for i in genome_info:
        aname   = str(i[2])
        achr    = int(i[0])
        astart  = int(i[3])
        aend    = int(i[4])
        astrand = str(i[1])
        geneCoords.append((aname,achr,astart,aend,astrand))

    ## Prepare SQLite DB - Remove old and make new
    ##############################################
    DB = 'tempdb'
    try:
        os.remove(DB)
    except OSError:
        pass
    conn = sqlite3.connect(DB)

    ## Merge TEs
    if teSubstract == 1:
        print("Merging genes and TE's to generate intergenic coords")
        teTable,teList = getDB(teFile,conn,1) ## 1 specifies that atable is of TEs is required
        # teTable = "teTable"
        mergedCoords = merger(geneCoords,teTable,teList,conn,"Te") ## merged coords is list with Name,chr,start,end,strand

    if mirSubstract == 1:
        if teSubstract == 1:
            geneCoords = mergedCoords
        print("Merging miRNAs and TE's to generate intergenic coords")
        mirTable,mirList = getDB(mirFile,conn,2) ## 2 specifies that atable is of miRs is required
        mergedCoords = merger(geneCoords,mirTable,mirList,conn,"Mi") ## merged coords is list with Name,chr,start,end,strand

    if phasiSubstract == 1:
        if teSubstract == 1 or mirSubstract == 1:
            geneCoords = mergedCoords
        print("Merging PHAS to generate intergenic coords")
        phasiTable,phasiList = getDB(phasiFile,conn,3) ### 3 specifies that atable is of PHAS is required
        mergedCoords = mergerNoStrand(geneCoords,phasiTable,phasiList,conn,"Ph") ## merged coords is list with Name,chr,start,end,strand

    return mergedCoords

def merger(coords,featureTable,featureList,conn,nameFlag):
    
    aindex          = 0 ## Index of entry to be checked next
    cur             = conn.cursor() ## Connection to SQLite
    mergedList      = [] ## List of merged coordinates
    overlapFeatures = [] ## Features (TE,miR,phasi) that overlapped; and those other then these will be added as is to mergedList

    while aindex < len(coords):
        queryCoord = coords[aindex]
        print("Checking overlap for:",queryCoord)
        
        ## Check for overlap
        # print("Checking overlap for %s" % (queryCoord[0]))
        ## Table Test
        # cur.execute("SELECT * FROM %s where chr = %s and strand = '%s' limit 10" % (table,queryCoord[1],queryCoord[4]))
        # test = cur.fetchall()
        # print(test)
        # sys.exit()

        cur.execute("SELECT * FROM %s where chr = %s AND strand = '%s' AND ((start between %s and %s) OR (end between %s and %s)) ORDER BY end asc" % (featureTable,queryCoord[1],queryCoord[4],queryCoord[2],queryCoord[3],queryCoord[2],queryCoord[3]))
        overlapList = cur.fetchall()

        if not overlapList:
            # print("No overlaps for %s" % (queryCoord[0]))
            mergedList.append(queryCoord)
            aindex+=1 ## Use next query coord
            pass
        else:
            print("There were overlaps",overlapList)
            overlapList.append(queryCoord)
            leftSort    = sorted(overlapList,key=operator.itemgetter(2)) ## sorted on Start
            rightSort   = sorted(overlapList,key=operator.itemgetter(3)) ## sorted on End

            mergedStart     = leftSort[0][2]
            mergedEnd       = rightSort[-1][3]
            mergedChr       = queryCoord[1]
            mergedStrand    = queryCoord[4]
            mergedName      = "%s_%s" % (nameFlag,aindex)
            aindex+=1
            print("Merged Coords- Start:%s | End:%s | Chr:%s | Strand:%s" % (mergedStart,mergedEnd,mergedChr,mergedStrand))
            # if mergedStrand == 'c':
            #     sys.exit()

            ## Add overlapping feature names to overlapFeatures - Fetures other then these will be added directly to merged list
            for i in  overlapList[:-1]: ## Last one is querycoord itself
                overlapFeatures.append(i[0])
            
            ## Continuous loop to check if there are any more overlaps to this new coords
            overlapFlag = 1
            while overlapFlag == 1:
                ## Now check if new merged coords overlap to next querycoord, if yes merge coords and update index
                nextCoord = coords[aindex] ## Update aindex only if used here
                if (mergedChr == nextCoord[1]) and (mergedStrand == nextCoord[4]) and (mergedEnd >= nextCoord[2]):
                    print("Next query coord overlaps with mergedCoords",nextCoord)
                    mergedEnd = nextCoord[3] ## Overlap extended because the merged end overlaps with nextCoord
                    aindex+=1

                    ## Now check if this merged coord further overlaps with TE and update coords, set overlapFlag =1 for further loop
                    cur.execute("SELECT * FROM %s where chr = %s AND strand = '%s' AND ((start between %s and %s) OR (end between %s and %s)) ORDER BY end asc" % (featureTable,queryCoord[1],queryCoord[4],mergedStart,mergedEnd,mergedStart,mergedEnd))
                    overlapList = cur.fetchall()
                    sys.exit()

                    if overlapList: ## Update coords
                        overlapList.append(mergedName,mergedChr,mergedStart,mergedEnd,mergedStrand) ## Add merged coordinates too, to get the end coordinate
                        rightSort       = sorted(overlapList,key=operator.itemgetter(3)) ## sorted on End
                        mergedEnd       = rightSort[-1][3]
                        overlapFlag = 1
                        ## Record overlapping TEs, Phasi, miRNAs - features other then overlapping will be added later
                        for i in  overlapList[:-1]:
                            overlapFeatures.append(i[0])

                    else:
                        print("This is the end of multiple overlaps")
                        overlapFlag = 0

                else:
                    print("No further overlap")
                    overlapFlag = 0

            ## Append final merged to mergedList
            mergedList.append((mergedName,mergedChr,mergedStart,mergedEnd,mergedStrand))

    ## Add those that did not overlapped
    for i in featureList:
        if i[0] not in overlapFeatures:
            print("%s didn't overlap with anybody" % (i[0]))
            mergedList.append(i) ## These are non-overlapping features, i.e not used for merging above
        else:
            print("%s already merged" % (i[0]))

    print("Merged list of len:%s generated" % (len(mergedList)))


    return mergedList

def mergerNoStrand(coords,featureTable,featureList,conn,nameFlag):
    '''This checks for overlap on both strands and provides back a list of merged coords'''
     
    aindex          = 0 ## Index of entry to be checked next
    cur             = conn.cursor() ## Connection to SQLite
    mergedList      = [] ## List of merged coordinates
    overlapFeatures = [] ## Features (TE,miR,phasi) that overlapped; and those other then these will be added as is to mergedList

    while aindex < len(coords):
        queryCoord = coords[aindex]
        print("Checking overlap for:",queryCoord)
        
        ## Check for overlap
        # print("Checking overlap for %s" % (queryCoord[0]))
        ## Table Test
        # cur.execute("SELECT * FROM %s where chr = %s and strand = '%s' limit 10" % (table,queryCoord[1],queryCoord[4]))
        # test = cur.fetchall()
        # print(test)
        # sys.exit()

        cur.execute("SELECT * FROM %s where chr = %s AND ((start between %s and %s) OR (end between %s and %s)) ORDER BY end asc" % (featureTable,queryCoord[1],queryCoord[2],queryCoord[3],queryCoord[2],queryCoord[3]))
        overlapList = cur.fetchall() ## This will include overlap with no strand info, virtual strands need to be considered, i.e. strand of query for merging and and other added as is

        if not overlapList:
            # print("No overlaps for %s" % (queryCoord[0]))
            mergedList.append(queryCoord)
            aindex+=1 ## Use next query coord
            pass
        else:
            print("There were overlaps",overlapList)
            overlapList.append(queryCoord)
            leftSort    = sorted(overlapList,key=operator.itemgetter(2)) ## sorted on Start
            rightSort   = sorted(overlapList,key=operator.itemgetter(3)) ## sorted on End

            mergedStart     = leftSort[0][2]
            mergedEnd       = rightSort[-1][3]
            mergedChr       = queryCoord[1]
            mergedStrand    = queryCoord[4] ## This is the entry from mergedList so it has strand, use this
            mergedName      = "%s_%s" % (nameFlag,aindex)
            aindex+=1
            print("Merged Coords- Start:%s | End:%s | Chr:%s | Strand:%s" % (mergedStart,mergedEnd,mergedChr,mergedStrand))

            ## Add overlapping feature names to overlapFeatures - Fetures other then these will be added directly to merged list
            for i in  overlapList[:-1]: ## Last one is querycoord itself
                overlapFeatures.append((i[0],mergedStrand))
            # print(overlapFeatures)
            # sys.exit()

            ## Continuous loop to check if there are any more overlaps to this new coords
            overlapFlag = 1
            while overlapFlag == 1:
                ## Now chack if new coords overlap to next querycoord, if yes merge coords and update index
                nextCoord = coords[aindex] ## Update aindex only if used here
                if (mergedChr == nextCoord[1]) and (mergedEnd >= nextCoord[2]):
                    print("Next query coord overlaps with mergedCoords",nextCoord)
                    mergedEnd = nextCoord[3] ## Overlap extended because the merged end overlaps with nextCoord
                    aindex+=1

                    ## Now check if this merged coord further overlaps with TE and update coords, set overlapFlag =1 for further loop
                    cur.execute("SELECT * FROM %s where chr = %s AND ((start between %s and %s) OR (end between %s and %s)) ORDER BY end asc" % (featureTable,queryCoord[1],mergedStart,mergedEnd,mergedStart,mergedEnd))
                    overlapList = cur.fetchall()

                    if overlapList: ## Update coords
                        overlapList.append(mergedName,mergedChr,mergedStart,mergedEnd,mergedStrand) ## Add merged coordinates too, to get the end coordinate
                        rightSort       = sorted(overlapList,key=operator.itemgetter(3)) ## sorted on End
                        mergedEnd       = rightSort[-1][3]
                        overlapFlag = 1
                        ## Record overlapping TEs, Phasi, miRNAs - features other then overlapping will be added later
                        for i in  overlapList[:-1]:
                            overlapFeatures.append((i[0],mergedStrand))

                    else:
                        print("This is the end of multiple overlaps")
                        overlapFlag = 0

                else:
                    print("No further overlap")
                    overlapFlag = 0

            ## Append final merged to mergedList
            mergedList.append((mergedName,mergedChr,mergedStart,mergedEnd,mergedStrand))

    ## Add those that did not overlapped, here main features list do not have starnds while overlap list has strands
    ## Make dict
    overlapDict = {} ## This stores those PHAS that were merged and used up
    for i in overlapFeatures:
        overlapDict[i[0]] = i[1] ## Name as value and strand as dict

    for i in featureList: ## This is all fetaures like PHAS with no strands
        aname,achr,astart,aend = i
        
        if aname not in overlapDict.keys(): ## 
            mergedList.append((aname,achr,astart,aend,'w') ) ## These are non-overlapping features, i.e not used for merging above
            mergedList.append((aname,achr,astart,aend,'c') )  ## These are non-overlapping features, i.e not used for merging above
        else:
            astrand = overlapDict[i[0]]
            print("Merged feature %s strand %s" % (i[0],astrand))
            # sys.exit()
            if astrand == 'w':
                mergedList.append((aname,achr,astart,aend,'c')) ## These are non-overlapping features, i.e not used for merging above
            elif astrand == 'c':
                mergedList.append((aname,achr,astart,aend,'w')) ## These are non-overlapping features, i.e not used for merging above
            else:
                print("No strand found - debug module mergerNoStrand")
                sys.exit()



    print("Merged list of len:%s generated" % (len(mergedList)))


    return mergedList

def writer(mergedCoords):
    '''writes file for reading by sPARTA modules'''

    mergedFile = "mergedCoords.txt"
    fh_out = open(mergedFile,'w')

    print("Number of entries:%s" % (len(mergedCoords)))
    acount = 0 ## Count of entries written
    for i in mergedCoords:  ### Expected out format is (1, 'c','AT1G01020', 5928, 8737, protein_coding)
        print(i)
        mergedName,mergedChr,mergedStart,mergedEnd,mergedStrand = i
        fh_out.write("%s\t%s\t%s\t%s\t%s\tmerged\n" % (mergedChr,mergedStrand,mergedName,mergedStart,mergedEnd))
        acount+=1

    print("Number of entries in mergedList:%s | and written:%s" % (acount,len(mergedCoords)))

    fh_out.close()

    return mergedFile

def extractCoords(con,db,mergedFile):
    '''This is the modified function from original sPARTA internal v1.17
    which fetched intergenic cooords given a set of genic coords'''

    cur= con.cursor()
    genomeFeature = 1 ## Since this function will be used to get intergenic coords, hard coding setting here
    
    ## Read merged coords for genes, TEs, miRs and PHAS
    fh_in = open(mergedFile,'r')
    mergedRead = fh_in.readlines()
    fh_in.close()

    genome_info = [] ## Empty list to store genome coords
    for i in mergedRead:
        mergedChr,mergedStrand,mergedName,mergedStart,mergedEnd,trash = i.strip("\n").split("\t")
        genome_info.append((int(mergedChr),str(mergedStrand),str(mergedName),int(mergedStart),int(mergedEnd),trash)) ## (1, 'c','AT1G01020', 5928, 8737, protein_coding)

    print("Genome Info prepared, here is example:")
    print(genome_info[:5])
    # sys.exit()

    ## Check if genome info is empty
    if not genome_info:
        print ('^^^Gene Coords query returned with an empty list..Exiting^^^')
        sys.exit()
    
    ## Find length of chromosomes to calculate intergenics
    ######################################################
    cur.execute('SELECT chr_id, length FROM %s.chromosome_master' % (db))
    chromo_len = cur.fetchall()
    chromo_dict = dict(chromo_len)  ### Made a dict so that chromosome numer could be searched to get chromosome length
    # print ('These are the chromosomes: %s and their length' % (chromo_dict))
    
    genome_info_inter = genome_info ## This list will also hold intergenics
    
    ### GET INTERGENIC REGIONS AND APPEND TO THE LIST WITH GENE COORDS ###
    ######################################################################
    alist = []###list maintained to check if first gene on chromosome and strand shows up than intergenic is just the start of gene
    for i in range(0, int(len(genome_info))+1): ## May 20-14 - Modifed from above after extensive trouble shooting - Now the last entry is read and both up and down calculated
        #print (i)
        gene1 = (genome_info[i])
        gene2 = (genome_info[i+1])
        gene_type = 'inter' ###set to intergenic by default
        #print(gene1,gene2)
        
        ##Remove/skip redundant genes with same start and end......What about Overlapping genes????
        if gene1[3] == gene2[3] and gene1[4] == gene2[4]:
            ##gene is same/overlapping consider next gene
            pass
        
        else:
            ##Calculate coordinates of intergenic regions
            ##Is this is first gene on chromosome and strand intergenic region is from position1 - Only chr_id and strand is checked. 
            if tuple(gene1[0:2]) not in alist:
                print ('\n------Caching gene coords for chromosome: %s and strand: %s------\n' % (gene1[0], gene1[1]))
                #print ('Gene1:%s\nGene2:%s' % (gene1,gene2))
                alist.append((gene1[0:2]))
                
                inter_start1 = 1
                inter_end1 = gene1[3]-1###1 nt before start of Gene1 a.k.a the first gene on chromosome in this case
                ##As two genes are read together, the upstream intergenic region gor gene2 must be calculated in same step
                
                ## If both the genes belong to same chr and strand i.e. chromosome has atleast two genes
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]: 
                    inter_start2 = gene1[4]+1##From end of first gene of chromosome
                    inter_end2 = gene2[3]-1###Till start of second gene
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene2[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])
                    
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                    genome_info_inter.append((gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                
                ## If the first two genes are not from same strand i.e. there is only one gene on this chromosme and strand- code added after Aspragus scaffolds had just one gene
                else: ## intergenic from end of chromosme/scaffold
                    inter_start2 = gene1[4]+1##From end of first gene of chromosome
                    inter_end2 = chromo_dict[gene1[0]]###Till end of chromosome
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_down' % gene1[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])

                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                    genome_info_inter.append((gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene typ
                
            else:
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]:###If chr_id and strands are equal than find intergenic. These are gene on same chromosme and strand
                    inter_start = gene1[4]+1###End of Gene 1
                    inter_end = gene2[3]-1 ###1 nt before start of gene 2
                    if gene2[1] == 'w': ##Positive strand
                        inter_name = ('%s_up' % (gene2[2]))
                    else:## reverse strand
                        inter_name = ('%s_up' % (gene1[2]))
                    #print ('\nLoop3 - Not the first gene on chr and strand')
                    #print (gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type)
                    genome_info_inter.append((gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type))
                
                else: ###That means gene1 is at end of one chromosome and gene 2 is begining of chromosome so we have to extract intergenic at end of one chromosome
                    inter_start = gene1[4]+1###End of gene1
                    inter_end = chromo_dict[gene1[0]]###End of chromosome searched using chromosome id of gene1 from chromosome dictionary
                    if gene1[1] == 'w':##Positive strand end of chromosme
                        inter_name = ('%s_down' % (gene1[2]))
                    else: ##Negative strand first intergenic of chromosme
                        inter_name = ('%s_up' % (gene1[2]))
                        
                    #print ('\nLoop4 - Not the first gene on chromosome and Strand AND the last gene on chromosome')
                    #print (gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)) ## Chr_id, strand


    ## Additional check for scaffolded genomes, if there are no genes in a scffold it's whole seqeunce will be fetched as intergenic
    if genomeFeature == 1:
        for i in chromo_dict.keys():
            alen = chromo_dict[i]
            # print("Chr:%s | Length:%s" % (i,alen))
            if tuple((i,'c')) in alist:
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_c_all' % (i))
                genome_info_inter.append((i,'c',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

            if tuple((i,'w')) in alist:
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_w_all' % (i))
                genome_info_inter.append((i,'w',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

    
    
    ###Sort the list after adding intergenic regions on on basis of chr_id and strand that is essential while caching chromosme during slicing sequences
    genome_info_inter_sort = sorted(genome_info_inter, key=operator.itemgetter(0,1))
    #print(genome_info_inter_sort)
    ###Write all cooords for troubleshooting
    all_coords_out = open('Allcoords', 'w')
    for i in genome_info_inter_sort:
        all_coords_out.write('%s,%s,%s,%s,%s,%s\n' % (i[0:]))
    all_coords_out.close()
            
    
    ###Filter list to remove unwanted types like miRNA,tRNA,rRNA,snoRNA,snRNA, short or no intergenic
    gene_coords_file    = './interCoords'####To check wheter coords are printed in chr_id and strand sorted or not
    fh_out              = open(gene_coords_file, 'w')
    gene_coords         = [] ## List that will hold genes to fetch, this removes unecessary RNAs and also fix miRNA double entry i.e as gene and miRNA
    
    for entry in genome_info_inter_sort:

        ## If RNA type is removed here than that region is not included in analysis but if RNA is removed in mySQL query than only gene is removed and region becomes intergenic
        if genomeFeature == 1: ## Inter
            if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA' or entry [5] == 'protein-coding' or entry [5] == 'protein_coding' or entry [5] == 'misc_RNA' or entry [5] == 'merged'):
                pass
            else:
                if entry[4]-entry[3] > 25:###If there is no intergenic region b/w genes or too short than filter
                    gene_coords.append(entry[0:])
                    fh_out.write('%s,%s,%s,%s,%s,%s\n' % (entry[0:]))
        else: ## Protein coding
            if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA' or entry [5] == 'inter'):
                pass
            else:
                if entry[4]-entry[3] > 25:###If there is no intergenic regon b/w genes or too short than filter
                    gene_coords.append(entry[0:])
                    fh_out.write('%s,%s,%s,%s,%s,%s\n' % (entry[0:]))
            
            
    fh_out.close()

    return gene_coords_file ### A list of selected gene coords

def ConnectToDB(server, infile):
    
    ##infile values are '0' when you dont want to uplaod data from local file and '1' when you wish to upload data by local file
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

def main():
    con = ConnectToDB(dataserver,0)
    # teTable,miRTable,phasiTable = getDB(afile)
    mergedCoords    = mergeAnnotations(con)
    mergedFile      = writer(mergedCoords)
    pureIntergenic   = extractCoords(con,genomeDB,mergedFile)

if __name__ == "__main__":
    main()
    sys.exit()

## getIntergenic v.01
## This script was written to merge coordinates for protein/genes, miRNAs, phasiRNAs and TEs
## This merged coords output is used to fetch intergenic regions without overlap to any of above.
## The output of this script can beused to fetch sRNAs from purely integenic regions

## v01 ->v02 [first release]
## Generates mergedCoords.txt with all genic,TEs,PHAS,miRs merged; and final_coords.txt with pure intergenic coords i.e no overlap with warlier features
## Fixed writer issue due to improper structure of non-verlapping feature list being added to merged list