#!/usr/local/bin/python3

## Written by ATUL to work on IR-based phased siRNAs

import os,sys,operator,time,datetime,string,subprocess,difflib,math
from collections import Counter
import mysql.connector as sql
from multiprocessing import Process, Queue, Pool
import os.path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

### The phasing results must have transcript name, not chr for this script 

### Three BLAST results required from transRNA.sh

runMode     = 1                                             ## 0: Full run - find candidates, validate candidates 1: finalPairs file exists from earlier run now lib fetching, mapping, and plot files need to be generated 3: All earlier file exists, only consensusCands are new
fetchMap    = 0                                             ## Fetch libs and map to fasta file, madatory for runMode 0 
phase       = 24
fastaFile   = '24Phas.fa'
blast_rc    = "Final_BLAST_RC.txt"
blast_nor   = "Final_BLAST_NOR.txt"
blast_comp  = "Final_BLAST_COMP.txt"
phasedRes   = "Final_PHAS_Loci_1e-07_ALL-manualcuration.txt"
clustfile   = "All.txt.score_p1e-07_sRNA_24_out.cluster"    ## Cluster file from phasing analysis, coul dbe for one lib or concatanated file for all libs
phasMode    = "O"                                          ## O: Observed phasiRNAs positions or G: guessed phasiRNA position begining from first phasiRNAs
db          = 'DAYLILY_priv_sRNA'               ## sRNA DB
fetchLibIDs = 'N'                               ## (Y): Get IDs for all libs in DB (N): If you want to run on specific libs than 'N' and specifiy libs below
userLibs    = [(5233,),]                        ## If fetchLib = 'N' then specifiy library IDs in the format [(lib1,),(lib2,),(lib3,),]

consensusRes = 'Y'
consensusCands = "IRcands.v3.list"
matchThres   = 0.25
libType      = 0                                                 ## 0: Lib_ids (4518) | 1: lib_code ('leaf_1')
excludeLibs  = [4518,5006,5003,5004,5005,4519,4521,4523,4525,4527,4520,4522,4524,4526,4528,5230,5232,5234,5235]   ##  Used if fetchMax == 'Y' | Libs you wish to exclude, your tag position summary table should have libs ids and not lib codes

######## Developer settings ####################
minNoiseLen = 18        ## Min noise length included in chart-input, consensus-input and consensusus file, ## Recommended: 18
maxNoiseLen = 32        ## Max noise length included in chart-input, consensus-input and consensusus file, ## Recommended: 24
minPlotLen  = 18        ## Min length considered for plotting
maxPlotLen  = 32        ## Max length used for plotting
gap         = 12        ## Default:12 (from emboss website)
match       = 3         ## Default:3
mismatch    = -4        ## Default:-4
threshold   = 100       ## Default:50
maxRepLen   = 5000      ## Default: 2000
server      = "tarkan.ddpsc.org"             ## Server to use to fetch library information and smallRNA libraries
numProc     = 58

######## FUNCTIONS #############################
def inferIRs(list_rc,header,phasSet):


    fh_out = open("uniqPairs.txt", 'w')
    fh_out.write("%s\ttype\n" % (header))

    rc_sort = sorted(list_rc,key=lambda x: float(x[14]),reverse=True) ## Sorted on bitscore

    for i in rc_sort[0:5]:
        print("Sorted example",i)
    
    ## Retain best pair
    list5           = [] ## Store assigned 5' pairs
    list3           = [] ## Store assgined 3' pairs
    candidateList   = [] ## Store both transcripts of pair
    uniqList        = [] ## Store uniq results

    set5        = set() ## Unique set of transcripts
    set3        = set() ## Unique set of those matched

    acount = 0  ## Total entries count
    for i in rc_sort:
        trans5      = i[0]
        trans3      = i[1]
        pid         = float(i[2])
        length      = int(i[3])
        bitscore    = i[14]
        hang5       = i[15]
        hang3       = i[14]
        match       = i[16]
        minLen      = min(int(len(trans5)),int(len(trans3))) * 0.70 ## 75% or more of length of smaller transcript

        # print("Interaction:%s-%s | PID:%s | Bitscore:%s | AlignLen:%s" % (trans5,trans3,pid,bitscore,length))

        set5.add(trans5)
        set3.add(trans3)
        acount+=1

        if (trans5 not in candidateList) and (trans3 not in candidateList) and length >= 180 and pid >= 85.00 and length >= minLen:
            uniqList.append((trans5,trans3,pid,bitscore,hang5,hang3,match,"IR",i))
            fh_out.write("%s\tIR\n" % ("\t".join(z for z in i)) )

            candidateList.append(trans5)
            candidateList.append(trans3)
        else:
            # print("either hang5 or hang 3 have a partner assigned - Better pair exists")
            pass

    print("\nTotal BLAST pairs:%s | Uniq inferred IR pairs:%s" % (acount,len(uniqList)))
    fh_out.close()

    ## Results
    allSet = set5.union(set3)           ## Total unique phased transcripts in BLAST Results
    candidateSet = set(candidateList)   ## All transcripts that have a pair
    print("Total phased transcripts:%s | Uniq transcripts in BLAST_RC:%s | Candidate paired Transcipts:%s | Uniq paired Transcripts:%s" % (len(phasSet),len(allSet),len(candidateList),len(candidateSet)))

    return uniqList,candidateSet,allSet

def isoformChecker(phasSet,unassignList,pairedSet,isoDict,list_rc,header,phasList):

    '''Checks if unassigned transcripts has paired isoform'''

    noIsoPair   = []    ## Final results storing those transcripts that neither have a confident pair phase transcripts not any of their isoform
                        ## Non-IR transcripts - These are most probably
    isoPair     = []    ## These are isoforms to paired transcipts - Write their best paired isoform from results

    wcount = 0 ## Count unassigned
    xcount = 0 ## Unassigned with some isoform transcript
    ycount = 0 ## Unassigned transcripts that have a paired isoform
    zcount = 0 ## Unassigned transcripts that have no isoform at all
    
    ### Find isoforms of paired transcripts, and keep those unassigned that hare not isoforms to paired for downstream checking
    for i in unassignList:
        isoforms = isoDict[i]
        wcount +=1
        # print("Unassigned Entry:",i,"Isoforms:",isoforms)
        
        if isoforms: ## There is at least one isoform, test if it is paired
            xcount+=1   ## Counts how many unassigned had isoforms reported
            npair = 0   ## Count how many of the isoforms are in paired set
            for ent in isoforms:
                # print("Unassigned:",ent)
                iso = ent[0]
                pid = ent[1]
                mapPerc = ent[2]

                ## Check if it's isoform has been assigned a pair - Then we already captured this loci as IR
                if iso in pairedSet:
                    print("+Isoform for this exists in paired-phas set")
                    # print(i,iso,pid,mapPerc)
                    npair+=1
                else:
                    # print("-This transcript might not be part of IR - lacks a confident pair, and it's putative isofrom too doen't have a pair")
                    # print(i,iso,pid,mapPerc)
                    pass

            ## Check how many of the isofrms of this unassigned were in paired set, if none then this transcript is not isoform of a paired-set
            if npair == 0:
                noIsoPair.append(i)
            else:
                ## This unassigned transcript is an isoform to paired set
                isoPair.append(i)
                ycount+=1
                pass

        else:
            # print("This unassigned transcript does not have any isoform at all")
            noIsoPair.append(i)
            zcount+=1
            pass

    print("Unassigned:%s | w/o isoform:%s | w. isoforms:%s | w. paired isoform:%s" % (wcount,zcount,xcount,ycount))


    print("\nTotal PHAS:%s | Paired PHAS:%s | Unassigned:%s | Non-IR:%s" % (len(phasSet),len(pairedSet),len(unassignList),len(noIsoPair)))

    
    #### Write Results ################
    fh_out      = open("candidateNonIRs.txt", 'w') ## These are transcripts that do not have any pair neither their isoforms have a pair ones
    fh_out2     = open("nonIRs.list", 'w') ## These are unassigned ones that simply do not have BLAST_RC results i.e. no possible IR pair
    
    summaryFile = "Summary_%s.txt" % (datetime.datetime.now().strftime("%m_%d_%H_%M"))
    fh_out3     = open(summaryFile,"w")

    fh_out.write("%s\ttype\n" % (header))
    
    ## Write results for non-IRs
    acount  = 0 ## Count if these non-IRs have some blast-result - These still could be canidate IRs, if we would have had a genome
    bcount  = 0 ## Count if these even have no BLAST-RC results - These could be canddates in which either 5' or 3' arm is not detected
    for trans in noIsoPair:
        blast_res = [] ### Store BLAST RC rsults for noIR
        
        for ent in list_rc:
            query = ent[0]
            sub   = ent[1]
            if query == trans or sub == trans:
                blast_res.append(ent)
                
            else:
                # print("No BLAST results for these final unassigned (Non_IRs):%s" % (trans) ) 
                pass

        if blast_res:
            # print(blast_res)
            for res in blast_res:
                fh_out.write("%s\tiso\n" % ("\t".join(i for i in res)))
            acount += 1
            
        else:
            # print("No BLAST reasults for this putative non-IR transcript:%s" % (trans))
            fh_out2.write("%s\n" % (trans))
            bcount += 1

    print ("Non-IRs with BLAST RC result:%s | Non-IRs with no BLAST results:%s" % (acount,bcount))


    ## Final Summary
    print("\n\n######## SUMMARY ############")
    print("Total PHAS:%s | Total Uniq PHAS:%s | Paired PHAS:%s | Unassigned PHAS:%s" % (len(phasList),len(phasSet),len(pairedSet),len(unassignList)))
    print("--Unassigned - Isoform to paired:%s | Unassigned - Not isoform to paired:%s" % (len(unassignList)-len(noIsoPair),len(noIsoPair)))
    print("----Not Isoform to paired - With BLAST-RC results:%s | Not Isoform to paired - Without BLAST-RC results:%s" % (acount,bcount))
    print("#################################")

    fh_out3.write("######## SUMMARY ############\n")
    fh_out3.write("Total PHAS:%s | Total Uniq PHAS:%s | Paired PHAS:%s | Unassigned PHAS:%s\n" % (len(phasList),len(phasSet),len(pairedSet),len(unassignList)))
    fh_out3.write("--Unassigned - Isoform to paired:%s | Unassigned - Not isoform to paired:%s\n" % (len(unassignList)-len(noIsoPair),len(noIsoPair)))
    fh_out3.write("----Not Isoform to paired - With BLAST-RC results:%s | Not Isoform to paired - Without BLAST-RC results:%s\n" % (acount,bcount))
    

    fh_out.close()
    fh_out2.close()
    fh_out3.close()

    return isoPair,noIsoPair,summaryFile

def isoformWriter(list_nor,isoPair,header,pairedSet):

    print("\nFUNCTION: isoformWriter")

    list_nor = sorted(list_nor,key=lambda x: float(x[14]),reverse=True) ## Sorted on bitscore
    # print("\nSnippet of sorted list:",list_nor[1:5])

    fh_out = open("isoformToPaired.txt", 'w')
    fh_out.write("%s\ttype\tpaired\tisoform\n" % (header))
    
    ### Extract BLAST NORMAL results for transcripts that are isoform to paired
    ############
    acount =    0 ## Count numer of results for isoform written, these should match number of identified isoform to piared
    for trans in isoPair:
        # print("check-0",trans)
        blast_res = [] ### Store BLAST RC rsults for noIR
        
        for ent in list_nor:
            # print("Check-1:",ent)
            query = ent[0]
            sub   = ent[1]
            if query.strip() != sub.strip(): ## It's not a self match
                if query == trans or sub == trans:
                    # print("check-2")
                    blast_res.append((ent))
            else:
                # print("Self match")
                pass

        #### Find the best paired isoform and write results
        ##########
        if blast_res:
            for res in blast_res:
                query   = res[0]
                sub     = res[1]
                if query in pairedSet:
                    fh_out.write("%s\tiso\t%s\t%s\n" % ("\t".join(i for i in res),query,sub))
                    acount += 1
                    break
                elif sub in pairedSet:
                    fh_out.write("%s\tiso\t%s\t%s\n" % ("\t".join(i for i in res),sub,query))
                    acount += 1
                    break ## Write just one best isoform and move to next transcript
                else:
                    pass
        else:
            print("This is assigned as isoform in upstream analysis - but no isoform found -Strange!!!")
            print("Check if BLAST-RESULTS are correct - Investigate")
            sys.exit()

    print("Isoform to paired:%s | Results written:%s" % (str(len(isoPair)),str(acount)))

    fh_out.close()

    return None

def collapseNoIsoforms(isoDict,noIsoPair,header,summaryFile):
    '''This function collapses the noIsoform set to isoform set of unpaired transcript to reduce transcripts to uniq clusters'''

    print("\nFUNCTION: collapseNoIsoforms\n")
    fh_out = open("noIsoToPairedClusts.txt",'w')
    fh_out.write("Trans\tClust_id\tClustFlag\n")
    fh_out2 = open(summaryFile,'a')
    
    assignedSet     = set()                 ## Those that have been assigned a cluster
    unpairedClust   = []
    unclustered     = []                ## Entries that were not asigned clusters, it's possible that their isoforms were already assigned a cluter, if these were matching other isoforms
    
    acount = 0                          ## Those with assigned isoforms
    bcount = 0                          ## Thos with no assigned isoforms
    for ent in noIsoPair:               ## For every unpaired transcript
        if ent not in assignedSet:      ## If no cluster has been assigned yet
            clust       = []                  ## Entry specific list of isoforms forming a cluster
            unclust     = []
            isoforms    = isoDict[ent]

            if isoforms:
                acount +=1
                # print("\nEntry:%s" % (ent))
                # print("Isoforms:",isoforms)

                for i in isoforms:
                    iso = i[0]
                    pid = i[1]
                    mapPerc = ent[2]

                    if iso in noIsoPair: ## If iso form is an unpaired transcipt then assign a cluster
                        # print("+Entry:%s | Isoform:%s" % (ent,iso))
                        acount +=1 ## This isoform is removed from list of unpaired, added to assigned and therefore will not be testedt further. so should be counted at least once
                        assignedSet.add(ent)
                        assignedSet.add(iso)
                        clust.append(ent)
                        clust.append(iso)
                    
                    else:
                        # print("-Entry:%s | Isoform:%s" % (ent,iso))
                        # print("-Isoform is not an noIsoPair member")
                        pass
            else:
                # print("%s - No isoform" % (ent))
                bcount+=1
                unclust.append(ent)
                pass

            unpairedClust.append((clust))
            unclustered.append(unclust)

    print("\nTotal noIsoPair Transcripts:%s | Trans with isoforms:%s | Without isoforms:%s" % (str(len(noIsoPair)),acount,bcount))
    print("Clusters identified:%s" % (str(len(unpairedClust)) ))
    fh_out2.write("----Not Isoform to paired - Trans with unpaired isoforms:%s | w/o unpaired isoforms:%s \n" % (acount,bcount))
    fh_out2.write("----Not Isoform to paired - collapsed to clusters:%s\n" % (str(len(unpairedClust)) ))
    fh_out2.write("#################################\n")

    ## Write a text file with transcripts that are not isoform to paired ones along with a cluster ID representing these are isoforms to others in same clusters
    acount = 0 ## Cluster number
    allClusts = unclustered+unpairedClust
    for clust in allClusts:
        acount += 1
        if len(clust) > 1: #### What if it's a list of (trans,isoform) but since length is one it will pass to next loop - Here we miss a few unpaired in noIsoToPairedClusts.txt
            ## Check if it was clustered and add a flag
            for i in clust:
                fh_out.write("%s\t%s\tY\n" % (i,acount))
        else:
            for i in clust:
                fh_out.write("%s\t%s\tN\n" % (i,acount))

    fh_out.close()
    fh_out2.close()

    return None

def validateIR (clustfile,uniqList,phasList,fastaFile,header,phasSet):

    ## Get clusters list
    fh_in       = open(clustfile,'r')
    clusters    = fh_in.read().split('>')
    fh_in.close()

    ## Get phased fasta list
    fastaList,fastaDict = fastaReader(fastaFile)

    ## Prepare output file
    resfile     = 'finalPairs.txt'
    fh_out      = open(resfile,'w')
    fh_out.write("%s\ttype\ttotalPhasi\tvalid5\tvalid3\tfinalStatus\n" % (header))

    clustFile    = "phased.clust"
    fh_out1     = open(clustFile,'w')


    ### Get the phasing information for pairs ###
    validList  = [] ## List that will store results of this module for further connectivity
    pairedList = [] ## List to store final IR pairs 
    for pair in uniqList:
        trans5,trans3,pid,bitscore,hang5,hang3,match,atype,blastent = pair
        print("\nFetching phased loci for both partners:%s-%s"% (trans5,trans3))

        ### Get phased entry
        trans5Phas = []
        trans3Phas = []
        
        if trans5 != trans3: ## Not self foldback i.e. not true IR
            for i in phasList:
                name,pval,trans,start,end,strand,lib = i
                if trans == trans5:
                    trans5Phas.append((i))
                    # aphasistrand = 'w'
                elif trans == trans3:
                    trans3Phas.append((i))
                    # bphasistrand = 'w'
                else:
                    # print("No match")
                    pass
        else: ## Self foldback, True IR
            for i in phasList:
                name,pval,trans,start,end,strand,lib = i
                if trans == trans5: ## Get all phased coords for this direct IT, which might be a single entry or two entries with different phasing positions
                    trans5Phas.append((i)) ## Same transcipt as trans5 and same as trans3 - Phased coords could be same or different
                    trans3Phas.append((i))
                    # aphasistrand = 'w' ## Need to catch both strands
                    # bphasistrand = 'c'

        print("Matching trans5 phas:%s" % (str(len(trans5Phas))))
        print("Matching trans3 phas:%s" % (str(len(trans3Phas))))
        # print("Phased loci info fetched for %s-%s pair" % (trans5,trans3))
        print("Phased info:%s:%s" % (trans5,trans5Phas))
        print("Phased info:%s:%s" % (trans3,trans3Phas))

        ## Test must be performed on every phased loci of a pair ##
        for aphas in trans5Phas:
            for bphas in trans3Phas:
                
                ## Get phased coordinates for both 
                acoords,aclust = getPhase(aphas,clusters)
                bcoords,bclust = getPhase(bphas,clusters)
                if acoords:
                    print("Coords for %s:%s" % (aphas,acoords))
                else:
                    print("Halt - No coordinates reported for %s" % (aphas))
                    sys.exit()

                if bcoords:
                    print("Coords for %s:%s" % (bphas,bcoords))
                else:
                    print("Halt - No coordinates reproted for %s" % (bphas))
                    sys.exit()

                ## Get foldback coordinates
                IRcoords,trans5len,trans3len = getIRcoords(fastaList,trans5,trans3)
                print("+These are einverted results:%s" % (IRcoords))
                print("+This trans5Len:%s | Trans3Len:%s" % (trans5len,trans3len))

                ## Convert trans3 coords to IR coords, filter valid coords that lie within IR regions, and tranform
                ## trans3 phased coords to refer from left hand side as in foldback struture

                ## Convert trans3 phased coords to correspond to IR with trans5+trans3 #############
                ####################################################################################
                # print(trans5len,bcoords)
                if trans5 != trans3: ## Two transcipts, not a direct IR
                    convBcoords = [int(x)+trans5len for x in bcoords]
                else: ## Self foldback i.e. direct IR
                    convBcoords = [int(x) for x in bcoords]

                print("+Trans5 phase coords: %s" % (acoords))
                print("+Trans3 converted coords: %s" % (convBcoords))

                ## Filter valid IR coords ##########################################################
                ####################################################################################
                start1      = int(IRcoords[0][6])
                end1        = int(IRcoords[0][7])
                start2      = int(IRcoords[0][8])
                end2        = int(IRcoords[0][9])
                print("+IR Coords - Start1:%s | End1:%s | Start2:%s | End2:%s" % (start1,end1,start2,end2))

                avalidcoords = [] ## Store coords that lie within 5'arm of IR
                bvalidcoords = [] ## Store coords that lie within 3'arm of IR
                for coord in acoords:
                    if coord >= start1 and coord <= end1:
                        avalidcoords.append(coord)

                for coord in convBcoords:
                    if coord >= start2 and coord <= end2:
                        bvalidcoords.append(coord)
                
                print("+Valid trans5 coords: %s" % (avalidcoords))
                print("+Valid trans3 coords: %s" % (bvalidcoords))

                if not avalidcoords or not bvalidcoords:
                    print("+Phased coords of either arm are not from foldback region - Not an IR-based PHAS")
                    totalPhasi  = len(acoords)+len(bcoords)
                    finalstatus = 'No'
                    validList.append((blastent,totalPhasi,len(avalidcoords),len(bvalidcoords),finalstatus))
                    fh_out.write("%s\tIR\t%s\t%s\t%s\t%s\n" % ("\t".join(z for z in blastent),totalPhasi,len(avalidcoords),len(bvalidcoords),finalstatus))
                    continue
                else:
                    pass

                ### Write the phased cluster to a file, in case of direct IR, phasiRNA and position is used to remove redundant entries
                ### If there are two diffrent phased positions on same transcript the, only the combination with valid 5' and 3' coords return,
                ### other three combinations will either have no valid coords, or empty 5' or 3' coords. So phased cluster file will have just one 
                ### cluster for this direct IR with all phasiRNAs with uniq positions
                phasedClusterWriter(trans5,trans3,aclust,bclust,fh_out1)

                ## Normalize coords ###############################################################
                ## trans5 and trans3 phased coords to correspond to left IR start
                anormcoords = [(int(x)-start1)+1 for x in avalidcoords] ## +1 to include first base of IR arm
                bnormcoords = [abs(int(x)-end2)+1 for x in bvalidcoords[::-1]] ## Coords flipped to correposnd to foldback trans5 from left to right, +1 to include first base of IR arm
                print("+Norm trans5 coords: %s" % (anormcoords))
                print("+Norm trans3 coords: %s" % (bnormcoords))

                ## Trim coords ####################################################################
                ## Choose overlapping coords between trans5 and trans3 ############################
                amin = min(anormcoords)
                bmin = min(bnormcoords)
                amax = max(anormcoords)
                bmax = max(bnormcoords)

                minTrim = int(max(amin,bmin))-3 ## Trick is to trim trans5 or trans3 coords to an overlapping left hand, and include any smaller matching +2nt phase from other
                maxTrim = int(min(amax,bmax))+3 ## Here trick is to trim to max overlapping phase, plus include +2nt phase from other arm

                atrimcoords = [] ## List to store trans5 coords overlapping with trans3
                btrimcoords = [] ## List to store trans3 coords overapping with trans5

                for coord in anormcoords:
                    if coord >= minTrim and coord <= maxTrim:
                        atrimcoords.append(coord)

                for coord in bnormcoords:
                    if coord >= minTrim and coord <= maxTrim:
                        btrimcoords.append(coord)
                print("+Trim trans5 coords: %s" % (atrimcoords))
                print("+Trim trans3 coords: %s" % (btrimcoords))

                if not atrimcoords or not btrimcoords:
                    print("+Phased coords from both arms do not overlap - Possible IR-based PHAS")
                    totalPhasi  = len(acoords)+len(bcoords)
                    finalstatus = 'Fo'
                    validList.append((blastent,totalPhasi,len(avalidcoords),len(bvalidcoords),finalstatus))
                    fh_out.write("%s\tIR\t%s\t%s\t%s\t%s\n" % ("\t".join(z for z in blastent),totalPhasi,len(avalidcoords),len(bvalidcoords),finalstatus))
                    continue
                else:
                    pass

                ## Test for overhang ##############################################################
                ## This considers possibility if afew phases are missing on either arm ############
                finalstatus = checkOverhang(atrimcoords,btrimcoords)
                print("+Final status of this pair: %s" % (finalstatus))

                ## This information will be used for further summarization of results
                if finalstatus == "IR":
                    pairedList.append(trans5)
                    pairedList.append(trans3)

                ## Write results #################################################################
                ## uniq pair entry, total phasiRNA positions, valid phasiRNAs 5' and 3', and status P or NP
                totalPhasi = len(acoords)+len(bcoords)
                validList.append((blastent,totalPhasi,len(avalidcoords),len(bvalidcoords),finalstatus))
                fh_out.write("%s\tIR\t%s\t%s\t%s\t%s\n" % ("\t".join(z for z in blastent),totalPhasi,len(avalidcoords),len(bvalidcoords),finalstatus))

                ## For further summarization of results
                pairedSet = set(pairedList) ## All transcripts that have a pair
                ## Identify unpaired from total phased transcripts
                unassignList = [] ## List to store those transcripts that are missing from pair-test or didn't had unoq pair
                for i in phasSet:
                    if i not in pairedSet:
                        unassignList.append(i)
                print("Total phased transcripts:%s | Unassigned transcripts:%s" % (len(phasSet),len(unassignList)))


    fh_out.close()
    fh_out1.close()
    print("These are 5 entries of validList:%s" % (validList[:1]))

    return validList,resfile,clustFile,pairedSet,unassignList

def prepareChartFiles(validList,clustFile,userLibs,fastaFile,con):
    '''
    This function will prepare chart files for every pair, which will include phasiRNAs
    and non-phased sRNAs
    '''

    cur= con.cursor()
    queryLibs,sumLibs,finalLibs = prepareQuery(excludeLibs,cur)

    print("\nFUNCTION: prepareChartFiles")
    ## Get phased fasta list
    fastaList,fastaDict = fastaReader(fastaFile)

    ## Fetch phased phasiRNAs from cluster file and make a dictionary with transcript name as key
    ## and phasiRNAs as values
    phasedDict = clust2Dict(clustFile)
    
    ## Fetch sRNA dict for specified library
    for lib in userLibs:

        ## Prepare for writing results
        chartInput      = "%s.chart.input.txt" % (lib)
        fh_out          = open(chartInput,'w')
        fh_out.write("pairname\tphasiname\tphasistrand\tphasiabun\tphasiseq\tphasilen\tphasihits\tphasiflag\tphasipos\tphasiarm\n") ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos,phasiarm

        consensusInput  = "%s.consensus.input.txt" % (lib)
        fh_out2          = open(consensusInput,'w')
        fh_out2.write("pairname\tphasiname\tphasistrand\tphasiabun\tphasiseq\tphasilen\tphasihits\tphasiflag\tphasipos\tphasiarm\tconsensusFlag\n") ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos,phasiarm

        ### Parse bowtie map file for this library, and prepare a dictionary with transcript name as key
        ## and phasiRNAs as values
        sRNADict = map2Dict(lib)

        ## Start the process of making a trans5cluster and trans3cluster list
        for pair in validList:
            # print(pair)
            astatus = pair[-1]
            if astatus == 'IR':
                trans5 = pair[0][0] ## ([blastent],totalPhasi valid5 valid3 finalStatus)
                trans3 = pair[0][1]
                print("\n+IR pair:%s-%s | status:%s being processed" % (trans5,trans3,astatus))
                atransclust = mergeRNAs(trans5,sRNADict,phasedDict,fastaDict,cur,finalLibs)
                btransclust = mergeRNAs(trans3,sRNADict,phasedDict,fastaDict,cur,finalLibs)

                ## Get IR coords ####################################################
                #####################################################################
                IRcoords,trans5len,trans3len = getIRcoords(fastaList,trans5,trans3)
                print ("+These are einverted results:%s" % (IRcoords))
                
                ## Convert coordinates to IR transcript ############################
                ####################################################################

                # print(trans5len,bcoords)
                convBclust      = []            ## List with converted coords for trans3
                if trans5       != trans3:      ## Two transcipts, not a direct IR
                    for phasi in btransclust:
                        apos    = phasi[7]      ## This is the position
                        convPos = int(apos)+trans5len
                        convBclust.append((phasi[0],phasi[1],phasi[2],phasi[3],phasi[4],phasi[5],phasi[6],convPos)) ## phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos

                else: ## Self foldback i.e. direct IR, no need to convert coordinates
                    convBclust = list(btransclust) ## Just copy the list to a new name

                # print("+trans5len:%s | trans3len:%s" % (trans5len,trans3len))
                # print("+Trans3 normal coords: %s" % (trans3Clust))
                # print("+Trans3 converted coords: %s" % (convBclust))


                ## Filter valid coords i.e. those fall between foldback ###############
                #######################################################################
                start1      = int(IRcoords[0][6])
                end1        = int(IRcoords[0][7])
                start2      = int(IRcoords[0][8])
                end2        = int(IRcoords[0][9])
                print("+IR Coords - Start1:%s | End1:%s | Start2:%s | End2:%s" % (start1,end1,start2,end2))

                avalidclust = [] ## Store phasiRNAs that lie within 5'arm of IR
                bvalidclust = [] ## Store phasiRNAs that lie within 3'arm of IR
                for phasi in atransclust:
                    apos = phasi[7]
                    if apos >= start1 and apos <= end1:
                        avalidclust.append(phasi)

                for phasi in convBclust:
                    apos = phasi[7]
                    if apos >= start2 and apos <= end2:
                        bvalidclust.append(phasi)
                
                # print("\n+Snippet of valid trans5 clust: %s" % (avalidclust))
                # print("\n+Snippet of valid trans3 clust: %s" % (bvalidclust))
                print("+All trans5 sRNAs:%s | valid trans5 sRNAs:%s" % (len(atransclust),len(avalidclust)))
                print("+All trans3 sRNAs:%s | valid trans3 sRNAs:%s" % (len(btransclust),len(bvalidclust)))

                ## Sanity check - #1
                if not avalidclust or not bvalidclust:
                    print("+Phased coords of either arm are not from foldback region - Not an IR-based PHAS")
                    continue
                else:
                    pass

                ## Normalize coords ###################################################
                ## trans5 and trans3 phased coords to correspond to left of IR start ###
                anormclust = [] ## Store phasiRNAs for trans5, with normalized positions that coorspond from left of IR
                bnormclust = [] ## Store phasiRNAs for trans3, with normalized positions that coorspond from left of IR
                for phasi in avalidclust:
                    apos    = phasi[7]
                    normpos = (int(apos)-start1)+1 ## +1 to include first base of IR arm
                    anormclust.append((phasi[0],phasi[1],phasi[2],phasi[3],phasi[4],phasi[5],phasi[6],normpos)) ## phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos

                for phasi in bvalidclust[::-1]: ## Iterate in reverse direction, so that print statment corresponds to IR left to right
                    apos        = phasi[7]
                    normpos     = abs(int(apos)-end2)+1 ## +1 to include first base of IR arm
                    bnormclust.append((phasi[0],phasi[1],phasi[2],phasi[3],phasi[4],phasi[5],phasi[6],normpos))    ## Coordinates flipping not required like validate 
                                                                                                                    ## phas because these will be automatically arranged while plotting
                # print("\n+Snippet of Norm trans5 clust: %s" % (anormclust))
                # print("\n+Snippet of Norm trans3 clust: %s" % (bnormclust))

                ## Normlize file to align multiple IRs
                ########################################################################
                if consensusRes == 'Y':
                    anormclust_s = sorted(anormclust,key=operator.itemgetter(6,7)) ## sorted on sRNAs/phased RNAs and then on Pos
                    bnormclust_s = sorted(bnormclust,key=operator.itemgetter(6,7)) ## sorted on sRNAs/phased RNAs and then on Pos
                    print("\n+anormclust_s",anormclust_s)
                    print("\n+bnormclust_s",bnormclust_s)

                    aflag,bflag                     =  flipflag(anormclust_s,bnormclust_s)
                    aconsensusclust,bconsensusclust =  consensusNorm(anormclust_s,bnormclust_s,aflag,bflag)
                    # print("\n+Snippet of Consensus trans5 clust: %s" % (aconsensusclust))
                    # print("\n+Snippet of Consensus trans3 clust: %s" % (bconsensusclust))

                ### Write merged normalized sRNAs results for this pair for plotting
                combName = "%s-%s" % (trans5,trans3)
                for phasi in anormclust :
                    # print(phasi)
                    phasiarm = "5"
                    fh_out.write("%s\t%s\t%s\n" % (combName,'\t'.join(str(x) for x in phasi),phasiarm)) ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos,phasiarm

                for phasi in bnormclust:
                    phasiarm = "3"
                    fh_out.write("%s\t%s\t%s\n" % (combName,'\t'.join(str(x) for x in phasi),phasiarm)) ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos.phasiarm

                ### Prepare abundance chart for this pair
                chart(trans5,trans3,anormclust,bnormclust)


                ## Write results that are normalized to align different IRs together
                if consensusRes == 'Y':
                    for phasi in aconsensusclust :
                        # print(phasi)
                        phasiarm = "5"
                        fh_out2.write("%s\t%s\t%s\t%s\n" % (combName,'\t'.join(str(x) for x in phasi),phasiarm,aflag)) ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos,phasiarm

                    for phasi in bconsensusclust:
                        phasiarm = "3"
                        fh_out2.write("%s\t%s\t%s\t%s\n" % (combName,'\t'.join(str(x) for x in phasi),phasiarm,bflag)) ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos,phasiarm


    fh_out.close()
    fh_out2.close()

    return chartInput,consensusInput

def mergeIRs(consensusInput,consensusCands):
    '''
    This function will read consesusInput file (".consensus.input.txt") and merge abundances for all IRs
    in mergeFile, the IRs' will be detected in combined name, if present then their abundance is added at position
    of either A-side or B-side depending upon the aflag and bflag

    consensusInput : File with coords normalized to alingn IRs and flip flags
    IRList : A list of transcripts whose IRs need to be used to generate consensus file
    '''
    print("\nFUNCTION:mergeIRs")

    ## Read inputs
    fh_in = open(consensusInput,'r')
    fh_in.readline() ## remove header
    consensusRead = fh_in.readlines()

    fh_in = open(consensusCands,'r')
    IRList = fh_in.read().split("\n")
    print(IRList)

    ## Prepare to write results ###############################################
    ###########################################################################
    mergedOut       = "consensus.merged.txt" ## All size abundances merged together
    mergedsizeOut   = "consensus.merged.size.txt" ## size-specific abundances merged together

    fh_out = open(mergedOut,'w')   
    fh_out.write("Pos\tAbun\tFlag\n")

    fh_out2 = open(mergedsizeOut,'w')
    fh_out2.write("Pos\tAbun\tFlag\tSize\n")


    ## Extract combined names, and get max position i.e highest coord for trans in IRList
    ##################################################################################
    combnameSet = set() ## List of comb-names to work on
    allPos      = []    ## Store all coords for consensus candidates 
    
    for i in consensusRead:
        ent             = i.strip('\n').split('\t') ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos,phasiarm,a/bflag
        # print("+ent",ent)
        combname        = ent[0]
        trans5,trans3   = ent[0].split("-")

        # acount = 0 ## Pairs collected for merging
        if trans5 in IRList or trans3 in IRList:
            # print("+Found pair %s for input candidate %s,%s" % (combname,trans5,trans3))
            flipflag        = ent[10]
            
            if flipflag != 'N': ## Which means assigned flipFlag is N(No), this pair cant be merged
                phasiflag       = ent[7]
                combnameSet.add(combname)
                # acount+=1
            
                if phasiflag == 'P': ### If the phasiRNAs is phased, no point in add sRNAs at the tail i.e. after phasing end
                    apos        = int(ent[8]) ## Record phased position to define the length of consensus positions
                    allPos.append(apos)
            else:
                # print("flipFlag for this couldn't be determined - So entry will be skipped")
                pass
    
    print("+Pairs collected for merging:%s | Pairs in input list:%s" % (len(combnameSet),len(IRList)))


    ## Prepare empty lists to store positional abundnaces, and size wise positional abundances
    ##########################################################################################

    maxLen = int(max(allPos)) ## Position with last phased loci recorded
    alist = [0]*int(maxLen+1) ## alist start from 0, but positions to be stored will start from 1, so one extra position in list 
    blist = [0]*int(maxLen+1) ## To account for zero offset of python, list is extended by 1
    anestlist = [[0,0,0,0] for x in range(maxLen+1)] ## Nested list with positions, and nest with 4 slots for 21,22,23,and 24nt size abundnaces
    bnestlist = [[0,0,0,0] for x in range(maxLen+1)] ## Nested list with positions, and nest with 4 slots for 21,22,23,and 24nt size abundnaces
    print("Length of positional lists:%s and %s" %  (len(alist),len(blist)))

    ## Merge abundances for each position, and size-wise for each position
    #######################################################################
    for i in consensusRead:
        ent         = i.strip("\n").split('\t') ## combname,phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos,phasiarm,a/bflag
        combname    = ent[0]
        apos        = int(ent[8])
        # print("+checking %s in canididate list" % (combname))
        
        if combname in combnameSet: ## If pair includes a transcript in user input IR.list
            # print("+found %s in canidate list" % (combname))
            if apos <= maxLen: ## Record only those sRNAs or phased abundances whose position is till last phased siRNAs, no need to record sRNA when there is no phasiRNA on precursor 
                aabun   = int(ent[3])
                alen    = int(ent[5])
                xflag   = ent[10]
                # print("+Position:%s | Abun:%s | Flag:%s" % (apos,aabun,xflag))
                if xflag    == 'A':
                    alist[apos] += aabun ## All size combined, the abun will be slightly more then sum of size wise abun, becuase there might be additional phased size other then those 21,22,23,24
                    if alen == 21:
                        anestlist[apos][0]+=aabun
                    elif alen == 22:
                        anestlist[apos][1]+=aabun
                    elif alen == 23:
                        anestlist[apos][2]+=aabun
                    elif alen == 24:
                        anestlist[apos][3]+=aabun
                
                elif xflag  == 'B':
                    blist[apos] += aabun ## All size combined, the abun will be slightly more then sum of size wise abun, becuase there might be additional phased size other then those 21,22,23,24
                    if alen == 21:
                        bnestlist[apos][0]+=aabun
                    elif alen == 22:
                        bnestlist[apos][1]+=aabun
                    elif alen == 23:
                        bnestlist[apos][2]+=aabun
                    elif alen == 24:
                        bnestlist[apos][3]+=aabun
                
                else:
                    print("The flag for side-A/B is not recognized:%s" % (xflag))
                    sys.exit()

    ## Write final file with virtual positions and abundances
    #########################################################
    acount = 0 ## Position count
    for i,j in zip(alist[1:],anestlist[1:]): ## Not include the 0 position of list, because there is no 0 position in IRs
        acount  +=1
        flag    = 'A' ## Since this is Alist
        fh_out.write("%s\t%s\t%s\n" % (acount,i,flag)) ## File with abundance for all sizes merged

        fh_out2.write("%s\t%s\t%s\t21\n" % (acount,j[0],flag)) ## File with abundance for all sizes merged
        fh_out2.write("%s\t%s\t%s\t22\n" % (acount,j[1],flag)) ## File with abundance for all sizes merged
        fh_out2.write("%s\t%s\t%s\t23\n" % (acount,j[2],flag)) ## File with abundance for all sizes merged
        fh_out2.write("%s\t%s\t%s\t24\n" % (acount,j[3],flag)) ## File with abundance for all sizes merged
    
    acount = 0
    for z,x in zip(blist[1:],bnestlist[1:]): ## Not include the 0 position of list, because there is no 0 position in IRs
        acount  +=1
        flag    = 'B' ## Since this is Blist
        fh_out.write("%s\t%s\t%s\n" % (acount,z,flag)) ## File with abundance for all sizes merged
        
        fh_out2.write("%s\t%s\t%s\t21\n" % (acount,x[0],flag)) ## File with abundance for all sizes merged
        fh_out2.write("%s\t%s\t%s\t22\n" % (acount,x[1],flag)) ## File with abundance for all sizes merged
        fh_out2.write("%s\t%s\t%s\t23\n" % (acount,x[2],flag)) ## File with abundance for all sizes merged
        fh_out2.write("%s\t%s\t%s\t24\n" % (acount,x[3],flag)) ## File with abundance for all sizes merged

    fh_out.close()
    fh_out2.close()

    return mergedOut,mergedsizeOut

###### Helper Functions ########################
###############################################
def chart(trans5,trans3,anormclust,bnormclust):
    '''
    This prepares line plot for pair
    '''

    print("\nFn: Chart")
    
    ## Prepare #####
    plotFile        = "%s-%s" % (trans5,trans3) 
    anormclust_s    = sorted(anormclust,key=lambda x: float(x[7]),reverse=False)
    bnormclust_s    = sorted(bnormclust,key=lambda x: float(x[7]),reverse=False)
    ### Prepare abundance and positions list for A
    #############################################
    print("Prepare list for 'A' for plotting")
    a21x = [] ## Record size specific positions
    a22x = []
    a23x = []
    a24x = []
    aOthersx = [] ## Other sizes in allowed in plotLen
    a21y = [] ## Record size specific abundances
    a22y = []
    a23y = []
    a24y = []
    aOthersy = [] ## Other sizes in allowed in plotLen
    aALLx = [] ## Collect psoitions for all sizes (sorted)
    aALLy = [] ## Collect abundances for all sizes (sorted on position)
    aphasposS  = set() ### Records all phased psoition to define x-axis
    for aphasi in anormclust_s:
        print("Ent",aphasi)
        phasipos    = int(aphasi[7])
        phasiabun   = round(math.log(int(aphasi[2]),2),2)
        phasilen    = int(aphasi[4])
        phasiflag   = aphasi[6]
        # print("Phasi Abundance:%s" % (phasiabun))
        if phasiflag == "P": ## Record the psoition to define x-axis
            aphasposS.add(phasipos)

        aALLx.append(phasipos)
        aALLy.append(phasiabun)

        if phasilen     == 21:
            a21x.append(phasipos)
            a21y.append(phasiabun)
        elif phasilen   == 22:
            a22x.append(phasipos)
            a22y.append(phasiabun)
        elif phasilen   == 23:
            a23x.append(phasipos)
            a23y.append(phasiabun)
        elif phasilen   == 24:
            a24x.append(phasipos)
            a24y.append(phasiabun)
        elif phasilen >= minPlotLen and phasilen <= maxPlotLen:
            aOthersx.append(phasipos)
            aOthersy.append(phasiabun)
        else:
            print("The sRNA length not allowed for plotting - will be trashed")
            pass
    print("Length of 21 list:%s | 22 list:%s | 23 list:%s | 24 list:%s" % (len(a21x),len(a22y),len(a23x),len(a24x)))
    print("phased positions on A:",aphasposS)


    ### Prepare abudnace and positions list for B
    ###############################################
    print("Prepare list for 'B' for plotting")
    b21x = []
    b22x = []
    b23x = []
    b24x = []
    bOthersx = [] ## Other sizes in allowed in plotLen
    b21y = []
    b22y = []
    b23y = []
    b24y = []
    bOthersy = [] ## Other sizes in allowed in plotLen
    bALLx = [] ## Collect psoitions for all sizes (sorted)
    bALLy = [] ## Collect abundances for all sizes (sorted on position)
    bphasposS  = set() ### Records all phased psoition to define x-axis
    for bphasi in bnormclust_s:
        print("Ent",bphasi)
        phasipos    = int(bphasi[7])
        phasiabun   = round(-1*(math.log(int(bphasi[2]),2)),2)
        phasilen    = int(bphasi[4])
        phasiflag   = bphasi[6]
        # print("Phasi Abundance:%s" % (phasiabun))
        if phasiflag == "P": ## Record the psoition to define x-axis
            bphasposS.add(phasipos)
            # maxabunS
        
        bALLx.append(phasipos)
        bALLy.append(phasiabun)

        if phasilen     == 21:
            b21x.append(phasipos)
            b21y.append(phasiabun)
        elif phasilen   == 22:
            b22x.append(phasipos)
            b22y.append(phasiabun)
        elif phasilen   == 23:
            b23x.append(phasipos)
            b23y.append(phasiabun)
        elif phasilen   == 24:
            b24x.append(phasipos)
            b24y.append(phasiabun)
        elif phasilen >= minPlotLen and phasilen <= maxPlotLen:
            bOthersx.append(phasipos)
            bOthersy.append(phasiabun)
        else:
            print("The sRNA length not allowed for plotting - will be trashed")
            pass
    print("Length of 21 list:%s | 22 list:%s | 23 list:%s | 24 list:%s" % (len(b21x),len(b22y),len(b23x),len(b24x)))
    print("phased positions on B:",bphasposS)



    #### Prepare chart ###############
    ##################################

    a21 = plt.scatter(a21x, a21y, s=[n**2 for n in a21y], c='#27408B', lw = 0, alpha=0.5)
    a22 = plt.scatter(a22x, a22y, s=[n**2 for n in a22y], c='#2E8B57', lw = 0, alpha=0.5)
    a23 = plt.scatter(a23x, a23y, s=[n**2 for n in a23y], c='#7D26CD', lw = 0, alpha=0.5)
    a24 = plt.scatter(a24x, a24y, s=[n**2 for n in a24y], c='#EE9A00', lw = 0, alpha=0.5)
    aOthers = plt.scatter(aOthersx, aOthersy, s=[n**2 for n in a24y], c='#383838', lw = 0, alpha=0.5)
    aAll = plt.plot(aALLx,aALLy, c='#828282', lw = 0.5)

    b21 = plt.scatter(b21x, b21y, s=[n**2 for n in b21y], c='#27408B', lw = 0, alpha=0.5) ### Size values need not be in negative like y-axis
    b22 = plt.scatter(b22x, b22y, s=[n**2 for n in b22y], c='#2E8B57', lw = 0, alpha=0.5)
    b23 = plt.scatter(b23x, b23y, s=[n**2 for n in b23y], c='#7D26CD', lw = 0, alpha=0.5)
    b24 = plt.scatter(b24x, b24y, s=[n**2 for n in b24y], c='#EE9A00', lw = 0, alpha=0.5)
    bOthers = plt.scatter(bOthersx, bOthersy, s=[n**2 for n in a24y], c='#383838', lw = 0, alpha=0.5)
    bAll = plt.plot(bALLx,bALLy, c='#828282', lw = 0.5)

    ##### Add phasing lines #############
    aphasmin    = min(aphasposS)
    aphasmax    = max(aphasposS)
    aphasx_marker = list(np.arange(aphasmin,aphasmax+phase,step=phase)) ## Phased position from first till last phase, plus an extra marker added
    bphasmin    = min(bphasposS)
    bphasmax    = max(bphasposS)
    bphasx_marker = list(np.arange(bphasmin,bphasmax+phase,step=phase)) ## Phased position from first till last phase, plus an extra marker added

    aphasy_marker = [max(aALLy)]*len(aphasx_marker) ## Most aundant value on y-axis multiples the times of x marker to get list of same size as aphasx_marker
    bphasy_marker = [min(bALLy)]*len(bphasx_marker) ## Most aundant value on y-axis multiples the times of x marker to get list of same size as bphasx_marker

    for i,j in zip(aphasx_marker,aphasy_marker):
        # plt.axvline(i,linewidth=0.5, color='#454545')
        plt.plot((i,i),(0,j),linewidth=0.3, color='#454545',linestyle='--') ## plt.plot((x1, x2), (y1, y2), 'k-')
    for x,y in zip(bphasx_marker,bphasy_marker):
        plt.plot((x,x),(0,y),linewidth=0.2, color='#454545',linestyle='--') ## plt.plot((x1, x2), (y1, y2), 'k-')


    atitle      = "%s-%s" % (trans5,trans3) 
    plt.axhline(linewidth=0.5, color='#454545')
    plt.ylabel('Abundance', fontproperties=font_manager.FontProperties(size=8))
    plt.xlabel('Fold-back normalized positions [phas-lines computed on first phasiRNA - for indicative purpose only]', fontproperties=font_manager.FontProperties(size=7))
    plt.title('%s | Phases: %s(5\') - %s(3\')' % (atitle,len(aphasx_marker)-1,len(bphasx_marker)-1),fontproperties=font_manager.FontProperties(size=9)) ## In figures 1 or 2 phases might be extra, these are noise sRNAs showing up in phass

    plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
    plt.clf() ## Clear figure so that next plot is different file

    print("\nFn: Exiting Chart")

    return None

def blastReader(blast_rc,blast_comp,blast_nor):
    '''Read all three blast results and give back a transcript'''

    fh_in = open(blast_rc,'r')
    header = fh_in.readline().strip("\n")
    read_rc = fh_in.readlines()
    fh_in.close()

    fh_in = open(blast_nor,'r')
    fh_in.readline()
    read_nor = fh_in.readlines()
    fh_in.close()

    fh_in = open(blast_comp,'r')
    fh_in.readline()
    read_comp = fh_in.readlines()
    fh_in.close()

    list_rc = []
    list_nor =[]
    list_comp =[]

    for i in read_rc:
        ent = i.strip("\n").split("\t")
        list_rc.append((ent))

    for i in read_nor:
        ent = i.strip("\n").split("\t")
        list_nor.append((ent))

    for i in read_comp:
        ent = i.strip("\n").split("\t")
        list_comp.append((ent))

    print("\nEntries in BLAST_RC:%s | BLAST_NOR:%s | BLAST_COMP:%s\n" % (len(read_rc),len(read_nor),len(read_comp)))

    return list_rc,list_nor,list_comp,header

def phasReader(phasedRes):
    '''Read phased transcripts file and provide a non-redundant set'''

    fh_in = open(phasedRes,'r')
    fh_in.readline()
    phasRead = fh_in.readlines()
    fh_in.close()

    phasList = [] ## List to store phased results
    phasSet  = set() ## Set to keep non-redundant phased transcripts
    for i in phasRead:
        name,pval,trans,start,end,strand,lib = i.strip("\n").split("\t")
        phasList.append((name,pval,trans,start,end,strand,lib))
        phasSet.add(trans)

    print("\nTotal phased loci:%s | Phased transcripts:%s" % (len(phasList),len(phasSet)))

    return phasList,phasSet

def isoformsDict(list_nor,phasSet,unassignList):
    ''' This function identifies phased transcripts that are isoforms of other phased transcripts
    and provide a non-redundant set of phased transcripts'''

    ### Make dictionary of isoforms
    isoDict = {}
    for trans in phasSet:
        # print("\nCaching isoforms for:%s" % (trans))
        isoforms = [] ### List to store trans wise isoforms
        for ent in list_nor:
            query   = ent[0]
            sub     = ent[1]
            pid     = float(ent[2])
            length  = float(ent[3]) ## If threated as integer the divison below show "truncation" and will report results in 0 or 1
            qlen    = float(ent[4])
            # print(length,qlen,(length/qlen))
            mapPerc = round(float(length/qlen),2) ##
            # print("BLAST ENTRY: %s\t%s\t%s\t%s" % (query,sub,pid,mapPerc))

            if trans.strip() == query.strip(): ## This transcript maps to another transcript, possible isoform
                if query.strip() != sub.strip():
                    if pid >= 99.0 and mapPerc >= 0.05:
                        print("+These are isoforms - Query:%s | Subject:%s" % (trans,sub))
                        print("+BLAST ENTRY: %s\t%s\t%s\t%s" % (query,sub,pid,mapPerc))
                        isoforms.append((sub,pid,mapPerc))
                    else:
                        print("-Possible isoforms- Query:%s | Subject:%s" % (trans,sub))
                        print("-BLAST ENTRY: %s\t%s\t%s\t%s" % (query,sub,pid,mapPerc))
                        pass
                else:
                    # print("-This is self match")
                    pass
            
        ## Append isoforms to dictionary for this trans
        isoDict[trans] = (isoforms) ## Key is transcripts and values are tuples with isoform,pid,mapPerc

    ### Compute how many of total phased transcripts have isoforms
    print("\nEntries in phasSet:%s | Isoform dictionary:%s" % (len(phasSet),len(isoDict)))
    acount = 0 ## count of transcripts that have isoforms
    bcount = 0 ## How many isoforms 
    ccount = 0 ## Transcripts with no isoforms
    for akey in isoDict.keys():
        avalue = isoDict[akey]
        # print("Key:%s | value:%s" % (akey,avalue))
        if avalue:
            acount+=1
            bcount+=len(avalue)
        else:
            # print("No isoforms for this transcript:%s" % (akey))
            ccount +=1
    print("Total Uniq phased trans:%s | Trans with Isoforms:%s | Trans with no Isoforms:%s\n" % (len(phasSet),acount,ccount))

    return isoDict

def getIRcoords(fastaList,trans5,trans3):
    
    '''
    This function will take pair name, fetch fasta, cat both and run inverted repeat, finally report back 
    start1,end1,start2,end2
    '''

    ## Fetch seq for trans5 and trans3 ##
    for aseq in fastaList:
        name,seq,alen = aseq
        if name == trans5:
            trans5seq = seq
            trans5len = alen
        elif name == trans3:
            trans3seq = seq
            trans3len = alen
        else:
            pass
    

    # print("+trans5:%s | len:%s" % (trans5,trans5len))
    # print("+trans3:%s | len:%s" % (trans3,trans3len))
    if trans5 != trans3: ## Two seprate transcripts that need to be combined
        if trans5seq != "" and trans3seq != "":
            combSeq     = trans5seq+trans3seq
            combName    = "%s-%s" % (trans5,trans3)

        else:
            if trans5seq == '':
                print("No seqeunces found for trans5:%s" % (trans5))
            else:
                print("No seqeunces found for trans3:%s" % (trans3))
                sys.exit()

    elif trans5 == trans3:
        if trans5seq != "":
            combSeq     = trans5seq
            combName    = "%s" % (trans5)
            trans3len   = trans5len ## Because in above loop if both trans5len and trans3len are same then only trans3len will be filled
        else:
            print("No seqeunces found self foldback trans5:%s" % (trans5))
            sys.exit()

    ## Run inverted repaet and fetch IR coords
    outseq,outinv   = IRchecker(combSeq,combName)
    IRcoords        = IRparser(outseq,outinv,combName)

    return IRcoords,trans5len,trans3len

def getPhase(aphas,clusters):

    '''
    Function to get phased position for a transcript. In guess mode - If the first sRNA
    is at positive stand then start position or sRNA position is correct phase,
    if first phasiRNA is on negative strand then sRNA position-21nt is the first
    phase correcponding to positive strand. This learned from pingchuans script
    which only uses positive strand as correct phase. In observed mode - phased positions 
    for only positive strand is recorded because because of assumption that phaiRNAs on other strand
    will map to other arm too and captured there. This avoids redundancy in phased positions.
    '''
    print("\nFUNCTION - getPhase")
    print("+Input:",aphas)

    ### Get the phased-entries for trans ###

    resList         = [] ## Store final results as (phas,[(phasiRNA),(PhasiRNA)],[extra info])
    resList2        = [] ## Store phasiRNAs from all clusters as (phas,[(phasiRNA),(PhasiRNA)],[extra info])
    
    phasCount       = 0  ## Total phased loci in file
    uniqMatchCount  = 0  ## Atleast one cluster present for one phased loci
    allMatchCount   = 0  ## Total number of matched cluster
                                                            
                                                                            
    phasID,pval,get_chr_id,get_start,get_end,trash,get_lib = aphas                        ## Given an entry in coords file
    get_value  = (list(range(int(str(get_start)),int(str(get_end)))))
    # print("This is the PhasId: %s | values:%s" % (phasID,get_value))
    print("+PhaseID being queried:%s ##############" % (phasID))
    phasCount +=1 
    # print("%s/%s phasID" % (phasCount,len(phasList)))


    ## Find matching cluster
    matchCount      = 0        ## Total maching clusters for a phased loci - if same cluster in multiple libraries
    # finalMatchList  = []       ## Holds best cluster, from multiple libraries
    finalMatchList  = []       ## Holds cluster for given PHAS locus, only one cluster is expected since we are checking for same transcript, start and stop
    tempAllList     = []     ## Hold phasiRNAs from all matching clusters, of use for phased transcripts to capture allphasiRNAs, must be used with low matchThres
    for aclust in clusters[1:]:
        tempMatchList   = [] ## To hold results of current matching cluster
        aclust_splt     = aclust.split('\n')
        header          = aclust_splt[0].split()
        clust_id        = header[2]
        chr_id          = header[6].replace("chr","").replace("Chr","")
        start           = header[10]
        end             = header[12] ##1 added because when opening a range using start and end, end number is not included in range
        value           = (list(range(int(str(start)),int(str(end)))))

        if get_chr_id == chr_id: ## Only chr_id i.e. transcript is checked to increase speed
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
                    phasistrand = str(phasient[2]).translate(str.maketrans("+-","wc"))
                    phasipval   = phasient[12]
                    # print(phasistrand)

                    # print("+phasiRNAs: %s,%s,%s,%s,%s,%s,%s"% (phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasipos))
                    tempMatchList.append((phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasipos))
                    tempAllList.append((phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasipos)) ## Records all phasiRNAs from all clusters
                    # sys.exit()

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
                    finalMatchList  = list(tempMatchList)
                    allphasiList    = list(tempMatchList) 

                    # print("\nFinal Match List:",finalMatchList)
        else:
            # print("No Match with this cluster")
            pass


    # print("Selected cluster:%s" % (finalMatchList))
    matchClust = finalMatchList[:-1]
    # print("+Found cluster:%s" % (matchClust))

    phasCoordList       = [] ## store guessed/observed phased coords
    ## Guess phased coords using the first phasiRNA position
    if phasMode == 'G':
    ## Make phase coords list
        firstPhasiStrand    = matchClust[0][1]
        firstPhasiPos       = matchClust[0][6]
        
        if firstPhasiStrand == 'w':
            firstphase = firstPhasiPos ## This is same as get_start position from phased file or start from cluster file
            print ("+First phase:%s | Phas end:%s " % (firstphase,str(get_end)))
            while firstphase <= int(get_end):
                phasCoordList.append(firstphase)
                firstphase+=phase
                # print("Current phase:",firstphase)

        elif firstPhasiStrand == 'c':
            firstphase = firstPhasiPos-(phase-3) ## Pingchuan substract 21nt if first phasi is from reverse strand
            print ("+First phase:%s | Phas end:%s " % (firstphase,str(get_end)))
            while firstphase <= int(get_end):
                phasCoordList.append(firstphase)
                firstphase+=phase

        else:
            print("+phasiRNA strand is not recongized:%s" % (phasistrand))
            sys.exit()

        print("+This is the guessed phased coords list:%s" % (phasCoordList))

    elif phasMode == "O":
        ## Fetch positions of phasiRNAs on 'w' strand, because in case of inverted repeat, 'c' strand phasiRNAs will
        ## also map to other arm
        for phasi in matchClust:
            aphasistrand = phasi[1] 
            if aphasistrand == 'w': ## These cordinates are used for validation and for pupose of simlicity on 'w' coords are used. For chart file both 'w' and 'c' phasiRNAs must be used
                aphasipos = phasi[6]
                if aphasipos not in phasCoordList:  ## To avoid adding same position twice, I noticed one case where same position appreared twice in one cluster
                                                    ## Also if there are multiple clusters from different libraries then phasiRNA positions from all included 
                    phasCoordList.append(aphasipos)
            else:
                # print("phasiRNA on other strand, must be mapped to other IR arm")
                pass

        print("+This is the observed phased coords list:%s" % (phasCoordList))

    return phasCoordList,matchClust            

def fastaReader(fastaFile):
    
    '''Cleans FASTA file - multi-line fasta to single line, header clean, empty lines removal'''

    print("\nFUNCTION - fastaReader")
    ## Read seqeunce file
    print ('-Reading "%s" FASTA file' % (fastaFile))
    fh_in       = open(fastaFile, 'r')
    fasta       = fh_in.read()
    fasta_splt  = fasta.split('>')
    acount      = 0 ## count the number of entries
    empty_count = 0

    fastaList = [] ## Stores name and seq for fastFile
    fastaDict = {}

    acount +=1
    for i in fasta_splt[1:]:
        acount  +=1
        ent     = i.split('\n')
        name    = ent[0].split()[0].strip()
        seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        alen    = len(seq)
        fastaList.append((name,seq,alen))
        fastaDict[name] = seq

    print("-Total entries in phased fastaFile:%s" % (str(acount)))
    print("-fastaList generated with %s entries\n" % (str(len(fastaList))))
    print("-Length for %s:%s" % (name,alen))

    return fastaList,fastaDict

def IRchecker(combSeq,combName):
    '''
    This function takes a seqeunce, and runs einverted and reports back
    coordinates
    '''

    print("\nFUNCTION - IRchecker")
    seq         = combSeq
    name        = combName
    tempInput   = "tempSeq.fa"
    fh_out      = open(tempInput,'w')
    fh_out.write('>%s\n%s\n' % (name,seq))
    fh_out.close()

    outseq = "%s.fa.temp" % name
    outinv = "%s.inv.temp" % name

    retcode = subprocess.call(["einverted", "-sequence", tempInput, "-gap", str(gap), "-threshold", str(threshold), "-match",str(match),"-mismatch", str(mismatch), "-maxrepeat",str(maxRepLen), "-outfile",outinv, "-outseq",outseq ])

    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        print('\n+einverted for %s\n' % (name) )
    else:
        print('Something wrong happened while running einverted for sequence: %s - - Debug for reason' % (name))
        sys.exit()

    ## Cleanup entry specifc FASTA file
    if os.path.exists(tempInput):
        os.remove(tempInput)

    return outseq,outinv

def IRparser(outseq,outinv,combName):
    '''
    parse results and delete input files
    '''
    print("\nFUNCTION - IRparser")
    print ("+Parsing %s results" % (outinv))

    ## File to record results
    # resOut= "RES_%s.csv" % (combName)
    # fh_out = open(resOut,'w')
    IRcoords = [] ## Store coordinates for this inverted repeat
    # fh_out.write("EntryName,Score,Matches,Perc,Gaps,AlignLen,5'Start,5'End,3'start,3'end,Loop\n")

    ## Get list of files
    afile = outinv
    name = combName
    # print("Seqeunce being parsed:%s" % (afile))
    try:
        if os.stat(afile).st_size > 0:
            print ("+Sequence %s - Results !!!!" % (afile.split(".")[0]))
            fh_in = open(afile,'r')
            invs = fh_in.read().split("\n\n")
            # print("Empty line splitted:",invs)
            
            for i in invs:
                # print('\nInverted:',i.strip('\n'))
                invLines = i.strip('\n').split('\n')
                # print(invLines)
                
                resBlock_splt = invLines[0].split(":")
                # print(resBlock_splt)
                score = resBlock_splt[1].split()[1]
                # print(resBlock_splt[1].split())
                matchesInfo,gapsInfo = resBlock_splt[2].split(",")
                # print(matchesInfo,gapsInfo)
                gaps = gapsInfo.split()[0]
                
                # print(matchesInfo.strip().split())

                matches = matchesInfo.strip().split()[0]
                matched,total = matches.split("/")
                alignLen = int(total)*2
                perc = round(int(matched)/int(total),2)

                # matches,garbage1,perc,garbage2 = matchesInfo.strip().split()
                # alignLen = int(matches.split("/")[1])*2
                # print(score,matches,gaps)

                arm5 = invLines[1].strip()
                start5,seq5,end5 = arm5.split(" ")

                arm3 = invLines[3].strip()
                end3,seq3,start3 = arm3.split(" ")

                loop = int(start3)-int(end5)

                IRcoords.append((name,score,matches,str(perc),gaps,str(alignLen),start5,end5,start3,end3,loop))
                # fh_out.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (afile.split(".")[0],score,matches,str(perc),gaps,str(alignLen),start5,end5,start3,end3,loop)) ## name of file,score,match,gaps

            # Clean up ### 
            garbage = [afile for afile in os.listdir('./') if afile.endswith (('.fa.temp','.inv.temp'))]
            for afile in garbage:
                print("+Deleting %s" % (afile))
                os.remove(afile)

        else:
            # print ("Sequence %s - No results" % (afile.split(".")[0]))
            pass

    except OSError:
        print("+No result file for sequence %s found - Please check" % (afile))
        print("+System will exit")
        sys.exit()


    return IRcoords

def checkOverhang(atrimcoords,btrimcoords):
    '''
    This function checks for overhang of phased coords by selecting either one of arm
    and comapring against phased coords from another arm - expected phase overhang 
    is either 2nt or (phase-2nt) i.e 19nt in case of 21mers or 21nt in case of 24mers.
    '''

    print("\nFUNCTION - checkOverhang")
    finalstatus = '' ## This will hold status of this pair P: phased or non-phased
    dominantstatus = [] ## Most dominant overhangcomputed for all 5' coords against all 3' coords

    ## There is considerable overlap, and any arm can used as query list - This scenario best suits case of guessed coordinates
    # if len(atrimcoords) >=3 and len(btrimcoords)>=3:
    if phasMode == 'G':
        ## Since both lists are of almost same size, no need to infer longer and shorter list, like
        ## required in observed mode below, so just rename lists to use with final test below
        xtrimcoordsList = list(atrimcoords)
        ytrimcoordsList = list(btrimcoords)

    elif phasMode == 'O':
        ## This scenario is for cases when real phasRNA coords are used instead of guessed phased coords
        ## that span from start to end of transcripts. 
        
        print("+Observed phasiRNA position will be used to infer overhang between arms ") ## Possible in case of real phases inferred from 24nt phasiRNA position, unlike guessed phases used above
        
        ### Choose shorter list and iterate it over longer list to find mode(overhang
        acoordslen          = len(atrimcoords)
        bcoordslen          = len(btrimcoords)
        coordsList          = [(atrimcoords,acoordslen),(btrimcoords,bcoordslen)]
        coordsList_sorted   = sorted(coordsList, key=operator.itemgetter(1),reverse=False)
        # print(coordsList_sorted)
        xtrimcoordsList     = coordsList_sorted[0][0] ## Shorter list - ([145, 187], 2) - Position and counts of position
        ytrimcoordsList     = coordsList_sorted[1][0] ## Longer list - ([121, 142, 184, 205, 226, 247], 6)] - Position and counts of position
        print("+This is the smaller list:%s" % (xtrimcoordsList))

    else:
        print("Choose correct mode to infer phasiRNA positions from combined cluster file")
        print("Available options - 'G' for guessed mode and 'O' for observed positions")
        sys.exit()


    ## Now use short and long phasiRNA position list to infer dominant status phasiRNAs on both arms i.e. 2nt overhang or not
    ###########################################################################################
    for x in xtrimcoordsList:
        tempstatus = [] ## List to store temporary status (phased or unphased) of one coords agains all coords on other strand
        
        for y in ytrimcoordsList:
            z = abs(abs(x)-abs(y)) ## 142-145 = 3nt - 1 (to avoid counting end/start of other phasi) = 2nt overhang

            if z > phase: ## Find distance greater then phase under investigation, then get reminder after dividng by phase len
                overhang = z % phase

                if overhang == 3 or overhang == phase-3: ## Overhang of 3nt = 2nt overhang with overlapping phasiRNA on other strand +1nt for start of other phasiRNA, because we are dealing with coordinates where sRNAs start, there is no end co-ordinate here. Practically, overlap should be phase-2nt, but we do not include 1nt lost in substraction, for eg 5' has phasiRNA mapped from 187 to 207 = 21nt, and 3' had overlapping phasiRNAs mapped from 205 to 185 = 21nt, but theoratical overlap b/w them using start paosition is 187-205 = 18nt, where as practical overlap is 19nt 
                    ## Perfectly phased
                    print("+coord:%s | coord:%s | overhang:%s | status = IR" % (x,y,overhang))
                    tempstatus.append('IR')
                
                elif overhang == 4 or overhang == 2 or overhang == phase-4 or overhang == phase-2:
                    ## Dicer offset
                    # print("Dicer got sloppy here")
                    print("+coord:%s | coord:%s | overhang:%s | status = IR - Dicer offset" % (x,y,overhang))
                    tempstatus.append('IR')
                
                else:
                    ## Out of phase
                    print("+coord:%s | coord:%s | overhang:%s | status = Non IR" % (x,y,overhang))
                    tempstatus.append('No')
  
            elif z < phase:
                ## Investigate the reminder
                if z == 3 or z == phase-3:
                    ## Perfectly phased
                    print("+coord:%s | coord:%s | overhang:%s | status = IR" % (x,y,z))
                    tempstatus.append('IR')
            
                elif z == 4 or z == 2 or z == phase-4 or z == phase-2:
                    ## Dicer offset
                    print("+coord:%s | coord:%s | overhang:%s | status = IR - Dicer offset" % (x,y,z))
                    tempstatus.append('IR')

                else:
                    ## Out of phase
                    # print("+Strange phase length:%s" % (z))
                    print("+coord:%s | coord:%s | overhang:%s | status = Non IR - weird phase" % (x,y,z))
                    tempstatus.append('No')
                    # time.sleep(1)


        statuscounts    = Counter(tempstatus) ## Create a dictionary/index of counts with most counts at the begining

        ## If two different status
        if len(statuscounts) > 1:
            status1,status2 = statuscounts.most_common(2) ## Top two most common overhang, e.x - [('atul', 2), ('setu', 2)], value and their counts, always ordered in decresing order
            print("+Two different local status reported:",status1,status2)
            print(status1[0],status1[1],status2[0],status2[1]) ## ('NP', 8) ('P', 7)
            # print(status1[1],status2[1])
            
            if int(status1[1]) > int(status2[1]): ## Which status (P or NP) has more counts
                print("+Clear local winner:%s" % (status1[0]))
                dominantstatus.append(status1[0])
                
            elif int(status1[1]) == int(status2[1]): ## If Tie, add both to dominant list 
                print("+Two overhang with similar count:%s and %s" % (status1[0],status2[0]))
                dominantstatus.append(status1[0])
                dominantstatus.append(status2[0])


        else:
            # print("+There one clear status for this query coord")
            status1 = statuscounts.most_common(1)
            astatus = status1[0]
            print("Single local status",astatus)
            # print(status1[0][0],status1[0][1])
            print("Clear local mandate:%s" % (astatus[0])) ## [('P', 6)]
            dominantstatus.append(astatus[0])
            

    ### Finally, check dominant status for every coord of query arm and reprt final status
    # print(dominantstatus)
    domstatuscounts         = Counter(dominantstatus)
    print("Dominant Status:",domstatuscounts)


    if len(domstatuscounts) > 1:
        domstatus1,domstatus2   = domstatuscounts.most_common(2)
        astatus = domstatus1 ## Tested OK
        bstatus = domstatus2 ## Tested OK
        print("++Mixed mandatae - adomstatus:%s | bdomstatus:%s" % (astatus,bstatus))
        
        if astatus[0] == 'IR' or bstatus[0] == 'IR':
            print("+The pair is true IR-based phased")
            finalstatus = 'IR'
            # time.sleep(1)
        else:
            finalstatus = 'No'

    else:
        domstatus1      = domstatuscounts.most_common(1)
        astatus         = domstatus1[0] ## Tested OK
        print("++ Clear mandate:",(astatus))
        
        if astatus[0] == 'IR':
            print("+The pair is true IR-based phased")
            finalstatus = 'IR'
        else:
            print("+This is not phased pair")
            finalstatus = 'No'

    return finalstatus

def phasedClusterWriter(trans5,trans3,aclust,bclust,fh_out1):
    
    '''
    This functions writes phased cluster file which is required to 
    prepare input for prepareChartFiles
    '''
    if trans5 != trans3: ## Two transcipts, not a direct IR
        ## Write trans5 coords
        fh_out1.write(">%s\n" % (trans5))
        for phasi5 in aclust:
            # print("phasi from 5arm:",phasi5)
            fh_out1.write("%s\n" % ('\t'.join(str(i) for i in phasi5)))    
        ## write trans3 coords 
        fh_out1.write(">%s\n" % (trans3))
        for phasi3 in bclust:
            # print("phasi from 3arm:",phasi3)
            fh_out1.write("%s\n" % ('\t'.join(str(i) for i in phasi3)))

    elif trans5 == trans3: ## Self foldback i.e. direct IR
        fh_out1.write(">%s\n" % (trans5))
        ## Write 5 arm coords
        writtenTags = []  ## If the tag is from same cluster, then it doesn't need to be written twice
        
        for phasi5 in aclust: ## aclust could be same or different for blcust - depends if phasing was close or far from each other on same transcript
            phasiseq = phasi5[3]
            phasipos = phasi5[6]
            uniqPos  ='%s-%s' % (phasiseq,phasipos) 
            writtenTags.append(uniqPos)
            print("phasi from 5arm:",phasi5)
            fh_out1.write("%s\n" % ('\t'.join(str(i) for i in phasi5)))
        
        for phasi3 in bclust:
            phasiseq = phasi3[3]
            phaspos  = phasi3[6]
            uniqPos  ='%s-%s' % (phasiseq,phasipos) 

            if uniqPos not in writtenTags:
                print("phasi from 3arm:",phasi3)
                fh_out1.write("%s\n" % ('\t'.join(str(i) for i in phasi3)))
    else:
        print("Check Trans5:%s and Trans3:%s names - You should have never reached here" % (trans5,trans3))
        pass

    return None

def TagAbundanceFile(con,db,libs):
    
    print("\nFUNCTION: TagAbundanceFile")
    
    fetchedLibs = []
    for alib in libs:##For all the libraries
        fetchedLibs.append(alib)
        
        ## Check if file already exsits in directory - This saves a lot of time downloading the same file
        filePath = '%s.fas' % (alib)
        if os.path.isfile(filePath) == False:
            print ('\nPreparing sRNA reads file for library: %s' % (alib[0]))
            #print (lib[0])
            #print ('Caching tag and count information from server for PARE alib %s' % (alib[0]) )
            cur = con.cursor()
            cur.execute("SELECT tag,norm from %s.run_master where lib_id = %s AND (hits between 0 and 20)" % (db,alib[0]))
            lib_info = cur.fetchall()
            #print('These are the tags:',lib_info[:10])
            
            fh_out = open('%s.fas' % (alib), 'w')##Naming file with lib_ids name
            print ('Library cached, writing abundance file')
            tag_num = 1
            for ent in lib_info:## All the PARE tags in a library
                #print (ent)
                fh_out.write('%s\t%s\n' % (ent[0],ent[1]))
                tag_num += 1
                
            fh_out.close()
        else:
            print('tag abundance file exists for library: %s' % (alib))
            pass
    
    return fetchedLibs

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

def finalPairs(resFile):
    '''
    parser for final pairs file, to enable direct run of 
    prepareChartFiles
    '''
    print("\nFUNCTION: finalPairs")
    fh_in = open(resFile,'r')
    fh_in.readline()
    fileRead = fh_in.readlines()

    validList = [] ## List to store final pair results in same format as came directly from script([blastent],totalPhasi,num avalidcoords,num bvalidcoords, finalstatus)
    for i in fileRead:
        # print(i)
        ent = i.split("\t")
        blastent        = ent[:22]
        totalPhasi      = ent[23]
        nvalidcoordsA   = ent[24]
        nvalidcoordsB   = ent[25]
        finalstatus     = ent[26].strip("\n")
        # print("This is blastent:",blastent)
        # print(totalPhasi,nvalidcoordsA,nvalidcoordsB,finalstatus)
        validList.append((blastent,totalPhasi,nvalidcoordsA,nvalidcoordsB,finalstatus))

    print("+Result list with %s final pairs prepared from 'finalPairs'" % (str(len(validList))))

    return validList

def tagCount2FASTA(inFile,Exprs):
    
    print("\nFUNCTION:tagCount2FASTA")

    fh_in=open(inFile, 'r')
    outFile = '%s.fa' % (inFile.rpartition('.')[0])
    fh_out =open(outFile, 'w')
    tag_num = 1 ### For naming tags

    if Exprs=='Y':  ### Write as raw sequencing file with tag repeate dnumber of times it appears in tag_count 
        ##Write to file
        print('\nWriting expression file for %s tagcount file' % (inp_file_name))
        print('\n---PLEASE BE PATIENT---')
        
        for ent in fh_in:##All the entries of the library
            #if len(ent[0]) == 20:
            ent = ent.split('\t')
            tag_count = int(ent[1])
            for count in range(tag_count):##Number of times the tag_count file
                fh_out.write('>Tag%s\n%s\n' % (tag_num, ent[0]))
                tag_num += 1
                
    else: ##Convert tag count to FASTA
        for i in fh_in:
            ent = i.strip('\n').split('\t')
            #print(ent)
            fh_out.write('>Tag%s_%s\n%s\n' % (tag_num,ent[1],ent[0]))
            tag_num += 1
            
    fh_in.close()
    fh_out.close()
    
    return outFile

def mapLibs(fastaFile,fetchedLibs):
    
    '''
    This module prepares bowtie maps for all the libraries
    '''

    print("\nFUNCTION: mapLibs")
    ## Prepare - Make index, and filenames
    fastaIndex = "fasta.index"
    retcode = subprocess.call(["bowtie-build",fastaFile,fastaIndex])
    if retcode == 0:
        print("+Index file prepared for %s" % (fastaFile))
    else:
        print("-There is some problem in preparing index")
        sys.exit()

    for alib in fetchedLibs:
        print("+Mapping library:%s" % (alib))
        inFile = '%s.fas' % (alib)
        print ('+Processing %s for mapping to genome' % (inFile))
        fastaFile = tagCount2FASTA(inFile,'N') ## Unique reads to FASTA format 

        mapFile = ('./%s.map' % (alib))
        print(fastaIndex,inFile,fastaFile,mapFile)

        ## Start mapping 
        print ('Mapping %s processed file to genome' % (alib))
        nproc2 = str(nproc)
        mismat = str(0)

        retcode = subprocess.call(["bowtie","-f","-n",mismat,"-p", nproc2,"-t" ,fastaIndex, fastaFile, mapFile])
        
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('+Bowtie mapping for %s complete' % (inFile) )
        else:
            print ("-There is some problem with mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()

    return None

def map2Dict(alib):
    '''
    parse the bowtie map, it assumes that map file have been generated already
    '''

    print("\nFUNCTION: map2Dict")

    ## Rad sRNA map file
    mapFile = './%s.map' % (alib)
    fh_in   = open(mapFile,'r')
    mapRead = fh_in.readlines()
    fh_in.close()

    ## Prepare list of mapped sRNAs, uniq trans names before dictionary prepration
    mapList     = [] ## List of all mapped sRNAs
    transSet    = set() ## Set of uniq transcripts name
    for i in mapRead:
        ent     = i.strip('\n').split("\t")
        aname,astrand,atrans,apos,phasiseq,trash1,hits,trash2, = ent ## Last column is Comma-separated list of mismatch descriptors. 
        phasiname,phasiabun = aname.split("_")
        phasistrand         = str(astrand).translate(str.maketrans("+-","wc"))
        phasilen            = int(len(phasiseq.strip()))
        phasihits           = int(hits)+1 ## This is not a good poxy of sRNA hits, more processing required to get hits, so this number can't be trusted
        phasiflag           = 'S' ## sRNA from bowtie map file
        phasipos            = int(apos)+1 ## Convert 0-based offset to 1-based offset like Pingchuans cluster file
        transname           = atrans.strip()
        transSet.add(transname)
        if phasistrand == 'w' and phasilen >= minNoiseLen and phasilen <= maxNoiseLen: ## sRNAs of only positive strand used - OK
            mapList.append((aname,phasistrand,transname,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos))
    
    mapList_s = sorted(mapList,key=operator.itemgetter(2))
    # print("-Snippet of sorted mapped list:",[mapList_s[:50]]) ## Tested OK
    print("-List of lib:%s sRNAs prepared with %s entries" % (alib,len(mapList)))
    print("-Total transcipts recorded:%s" % (str(len(transSet))))

    ### Now prepare the sRNA dictionary
    print("-Preparing dictionaries of sRNAs")
    
    srnaDict = {}
    for atrans in transSet:
        value       = []    ## Store values for key
        matchFlag   = 0     ## This jumps the loop once match is done
        # print("\n-Caching sRNAs for trans:%s" % (atrans))
        
        for i in mapList_s: ## For entry in mapped list
            # print("-sRNA entry:",i)
            aname,phasistrand,transname,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos = i
            # print(i)
            
            if transname == atrans and matchFlag == 0: ## This is the first phasiRNA for this transcript
                matchFlag = 1 ## Match detected, now followinf sRNAs belong to this transcripts
                value.append((phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos))

            elif transname == atrans and matchFlag == 1: ## These are phasiRNAs for this transcript, other then the first one
                value.append((phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos))

            elif transname != atrans and matchFlag == 1: ## phasiRNAs for this transcripts have been cached already, jump to next transcript
                # print("-phasiRNAs already cached for trans:%s - Jumping the loop" % (transname))
                continue

            else:
                # print("phasiRNAs from %s do not belong transcript %s being cached" % (transname,atrans))
                pass


        ## Add to dictionary 
        # print("Key:%s | values:%s" % (atrans,value))
        srnaDict[atrans] = value

    print("-Dictionary with entries for %s transcipts prepared" % (str(len(srnaDict))))

    return srnaDict

def clust2Dict(clustFile):
    '''
    This module reads cluster file and gives back a dictionary of coords. In case of 
    directIRs only one phased cluster, that had valid coords is expected with phasiRNAs
    uniuqe on seq and positions.
    '''

    print("\nFUNCTION: clust2Dict")
    # clustFile = 'phased.clust'
    fh_in       = open(clustFile,'r')
    clustRead   = fh_in.read().split('>')
    fh_in.close()

    transSet = set() ## Empty set for phased transcripts
    
    ## Make a set of transcript names
    for i in clustRead [1:]:
        aclust      = i.split("\n")
        transname   = aclust[0]
        transSet.add(transname)
        print(transname)
    print("-Trascripts set prepared with %s entries" % (len(transSet)))
    # sys.exit()

    ## Search clusters with transcript names, and append their phasiRNAs
    phasedDict  = {} ## Transcripts wise phasiRNAs
    for gettrans in transSet:
        value       = [] ### List to hold phaseRNAs for a transcript
        print("-\nCaching phasiRNAs for transcript:%s" % (gettrans))
        # phasiset    = set() ## Set of uniq phasiRNAs for this transcript
        
        for ablock in clustRead [1:]:
            aclust      = ablock.split("\n")
            transname   = aclust[0]

            if transname == gettrans:
                print("-Matching block found with transcript:%s" % (transname))
                # print(ablock)
                
                for x in aclust[1:-1]: ## First entry was header and last is always empty
                    ent = x.split("\t")
                    # print("-ent:%s" % (ent))
                    if ent:
                        # print(ent)
                        phasiname   = ent[0]
                        phasistrand = ent[1]
                        phasiabun   = int(ent[2])
                        phasiseq    = ent[3]
                        phasilen    = ent[4]
                        phasihits   = ent[5]
                        phasipos    = ent[6]
                        phasiflag   = 'P' ## Phased from cluster file
                        if phasistrand == "w" or phasistrand == "c": ## Only one strand is used because these are transcripts and not genome, and in case of IR-based, 'c' strand maps to other arm so kind of redundant
                            value.append((phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos))
                            # phasiset.add(phasiseq)
            else:
                # print("This cluster corresponds to different transcipt %s" % (transname))
                pass

        print("-%s phasiRNA cached for transcript:%s" % (len(value),gettrans))
        ## Results from multiple libraries will add phasiRNAs, 
        ## Filter most abundant uniq phasiRNAs. and add to dictionary
        # value_s = sorted(value, key=operator.itemgetter(2),reverse=True) ## Sorted on Abundance from high to low
        value_s = sorted(value, key=operator.itemgetter(1),reverse=True) ## Sorted on strand from w to c, we will fetch abundances later
        # print("\n-Snippet of abundance sorted phasiRNA list:",value_s)
        
        checklist   = [] ## List to check if most abundantinstance of this phasiRNA have been recorded
        valueuniq   = [] ## List to hold uniq phasiRNAs from multiple libraries, sorted on abundance
        for anent in value_s:
            aphasiseq = anent[3]
            aphasilen = int(anent[4])
            # print("-phasiRNA beind added",anent)
            if aphasilen == phase and aphasiseq not in checklist:
                checklist.append(aphasiseq) ## add to checklist, which means considered once
                valueuniq.append(anent)
            else:
                # print("This phasiRNA %s has been recorded before" % (aphasiseq))
                pass

        ## Finally add unique values to dict        
        print("-\ntrans name:%s | value:%s" % (gettrans,valueuniq))
        phasedDict[gettrans] = valueuniq


    print("-phasedDict with %s elements prepared\n" % (str(len(phasedDict))))
    # print("Length of dictionary for 'TRINITY_DN279522_c0_g1_i2' is %s and values:" % (len(phasedDict['TRINITY_DN279522_c0_g1_i2'])))
    # print("And values:",phasedDict['TRINITY_DN279522_c0_g1_i2'])
    # sys.exit()

    return phasedDict

def mergeRNAs(trans,sRNADict,phasedDict,fastaDict,cur,finalLibs):
    '''
    This function merges sRNAs and phasiRNAs for trans5 and same for trans3. In case of direct IRs
    phasiRNAs would be uniq in seq and position will be retained. And for 'c' strand mapping for which no 'w' mapping found in cluster (clust2dict), it
    will be retrived from sRNA mappings
    '''

    print("\nFUNCTION:mergeRNAs")

    ## Output ##
    transClust = [] ## Store merged sRNAs and phasiRNA entries

    ## Prepare trans5/3 cluster ###########################################
    ######################################################################
    
    #### Fetch phasiRNAs and sRNAs, and prepare a reference list of sRNA tags
    print("+Preparing sRNAs for trans: %s" % (trans))
    phasivalues = phasedDict[trans]
    try:
        sRNAvalues  = sRNADict[trans]
    except KeyError:
        sRNAvalues = [] ## This is possible if no sRNAs from specific library mapped to this transcript, in that case there are no siRNA, it would be best to avoid such libs,
                        ## because on phasiRNAs from phased lcuster will be plotted
    print("+Number of trans phasiRNA:%s | srnas:%s" % (len(phasivalues),len(sRNAvalues)))

    localsrnadict = {} ## dict of sRNA tags and abun for comaprison with phasiRNAs - This will make them unique here
    for i in sRNAvalues:
        aseq    = i[3]
        aabun   = i[2]
        localsrnadict[aseq] = aabun
    print("+Temp sRNA abundance dictionary prepared with values:%s" % ((str(len(localsrnadict)))))

    #### Fetch phasiRNAs abundance, find 'w' psoitions for 'c' mapped phasiRNAs and prepare list to substract with sRNAs
    ####################################################################################################################
    phasiList = [] ## List to store phased tags, so that sRNAs can be filtered out that match to phasiRNAs
    for ent in phasivalues:
        print("-phasiRNAs",ent)
        phasiseq    = ent[3]
        phasistrand = ent[1]
        phasipos    = ent[7]
        ## Fetch abundances for library(s) if interest
        atag,phasilibabun,lib_abun = getAbundance(cur,phasiseq,finalLibs)
        ### Sanity check
        if not lib_abun:
            print("No abundance found for this tag:%s" % (phasiseq))
            print("Debug")
            sys.exit()
        elif int(phasilibabun) == 0:
            phasilibabun = ent[2] ## Sum of abundance in requested libraries is zero, using the abundance from phas cluster
        else:
            # print(atag,phasilibabun,lib_abun)
            pass

        
        if phasiseq in localsrnadict.keys():    ## Then record to remove from sRNA list and correct 'c' strand if present
            phasiList.append(phasiseq)          ## Add to phasilist that will be used to avoid adding srna that matched phasiRNAs
            
            ### If phasiRNA has 'c' strand then the only way to find it 'w' position by mapping sRNA to precursor
            if phasistrand == 'c':
                ## Fetch position for 'w' strand from sRNA mapping
                print("Position of %s on 'c' strand:%s" % (phasiseq,phasipos))
                aseq = fastaDict[trans]         ## get the transcript
                # print(aseq)
                startindex = aseq.find(phasiseq)## Get the start site in python '0' format
                if startindex > 0:
                    print("This is the new position on 'w':%s" % (startindex+1))
                    phasiposw       = startindex+1
                    phasistrandw    = "w" 
                    print("-Updated phasi entry:",ent[0],phasistrandw,int(phasilibabun),ent[3],ent[4],ent[5],ent[6],int(phasiposw))
                    transClust.append((ent[0],phasistrandw,int(phasilibabun),ent[3],ent[4],ent[5],ent[6],int(phasiposw))) ## phasiname,phasistrand,lib-specific phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos
                else:
                    print("Phasi not found on 'w' strand")
                    print('''This phasiRNA from 'c' strand have no 'w' strand mapping - 
                        Perhaps because the arm on 'w' strand has mismatches in loop or 
                        because this tag is from processed end where only one arm continues so 
                        this tag is from processed other arm''')
                    print("Such strands are not part of foldback and must be skipped")
                    pass
            else:
                ### Expected 'w' strand, fetch abundances for library if interest
                # print("-Updated phasi entry:",ent[0],ent[1],int(phasilibabun),ent[3],ent[4],ent[5],ent[6],int(ent[7]))
                transClust.append((ent[0],ent[1],int(phasilibabun),ent[3],ent[4],ent[5],ent[6],int(ent[7]))) ## phasiname,phasistrand,lib-specific phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos
        else:
            ### Even if phasi is not found in this library(s), it's part of PHAS locus from other library, include and fetch abundances
            phasiList.append(phasiseq)
            # print("-Updated phasi entry:",ent[0],ent[1],int(phasilibabun),ent[3],ent[4],ent[5],ent[6],int(ent[7]))
            transClust.append((ent[0],ent[1],int(phasilibabun),ent[3],ent[4],ent[5],ent[6],int(ent[7])))

    #################################################################################################################################
    ### Now add those sRNAs which are not in phasiList i.e. noise or sRNAs of different sizes from phased positions - I choose to not 
    ### fetch their abundances from multiple libraries instead use just the library-specific abudnance form file
    ##################################################################################################################################
    for ent in sRNAvalues:
        # print("sRNAs",ent)
        srnaseq = ent[3]
        if srnaseq not in phasiList:
            transClust.append((ent[0],ent[1],ent[2],ent[3],ent[4],ent[5],ent[6],int(ent[7])))
        else:
            # print("This is a phasiRNA %s added already to trans cluster" % (srnaseq))
            pass

   ### Final results
    print("+Total entries in finacluster:%s" % (len(transClust)))

    return transClust

def consensusNorm(anormclust_s,bnormclust_s,aflag,bflag):
    '''
    This function normalized the IR coords so that multiple IRs can be alinged together,
    basically it will set IR coords to 1 from left plus add buffer that equals phase
    '''
    print("\nFUNCTION: consensusNorm")
    afirstcoord     = anormclust_s[0][7] ## Sorted to phased loci at begining and then on Pos, this is start for sRNA with end further in right
    bfirstcoord     = bnormclust_s[0][7] ## Sorted to phased loci at begining and then on Pos, this is start for sRNA with end towards left
    decisionList    = [(afirstcoord,aflag),(bfirstcoord,bflag)]
    print("+++++++++++",decisionList)
    decisionList_s  = sorted(decisionList,key=operator.itemgetter(1)) ## sorted on sRNAs/phased RNAs and then on Pos
    mincoord,mincoordflag = decisionList_s[0]
    print("+++++++++++",mincoord,mincoordflag)

    if mincoordflag   == 'A':
        xnormfactor     = 1-(int(mincoord)) + phase + 3 ## Strand offset = 3nt, and buffer = phase length -OK
    elif mincoordflag == 'B':
        xnormfactor     = 1-(int(mincoord)) + phase ## This is lower strand offset = 0, and buffer = phase length - OK
    else:
        print("-This flipflag is not recognized:%s" % (mincoordflag))
        print("-Treating entry as ideal 5' and 3' pair")
        print("-Skip this entry in mergeIR")
        xnormfactor     = 1-(int(mincoord)) + phase ## This is lower strand offset = 0, and buffer = phase length - OK

    aconsensusclust = [] ## Cluster to store consensus coords for trans5
    bconsensusclust = [] ## Clusters to store consensus coords for trans3
    for phasi in anormclust_s:
        apos    = phasi[7]
        normpos = xnormfactor + int(apos) ## xnormfactor has + or - polarity, so keep calculation according to it's polarity we add apos to x norm factor, if xnormfactor is - then it will substract automatically, earlier we had apos - abs(anormfactor) which disregards it polarity hence wrong - OK
        if normpos >= 1:
            aconsensusclust.append((phasi[0],phasi[1],phasi[2],phasi[3],phasi[4],phasi[5],phasi[6],normpos)) ## phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos

    for phasi in bnormclust_s:
        apos    = phasi[7]
        normpos = xnormfactor + int(apos) ## +phase is added to include a buffer -OK
        if normpos >= 1:
            # print("+")
            bconsensusclust.append((phasi[0],phasi[1],phasi[2],phasi[3],phasi[4],phasi[5],phasi[6],normpos)) ## phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos

    return aconsensusclust,bconsensusclust

def flipflag(anormclust_s,bnormclust_s):
    '''
    This functions assigns flags (A/B) to 5' and 3'. A would mean that this arm is upper strand, and B 
    would mean that this is lower arm. So that all IR can be aligned with +/-2 nt overhang on same side
    '''
    print("\nFUNCTION: flipflag")
    afirstcoord     = anormclust_s[0][7]
    bfirstcoord     = bnormclust_s[0][7]

    ## Test ###
    overhang     = int(afirstcoord)-int(bfirstcoord) ## 5'- 3' = overhang
    print(overhang)

    if abs(overhang) >= phase:
        xfactor = int( (abs(overhang) % phase) * abs(overhang)/overhang) ## % outputs polarity of the denominator, so we alternatively retain the polarity here, becuase if crucial for test below
    else:
        xfactor = overhang

    print("+This is 5' arm snippet:%s" % (anormclust_s[:2]))
    print("+This is 3' arm snippet:%s" % (bnormclust_s[:2]))
    print("+This is afirstcoord:%s | bfirstcoords:%s | xfactor:%d" % (afirstcoord,bfirstcoord,xfactor))

    if xfactor == 3 or xfactor == -(phase-3): ## ## +3 or -18 for 21mers
        aflag = 'A' ## This should be plotted at upper side i.e. no flip
        bflag = 'B' ## This should be plotted on lower side i.e. no flip
        print("+This is aflag:%s | bflag:%s" % (aflag,bflag))

    elif xfactor == 2 or xfactor == 4 or xfactor == -(phase-2) or xfactor == -(phase-4): ## ## +2/4 or -17/-19 for 21mers
        aflag = 'A' ## This should be plotted at upper side i.e. no flip
        bflag = 'B' ## This should be plotted on lower side i.e. no flip
        print("+DicerOffset - This is aflag:%s | bflag:%s" % (aflag,bflag))
    
    elif xfactor == -3 or xfactor == phase-3: ## -3 or +18 for 21mers
        aflag = 'B' ## This should be plotted at lower side
        bflag = 'A' ## This should be plotted on upper side
        print("+This is aflag:%s | bflag:%s" % (aflag,bflag))

    elif xfactor == -2 or xfactor == -4 or xfactor == phase-2 or xfactor == phase-4: ## -3 or +18 for 21mers
        aflag = 'B' ## This should be plotted at lower side
        bflag = 'A' ## This should be plotted on upper side
        print("+DicerOffset - This is aflag:%s | bflag:%s" % (aflag,bflag))

    else:
        print("This is unrecognized xfactor to identify flipFlag:%d" % (xfactor))
        aflag = 'N'
        bflag = 'N'
        pass

    return aflag,bflag

def getAbundance(cur,tag,finalLibs):
    '''Input is tag for each loci and out put is tag with maximmum abumdance and sum of phasiRNAs - 
    rewritten in v1.0 for fetching tags from run master'''

    lib_abun = [] ## list to hold lib-wise abudnances
    
    for alib in finalLibs:
        # print("Lib:",alib)
    
        cur.execute("SELECT tag,norm FROM %s.run_master where tag = '%s' and lib_id = %s" % (db,tag,alib))### Convert intergenic to gene name so as to get strand
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
    # print("--Lib-wise abundances",lib_abun)
    # print("--Sum of abundances:%s\n" % (abun_sum))

    # sys.exit()

    return tag,abun_sum,lib_abun
            
def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### Get lib names
    columns = [] ## empty list
    cur.execute("SELECT DISTINCT(lib_id) FROM %s.run_master" % (db))
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

def main():
    if runMode == 0:
        ## Prepare
        phasList,phasSet                        = phasReader(phasedRes)
        list_rc,list_nor,list_comp,header       = blastReader(blast_rc,blast_comp,blast_nor)
        ## Infer
        uniqList,candidateSet,allSet            = inferIRs(list_rc,header,phasSet)
        validList,resfile,clustFile,pairedSet,unassignList             = validateIR(clustfile,uniqList,phasList,fastaFile,header,phasSet)
        isoDict                                 = isoformsDict(list_nor,phasSet,unassignList)
        ## Summarize
        isoPair,noIsoPair,summaryFile           = isoformChecker(phasSet,unassignList,pairedSet,isoDict,list_rc,header,phasList)
        isoformWriter(list_nor,isoPair,header,pairedSet)
        collapseNoIsoforms(isoDict,noIsoPair,header,summaryFile)
        ## Prepare chart input
        con         = ConnectToDB(server)
        fetchedLibs = TagAbundanceFile(con,db,userLibs)
        mapLibs(fastaFile,fetchedLibs)
        chartInput,consensusInput = prepareChartFiles(validList,clustFile,userLibs,fastaFile,con)
        mergedOut,mergedsizeOut = mergeIRs(consensusInput,consensusCands)

    elif runMode == 1:
        
        if fetchMap == 1:
            print("+Preparing map files - NEED TO BE DONE ONLY ONCE")
            ### Fetch Libs
            con         = ConnectToDB(server)
            fetchedLibs = TagAbundanceFile(con,db,userLibs)
            
            ### Map libs
            mapLibs(fastaFile,fetchedLibs)
            
            print("You can turn off fetchMap now")
            time.sleep(3)
            pass
        else:
            print("+Bowtie maps are supposed ro be present for libs %s specified in settings,else turn ON fetchMap" % (userLibs))
            pass

        ## Prpeare chart files
        resFile     = "finalPairs.txt"      ## This must be prepared by validateIR function in ealier run
        clustFile   = "phased.clust"        ## This must be prepared by validateIR function in ealier run
        validList   = finalPairs(resFile)   ## List of results from old run of validateIR function
        con         = ConnectToDB(server)
        chartInput,consensusInput  = prepareChartFiles(validList,clustFile,userLibs,fastaFile,con)
        mergedOut,mergedsizeOut = mergeIRs(consensusInput,consensusCands)

    elif runMode == 2:
        ## valid pairs identified already, and their chart files are ready
        ## Provided the consensusCands list, these will be merged
        for lib in userLibs:

            consensusInput  = "%s.consensus.input.txt" % (lib)
            mergedOut,mergedsizeOut = mergeIRs(consensusInput,consensusCands)

    else:
        print("Please input correct mode - Script will exit now")
        sys.exit()

if __name__ == '__main__':
    #### Assign Cores
    if numProc == 0:
        nproc = int(multiprocessing.cpu_count()*0.95)
    else:
        nproc = int(numProc)
    main()
    sys.exit()


## v02 -> v03
## Fixed isoform catching, earlier same unassigned transcript would have been added multiple times to noIsoform set exaggarating its number
## because if multiple isforms same transcript is added at every loop

## v03->v04
## Added functinality to capture isoforms to paired results
## Added functionality to write a summary file
## Added functionality to collapse no isoforms to paired

## v04 -> v05
## Improved summary for clarity

## v05 -> v06
## Added pid threshold while inferring IRs

## v06 -> v07
## Fixed a counting bug, which missed noIsopair isoforms from being counted in final summary - once they have been added to "assigned" bin. See module collapseNoIsoforms

## v07 -> v09[major]
## Added functionality to extract phasiRNA positions for both partners from combined clusters file, extract FASTA from fasta file, perform einverted analysis, and finaly map the phased position to both arms
## Ultimately decide whther candidate pair based on blast results have phasiRNAs with 2nt (3nt) offsets
## Add functionality to account for sloppy dicer in checkOverhang function
## Added functionality to catch if there are no valid coords between pairs - This can be tweaked to guess possibel phase in IR regions and then check for overlap

## v09 -> v09b[major]
## Added functionality to handle self foldback transcipts i.e. direct loop, both trans5 and trams3 is same
## In case of direct IR - einverted is run on just one transcripts
## converted coords are normal phasing coords and no length of trans5 added to trans3 coords
## Contart to planned that in case of direct-repaeat phasiRNAs on 'w' and 'c' strand are required, only 'w' strand is fetched because when checked 'c' strands are redundant to 'w' strands
#### i.e. it's reproted twice on both strands
## In addition, if multiple phase clustters in clusters file, then uniq phased locaition from all are used
## Tested the localstatus and domstatus counter slicing of status

## v09b -> v1.0[stable]
## Testing stable

## v1.0 -> 1.1 [very major]
## Added functionalty to use final pairs and generate a chart input containing phasiRNAs and other sRNAs
## Fixed couple of bugs
## Added one more condition for candidate IR, alingment length should be geater then 0.75% of smaller transcript i.e. trans5 or trans3

## v1.1 -> v1.2[stable]
## Fixed running bugs i.e. tested on 24Phased transcripts
## Modified summary variables, moved them from infer IR to validated IR
## Re-organized the code for clarity

## v1.2 -> v1.3b [major]
## Fixed sRNA parsed from map file positions by converting 0-based offset to 1-based offset
## Fixed coding bugs and major normlization bug in consensusNorm function
## Fixed IRmerge issue to consider positions till last phased siRNAs and not include sRNAs beyond that

## v1.3b - v1.4 [major][stable]
## flipflags for a few couldn't be determined, these are assigned 'N' flags and normalized for merging considering
## trans5 as A and trans3 as B - Which can be in accurate, so must not be used for merging
## These entries with "N" flipFlags are not used for merging in mergeIRs function and skipped
## Added a new mode, if final pairs have been detcted, abd chart input file is done by preparechartInput function,
## then runMode == 2, allows to just input the candidatelist for which consensus report is required. Saves a hell
## lot of time 
## Improved map2Dict for faster performance - now a list is prepared before hand of only valid sRNAs, and later when
## transcript specific iteration is done, its terminated once phasiRNAs have been cached, i.e. search doesnt goes through
## complete list
## Added size specific abundance merged results file in mergeIR module

## v1.4 -> v1.5
## Rewritten getclust module, because earlier it was itenerating over cluster file for a transcipt, and adding phasiRNAs
## from block as value to dictioary. But cluster file is made up of multiple libs and therefore it has multiple lib-
## specific clusters of same transcript. In such case, dictionary was being updated with last cluster for this transcript
## which could from any library. Now for every transcript cluster file is scanned and phasiRNA sorted on abundace and uniq.

## v1.5 -> 1.9[major][stable]
## Added fuzzy phasiRNA extract for preparing chart input and connsensus input. The validation part still uses phasiRNAs for specific loci and not fuzzy matching
## Added functionality to fetch abudnaces for phasiRNAs for libraries of interest
## Added a chart module to create pair-specific abundnace charts
## Added functionality to include other noise sizes through min/maxNoiseLen paerameter added to settings


## V1.9 -> v2.0
## Added functionality to define thelength used for plotting by min/maxPlotLen parameters, and added them to plots as grey color


### TO DO
### 1. Fix bug where a few noIsoPair transcripts were not being written to "noIsoToPairedClusts.txt" file. Problem seems in evaluating length of clust (>1), if its one list or one entry, in
### 'collapseNoIsoforms' module. For now these misisng ones can be treated as unpaired with one isoforms.