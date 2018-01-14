#!/usr/bin/python

import os,sys,operator

### Three BLAST results required from transRNA.sh

blast_rc    = "Final_BLAST_RC.txt"
blast_nor   = "Final_BLAST_NOR.txt"
blast_comp  = "Final_BLAST_COMP.txt"
phasedRes   = "Final_PHAS_Loci_5e-07_ALL.csv"

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

def inferIRs(list_rc,header,phasSet):


    fh_out = open("uniqPairs.txt", 'w')
    fh_out.write("%s\ttype\n" % (header))

    rc_sort = sorted(list_rc,key=lambda x: float(x[14]),reverse=True) ## SOrted on bitscore

    for i in rc_sort[0:5]:
        print("Sorted example",i)
    
    ## Retain best pair
    list5 = [] ## Store assigned 5' pairs
    list3 = [] ## Store assgined 3' pairs
    pairedList = [] ## Store both transcripts of pair
    uniqList = [] ## Store uniq results

    set5 = set() ## Unique set of transcripts
    set3 = set() ## Unique set of those matched

    acount = 0  ## Total entries count
    for i in rc_sort:
        trans5      = i[0]
        trans3      = i[1]
        pid         = i[2]
        length      = i[3]
        bitscore    = i[14]
        hang5       = i[15]
        hang3       = i[14]
        match       = i[16]

        # print("Interaction:%s-%s | PID:%s | Bitscore:%s | AlignLen:%s" % (trans5,trans3,pid,bitscore,length))

        set5.add(trans5)
        set3.add(trans3)
        acount+=1

        if (trans5 not in pairedList) and (trans3 not in pairedList) and int(length) >= 180:
            uniqList.append((trans5,trans3,pid,bitscore,hang5,hang3,match,"IR",i))
            fh_out.write("%s\tIR\n" % ("\t".join(z for z in i)) )

            pairedList.append(trans5)
            pairedList.append(trans3)
        else:
            # print("either hang5 or hang 3 have a partner assigned - Better pair exists")
            pass

    print("\nTotal BLAST pairs:%s | Uniq inferred IR pairs:%s" % (acount,len(uniqList)))
    fh_out.close()

    ## Results
    allSet = set5.union(set3) ## Total unique phased transcripts in BLAST Results
    pairedSet = set(pairedList) ## All transcripts that have a pair
    print("Total phased transcripts:%s | Uniq transcripts in BLAST_RC:%s | Paired Transcipts:%s | Uniq paired Transcripts:%s" % (len(phasSet),len(allSet),len(pairedList),len(pairedSet)))

    ## Identify unpaired from total phased transcripts
    unassignList = [] ## List to store those transcripts that are missing from pair-test or didn't had unoq pair
    for i in phasSet:
        if i not in pairedSet:
            unassignList.append(i)
    print("Total phased transcripts:%s | Unassigned transcripts:%s" % (len(phasSet),len(unassignList)))

    return uniqList,pairedSet,allSet,unassignList

def isoformChecker(phasSet,unassignList,pairedSet,isoDict,list_rc,header,phasList):

    '''Checks if unassigned transcripts has paired isoform'''

    noIsoPair = []  ## Final results storing those transcripts that neither have a confident pair phase transcripts not any of their isoform
                    ## Non-IR transcripts - These are most probably

    wcount = 0 ## Count unassigned
    xcount = 0 ## Unassigned with some isoform transcript
    ycount = 0 ## Unassigned transcripts that have a paired isoform
    zcount = 0 ## Unassigned transcripts that have no isoform at all
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
    fh_out = open("candidateNonIRs.txt", 'w') ## These are transcriptsthe do not have any pair neither their isoforms have a pair ones that 
    fh_out2 = open("nonIRs.list", 'w') ## These are unassigned ones that simply do not have BLAST_RC results i.e. no possible IR pair

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

    fh_out.close()
    fh_out2.close()

    return None

def isoformsDict(list_nor,phasSet,allSet,unassignList):
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
                    if pid >= 95.0 and mapPerc >= 0.05:
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

def main():
    phasList,phasSet                        = phasReader(phasedRes)
    list_rc,list_nor,list_comp,header       = blastReader(blast_rc,blast_comp,blast_nor)
    uniqList,pairedSet,allSet,unassignList  = inferIRs(list_rc,header,phasSet)
    isoDict                                 = isoformsDict(list_nor,phasSet,allSet,unassignList)
    isoformChecker(phasSet,unassignList,pairedSet,isoDict,list_rc,header,phasList)

if __name__ == '__main__':
    main()
    sys.exit()



## v02 -> v03
## Fixed isoform catching, earlier same unassigned transcript would have been added multiple times to noIsoform set exaggarating its number
## because if multiple isforms same transcript is added at every loop