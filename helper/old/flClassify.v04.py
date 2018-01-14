#!/usr/local/bin/python3


### This script accepts 1) PHAS file from (PHASredundant script), 
#### 2) a GTF file and 
#### 3) target prediction file from sPARTA to identify transcripts that are FL or non-FL

import os,sys,time,sqlite3,operator

######### USER INPUTS ###############
featureFile     = "Trinity.collapsed.gtf"
gtfmode         = 3                                                         ## 1: PASA feature (GTF) file | 2: Trinity | 3: Rocket feature (GTF) file or generic gtf file with 
                                                                            ## missing "transcript" entry | 4: Both rocket and trinity feature files | 5: PacBio GTF | 6: Comprehensive 
                                                                            ## transciptome (from PASA) - NOTE: New mode should be registered in overlapchecker function

phasFile        = "Final_24PHAS_Loci_5e-07_ALL.csv"
delim           = "\t"                                                      ## Delimiter for summary file
head            = 'Y'                                                       ## Header is summary file: 'Y' else: 'N'
phasName        = 1
chrpos          = 3
startpos        = 4
endpos          = 5


tarFile         = "All.2118-2275.genic.inter.targs.parsed.csv_revmapped.csv"
PAREresType     = 'T'
paremode        = 0
PAREpval        = 0.25
tarScore        = 7
revmapped       = 1                                         ## Fixes 13th position observed in 'C' strand

prmLength       = 1000                                      ### Length of upstream promoter required

######## Functions ###################
######################################

def findFL(featureFile,phasFile):
    '''Find which of the phased trancript has miR2118.miR2275 site internal to transcript i.e. it's not a cleaved trancript'''

    ### parse the GTF, PHAS 
    gtfList     = gtfParser2(featureFile,gtfmode)
    phasList    = flPhasParser(phasFile)
    ### Parse targets file 
    tarList,header = PAREreader(tarFile,PAREpval,paremode)

    ##########################################################
    DB = 'tempdb'
    try:
        os.remove(DB)
    except OSError:
        pass
    ## Prepare features table - This might include entries from two GTF files (mode-2), identified by the flags
    conn            = sqlite3.connect(DB)
    featureTable    = tableMaker(gtfList,conn) ## TURN ON ****

    ### Work on each PHAS and decide whther it's FL or not #########
    ################################################################
    fh_out = open("fl.list","w") ### Just to record names of FL ones
    flDict = {} ### Dict of transcripts that are assigned FL status
    acount = 0  ### Total PHAS
    bcount = 0  ### Total FL
    for phas in phasList:
        acount += 1
        print("\n\nAnalyzing PHAS:",phas) ### aname,achr,astart,aend,astrand
        pname,pchr,pstart,pend,pstrand = phas

        ### Get GTF entries for PHAS
        phastrans,phasexons = getFeatures(phas,featureTable,conn) ## fetches transcipt and exons entry from gtf file

        ### Check for targets in the transcript
        gname,tname,achr,astrand,astart,aend,finalstatus = targetChecker(tarList,phastrans,phasexons)
        print("FinalStatus:%s" % (finalstatus))
        if finalstatus == "fl":
            bcount += 1
            fh_out.write("%s\t%s\n" % (gname,tname))
            flDict[pchr] = (achr,astrand,astart,aend) ### Becuase chr columsn in PHAS lict is the transcript name
        else:
            pass

    ### Get promoter coords
    print("Total PHAS:%s | FL PHAS:%s" % (acount,len(flDict)))
    longPrm,uniqPrm     = uniqpromoter(flDict)
    longPrmFa,uniqPrmFa = writer(longPrm,uniqPrm)

    fh_out.close()

    return flDict

def writer(longPrm,uniqPrm):
    '''
    writes coord file for uniq promoter i.e. one per loci and merged promoter 
    '''

    ### Write long promoter
    outfile1 = "%s.longprm.fa" % (phasFile.rpartition(".")[0])
    fh_out1  = open(outfile1,"w")
    acount  = 0 ### Number of entries written
    for i in longPrm:
        acount+=1
        agene,achr,astrand,prm1start,prm1end = i
        fh_out1.write("%s\t%s\t%s\t%s\t%s\n" % (agene,achr,astrand,prm1start,prm1end))

    ### Write unique promoter
    outfile2 = "%s.uniqprm.fa" % (phasFile.rpartition(".")[0])
    fh_out2  = open(outfile2,"w")
    for i in uniqPrm:
        agene,achr,astrand,prm2start,prm2end = i
        fh_out2.write("%s\t%s\t%s\t%s\t%s\n" % (agene,achr,astrand,prm2start,prm2end))

    print("Entries written:%s" % (acount))

    fh_out1.close()
    fh_out2.close()

    return outfile1,outfile2

def uniqpromoter(flDict):
    '''
    script takes a list of fl transcripts, and computes promoter regions by including upstream length and distance from rightest
    transcript till leftest transcript i.e. if three transcipts for same gene exist then one with 5' start at extreme right will be used to first compute 
    region till 5' start site of the leftest transcript amd then addition upsteeam promoter region will be extracted
    '''
    analyzedList = [] ### List of genes already processed

    ### Gene by gene find the promoter region
    longPrm = [] ## Promoter coords for all uniq loci, promoter includes gene overlapping region from rightest 5' end 
    uniqPrm = [] ## Promoter coords for all uniq loci, promoter is computed from the leftest 5' start i.e. no overlapping coding region  
    for atrans in flDict.keys():
        agene                       = atrans.rpartition(".")[0]
        achr,astrand,astart,aend    = flDict[atrans]

        #### If it's a new loci then analyze ###
        ########################################
        if agene not in analyzedList:
            print("\nComputing promoter for gene:%s | trans:%s" % (agene,atrans))
            analyzedList.append(agene)
            start5List = [] ### List to store 5' ends of all transcripts from this gene

            ## Find all the transcipts for this gene and record their 5' ends 
            ################################################################
            bcount = 0 ### Transcripts searched for same loci
            for i in flDict.keys():
                bcount+=1
                bgene = i.rpartition(".")[0]
                bchr,bstrand,bstart,bend = flDict[i]
                if agene == bgene:
                    print("Loci found:",i,bchr,bstrand,bstart,bend)

                    if astrand == bstrand and bstrand == "w": ## Check to make sure that both are on same strand and if "w"
                        start5 = int(bstart)
                        start5List.append(start5)

                    elif astrand==bstrand and bstrand == "c":
                        start5 = int(bend)
                        start5List.append(start5)

                    else:
                        print("That means current gene had transcripts from different strands - need a new strategy now")
                        pass
                else:
                    # print("%s doesn't belong to loci under scan:%s" % (i,agene))
                    pass


            ### Find extreme 5' ends from all transcripts of this gene  #####
            print("%s entries searched for transcipts to genes:%s" % (bcount,agene))
            print("Start Sites for phased transcipts in this gene:",start5List)
            min5    = min(start5List)
            max5    = max(start5List)

            if astrand      == "w":
                prm1start   = min5-prmLength
                prm1end     = max5
                prm2start   = min5-prmLength
                prm2end     = min5
                longPrm.append((agene,achr,"w",prm1start,prm1end))
                uniqPrm.append((agene,achr,"w",prm2start,prm2end))
            elif astrand    == "c":
                prm1start   = min5
                prm1end     = max5+prmLength
                prm2start   = max5
                prm2end     = max5+prmLength
                longPrm.append((agene,achr,"c",prm1start,prm1end))
                uniqPrm.append((agene,achr,"c",prm2start,prm2end))
            else:
                print("Unknown strand encountered:%s" % (astrand))
                print("Debug - script will exit")
                sys.exit()

        else:
            print("This loci already processed")
            pass

    return longPrm,uniqPrm

def targetChecker(tarList,phastrans,phasexons):
    '''Takes a feature list and returns back if target sites overlapping features'''
    
    for atrans in phastrans: ## For each transcript
        statusList = [] ## This holds status list for all miRNA target sites - there could be multiple target site sin one trans

        gname,tname,xchr,xstart,xend,astrand,atype,aflag = atrans
        achr    = int(xchr)
        astart  = int(xstart)
        aend    = int(xend)
        print("-Feature being scanned for targets: Gene:%s | Trans:%s | Chr:%s | Start:%s | End:%s | Strand:%s | Type:%s | Flag:%s" % (gname,tname,achr,astart,aend,astrand,atype,aflag))
        
        # ### Strand aware transcript coords ####
        # if astrand == 'c':
        #     astart  = xend
        #     aend    = xstart
        # elif astrand == 'w':
        #     astart  = xstart
        #     aend    = xend
        # else:
        #     print("unknown strand encountered:%s" % (astrand))


        
        for i in tarList: ## For every target in list
            # print("Target Site:",i) #### (mirName,tarName,cleaveSite,chrid,strand,i.strip('\n')
            mirname     = i[0].replace(">","")
            chrid       = int(i[3])
            strand      = i[4]
            ent_splt    = i[5].split(',')
            bindstart   = int(ent_splt[10])
            bindend     = int(ent_splt[11])
            # print("mir:%s | chr:%s | strand: %s | bindstart:%s | bindend:%s" % (mirname,chrid,strand,bindstart,bindend))

            # #### Strand aware bind sites ########################
            # ########################################
            # ent_splt    = i[5].split(',')
            # if strand == 'w':
            #     bindend     = int(ent_splt[10])
            #     bindstart   = int(ent_splt[11])
            # elif strand == 'c':
            #     bindstart   = int(ent_splt[10])
            #     bindend     = int(ent_splt[11])
            # else:
            #     print("please check the strand for target:%s" % (strand))
            #     pass

            ### Check for overlap with transcript #####
            ###########################################
            if achr == chrid and astrand == strand and bindstart > astart and bindend < aend:
                print("\n-Target site on transcript - mir:%s | chr:%s | strand: %s | bindstart:%s | bindend:%s" % (mirname,chrid,strand,bindstart,bindend))
                # sys.exit()

                if astrand      == 'w':
                    reminder5   = bindstart - astart
                    print("-Trans reminder:%s" % (reminder5))
                elif astrand    == 'c':
                    reminder5   = aend - bindend
                    print("- Trans. reminder:%s" % (reminder5))
                else:
                    print("-Stupid fuck - The strand is unknown:%s" % (astrand))

                ### Check if target site is on exon ###
                #######################################
                for aexon in phasexons:
                    gname,tname,zchr,zstart,zend,bstrand,btype,bflag = aexon
                    bchr    = int(zchr)
                    bstart  = int(zstart)
                    bend    = int(zend)
                    print("-Exon being scanned for targets: Gene:%s | Trans:%s | Chr:%s | Start:%s | End:%s | Strand:%s | Type:%s | Flag:%s" % (gname,tname,bchr,bstart,bend,bstrand,atype,aflag))

                    if bchr == chrid and bstrand == strand and bindstart > bstart and bindend < bend:
                        print("-Target site on exon")

                        if astrand      == 'w':
                            reminder5   = bindstart - astart
                            astatus     = "FL-Ex" 
                            print("-Exon reminder:%s" % (reminder5))
                            statusList.append(astatus)
                        elif astrand    == 'c':
                            reminder5   = aend - bindend
                            astatus     = "FL-Ex"
                            print("-Exon reminder:%s" % (reminder5))
                            statusList.append(astatus)
                        else:
                            print("Stupid fuck - The strand is unknown:%s" % (astrand))

                    else:
                        print("-Target site not on exon")
                        astatus     = "FL-In"
                        statusList.append(astatus)
                        pass

            else:
                # print("target site not on transcript")
                pass


        ### Count how many miRNAs have target site that suggest trans is FL
        flcount     = statusList.count("FL-Ex")
        flcount2    = statusList.count("FL-In")
        if flcount > 0:
            finalstatus = "fl"
        else:
            finalstatus = "no"

    return gname,tname,achr,astrand,astart,aend,finalstatus

def getFeatures(phas,featuretable,conn):
    '''
    Reports back the transcripts and exons for PHAS transcript
    '''

    transList = [] ## Store the transcipts that overlap

    aname,achr,astart,aend,astrand = phas
    cur = conn.cursor()

    # print("-FLag being queried:%s" % (aflag))
    print("-Query entry:",phas[1])

    ### Test 
    cur.execute("SELECT * FROM %s limit 10" % (featuretable))
    testquery = cur.fetchall()
    # print("Results from test query:",testquery)
    
    ## trans flanking phas or enclaved in phas
    cur.execute("SELECT * FROM %s where trans = '%s' AND type = 'transcript'" % (featuretable,str(achr)))
    phastrans = cur.fetchall()

    cur.execute("SELECT * FROM %s where trans = '%s' AND type = 'exon'" % (featuretable,str(achr)))
    phasexons = cur.fetchall()

    print("-Number of PHAS trans:%s" % (len(phastrans)))
    # print("-transList",phastrans)

    print("-Number of PHAS Exons:%s" % (len(phasexons)))
    # print("-transList",phasexons)

    return phastrans,phasexons

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

def gtfParser2(afeatureFile,gtfmode):

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
            gchr    = ent[0].replace("chromosome","").replace("chr","")
            gstart  = int(ent[3])
            gend    = int(ent[4])
            gstrand = ent[6].translate(str.maketrans("+-","wc"))
            info    = ent[8].strip("\n").split(";")
            # print(info)
            
            ## Parse the info and add the exon entries  #################
            #############################################################

            if gtfmode == 2 or (gtfmode == 4 and len(info)==4): ## This is  trinity GTF info
                ## Protein coding gene with a version number
                tid     = info[0].split()[1].replace('"','').split(".")[0] ## Transcript ID
                gid     = info[2].split()[1].replace('"','').rpartition("_")[0] ## Gene ID
                aflag   = 'T' ## Trinity
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))

            elif gtfmode == 3 or (gtfmode ==4 and len(info) >= 7): ## This is rocket GTF info, with and w/o p_id "P1" in the end
                gid     = info[1].split()[1].replace('"','') ## Gene ID
                tid     = info[0].split()[1].replace('"','') ## Transcript ID
                aflag   = 'R' ## Rocket
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))
            
            elif gtfmode == 5:
                tid     = info[0].split()[1].replace('"','').split(".")[0] ## Transcript ID
                gid     = info[1].split()[1].replace('"','').split(".")[0] ## Gene ID
                aflag   = 'PB' ## PacBio
                # print('-',gid,tid,gchr,gstart,gend,gstrand,gtype,aflag)
                gtfList.append((gid,tid,gchr,gstart,gend,gstrand,gtype,aflag))

            elif gtfmode ==6:
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
        print("First 10 entries of gtfList list: %s" % (gtfList[:20]))

    print("Total entries fetched from GTF file:%s" % (str(len(gtfList))))

    print ("Exiting function - gtfParser2\n")
    # sys.exit()
    time.sleep(1)
    
    return gtfList

def flPhasParser(phasFile):
    '''parses input PHAS file for collasped transcripts with lot of BS'''

    print("\nFunction: phasParser")

    fh_in = open(phasFile,'r')
    if head == 'Y':
        fh_in.readline()
    phasRead = fh_in.readlines()

    phasList = []
    for i in phasRead:
        ent = i.split(delim)
        aname   = ent[phasName-1]
        achr    = ent[chrpos-1].split("|")[0]
        astart  = int(ent[startpos-1])
        aend    = int(ent[endpos-1])
        # astrand = ent[strandpos-1]
        astrand = "na"
        print("Phas Name %s | chr:%s | start:%s | end:%s | astrand:%s" % (aname,achr,astart,aend,astrand))
        phasList.append(((aname,achr,astart,aend,astrand)))

    print("A list of query features prepared with entries:%s" % len(phasList))

    print ("Exiting function - phasParser\n")
    return phasList

def PAREreader(PAREres,PAREpval,mode):
    
    fh_in = open(PAREres, 'r')
    header = fh_in.readline() ## Header
    entries = fh_in.readlines()
    acount = 0  ## Total entries
    bcount = 0  ## Filtered entries on targetScore (predicted) or p-value (validated)
    resList = [] ## miRNA,Target,cleavesite,whole entry
    # print ("\n\n\nCreating list of cleave sites\n")

    
    if PAREresType == 'S' and paremode == 0: ## sPARTA server
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.split(',')
            # print("\nThis is target entry p-val:",ent_splt[13])
            if float(ent_splt[13]) <= PAREpval:
                cleaveSite = int(ent_splt[18])
                #print(cleaveSite)
                mirName = ent_splt[0]
                tarName = ent_splt[1]
                chrid   = ent_splt[15]
                strand  = ent_splt[16]
                resList.append((mirName,tarName,cleaveSite,chrid,strand,i.strip('\n')))
                # print("\nList values for this entry",mirName,tarName,cleaveSite,chrid,strand,i.strip('\n'))
                bcount += 1

    elif PAREresType == 'S' and paremode == 1: ## sPARTA Local
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.split(',')
            #print("\nThis is target entry",ent_splt[13])
            if float(ent_splt[14]) <= PAREpval:
                cleaveSite = int(ent_splt[8])
                #print(cleaveSite)
                tarName = ent_splt[1]
                mirName = ent_splt[0]
                chrid   = ent_splt[1] ## In local mode there is no chr_id, only transcripts
                strand  = 'None'
                resList.append((mirName,tarName,cleaveSite,chrid,strand,i.strip('\n')))
                # print(mirName,tarName,cleaveSite)
                bcount += 0

    elif PAREresType == 'T' and paremode == 0: ## revFerno targets i.e. miRferno revmapped using revFerno v02 or above
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.strip('\n').split(',')
            chrid       = ent_splt[8] ## In local mode there is no chr_id, only transcripts
            strand      = ent_splt[9]
            # print("\nThis is target entry",ent_splt)
            if float(ent_splt[5]) <= tarScore:
                # print("pass")
                bindStart   = ent_splt[10]
                bindEnd     = ent_splt[11]
                if revmapped == 1 and strand == "c":
                    ## This is to accomodate 13th psoition validation observed for 'c' strand
                    cleaveSite3 = int(bindEnd)-12+1 ## 12th pos +1  - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct ## Added after observing that 2275 cleaves a major portion at 13th pos
                    cleaveSite4 = int(bindEnd)-13+1 ## 13th pos +1 - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct
                    cleaveSite  = (cleaveSite3,cleaveSite4)
                else:
                    ## Watson strand
                    cleaveSite1 = int(bindEnd)-10+1 ## 10th pos +1  - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct
                    cleaveSite2 = int(bindEnd)-11+1 ## 11th pos +1 - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct
                    cleaveSite  = (cleaveSite1,cleaveSite2)
                # print(cleaveSite)
                mirName     = ent_splt[0]
                tarName     = ent_splt[1]
                resList.append((mirName,tarName,cleaveSite,chrid,strand,i.strip('\n')))
                # print("List values for this entry",mirName,tarName,cleaveSite,chrid,i.strip('\n'))
                bcount += 1

    elif PAREresType == 'T' and paremode == 1: ## sPARTA local targets only, no PARE
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.strip('\n').split(',')
            # print("\nThis is target entry",ent_splt)
            if float(ent_splt[5]) <= tarScore:
                bindStart,bindEnd = ent_splt[2].split("-")
                cleaveSite1 = int(bindEnd)-10+1 ## 10th pos +1  - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct
                cleaveSite2 = int(bindEnd)-11+1 ## 11th pos +1 - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct
                # cleaveSite3 = int(bindEnd)-12+1 ## 12th pos +1  - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct ## Added after observing that 2275 cleaves a major portion at 13th pos
                # cleaveSite4 = int(bindEnd)-13+1 ## 13th pos +1 - tested and corrected - Draw a miRNA-target interaction and try computing - This is correct
                cleaveSite  = (cleaveSite1,cleaveSite2,cleaveSite3,cleaveSite4)
                print(cleaveSite)
                mirName     = ent_splt[0]
                tarName     = ent_splt[1]
                chrid       = ent_splt[1] ## In local mode there is no chr_id, only transcripts
                strand      = 'None'
                resList.append((mirName,tarName,cleaveSite,chrid,strand,i.strip('\n')))
                # print("List values for this entry",mirName,tarName,cleaveSite,chrid,i.strip('\n'))
                bcount += 1
        
    ### Deprecated #### NON functional
    elif PAREresType == 'C': ## CL3
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.split(',')
            # print("\nThis is target entry",ent_splt)
            if float(ent_splt[9]) <= PAREpval:
                mirName     = ent_splt[0]
                tarName     = ent_splt[1]
                cleaveSite  = int(ent_splt[5])
                print (cleaveSite)
                resList.append((mirName,tarName,cleaveSite,i.strip('\n')))
                bcount += 1

    print("Input file:%s | File entries:%s | Passed p-val/score:%s | List length:%s" % (PAREres,acount,bcount,len(resList)))
    
    return resList,header

def main():
    findFL(featureFile,phasFile)
    pass



if __name__ == "__main__":
    main()
    print("Script finished sucessfully")
    sys.exit()