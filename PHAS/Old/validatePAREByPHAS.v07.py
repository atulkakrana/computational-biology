#!/usr/local/bin/python3
## v02 has modification for new MPPP format i.e with category column
## v07 is modified to use with sPARTA v1.16 or above and revFerno v02 or above

import sys,time,os,glob
import mysql.connector as sql

#### Settings
genomeDB    = 'ASPARAGUS_UGA2_genome' ## Used to identify strand of target
dataserver  = 'raichu.dbi.udel.edu'

mode        = 0                                         ## 0: Both target and PHAS coordinates match genome DB 
                                                        ## 1: PHAS prediction and target analysis done on transcripts or non-genome-DB targets (local analysis) like in case of PacBio transcripts
                                                        ## Files for these two modes can't be fixed

target      = 'C'                                       ## 'C' if genomic coordinates (mode 0) 'T' is transcripts used for target validation/prediction (mode 1)  
                                                        ## i.e. what is in the chromosome column of phasing results - transcript name or chromosome

PAREres     = 'genic_inter_All.targs.parsed.aspa.v2.csv_revmapped.csv'             ## Mode = 0 revmapped file with 14th column as cleave site | Mode: 1 -local validation file with cleave site at column 9
PAREresType = 'T'                                       ## 'S' for sPARTA PARE validated, 'C' for CL3 PARE Validated and 'T' for just sPARTA target prediction file
tarScore    = 5                                         ## Only if PAREresType == 'T' - Score cutoff
PAREpval    = 0.25                                      ## Only if PAREresType == 'S' - Cutoff of corrected p-value to use to filter out non-relevant predictions

phase       = 21                                        ## Phase length that you are using for double validation - Give corresponding file below
PHASresType = 'PP'                                      ## J: jixian results format | P:Ping list file | PP: Phas redundant script generated *.ping format
PHASres     = 'Final_21PHAS_Loci_ALL_v6.ping'
pVal        = '0.005'                                   ## P-value of phased loci to test for

offset      = 1                                         ## 0: No +1/-1 phases scanned 1: In addition to phase +1/-1 sites also checked for trigger
 

##Read phased loci file and make dictionary
def CreatePHASdict(PHASres,phase,pVal):
    ''' Reads the phased result file
    and creates a dictionary of -5/+7 phased locations as value
    '''
    
    fh_in = open(PHASres,'r')
    #fh_out = open('%s_Dict' % (PHASres),'w')
    PHASdict = {}
    entries = fh_in.readlines()
    
    print ("Creating dictionary of phased loci\n\n\n")
    if PAREresType == 'S' or PAREresType == 'C' or PAREresType == 'T':
        if phase == 21:
            phaseList = [-105,-84,-63,-42,-21,0,21,42,63,84,105,126,147]
        elif phase == 24:
            phaseList = [-120,-96,-72,-48,-24,0,24,48,72,96,120,144,168]
        else:
            print('What the phase for phased loci - Please specify correctly')
            sys.exit()
    
    if PHASresType == 'P':
        for i in entries:
            if i.strip(): ## Remove an empty line from end file
                ent_splt = i.strip('\n').split('=')
                #print(ent_splt[0].split('|'))
                pval,phase,trash = ent_splt[0].split('|')
                chromo_start,end = ent_splt[1].split('..')
                chromo,start = chromo_start.split(':')
                
                if target == 'T':
                    akey ='%s-%s-%s' % (chromo.strip().split("|")[0],start,end) ## Since its local analysis, transcrpt name would be 
                else:
                    akey ='%s-%s-%s' % (chromo.strip().replace('chr',''),start,end)
                
                aval = [sum(i) for i in zip([int(start)]*11,phaseList)]
                print ('Key:%s | Value:%s' % (akey,aval))
                if pval == pVal:
                    PHASdict[akey] = aval
                    #print ('Key:%s | Value:%s' % (akey,aval))
                    #fh_out.write('%s\t%s\t%s\t%s\t%s\tNONE\tNONE\n' % (phase,pval,chromo.strip(),start,end))##Chromosome has space before it which later gives error while key matching
                    #fh_out.write('%s\t%s\n' % (str(akey),str(aval).strip('[]')))
    
    elif PHASresType == 'PP': ##For processes/merged scenario - use *.ping file
        for i in entries:
            if i.strip(): ## Remove an empty line from end file
                phase,pval,chromo,start,end,trash1,trash2 = i.strip('\n').split('\t') ## In case of local analysis chromo is full transcript name
                print("\nPhased entry:",phase,pval,chromo,start,end,trash1,trash2)
                if target == 'T':
                    akey ='%s-%s-%s' % (chromo.strip().split("|")[0],start,end) ## Since its local analysis, transcrpt name would be 
                else:
                    akey ='%s-%s-%s' % (chromo.strip().replace('chr',''),start,end)
                
                aval = [sum(i) for i in zip([int(start)]*11,phaseList)]
                print("This is the PHAS key: %s | value: %s" % (akey,aval))
                
                if pval == pVal:
                    print("PHAS entry being added to list for validation\n")
                    PHASdict[akey] = aval
    

    elif PHASresType == 'J':
        for i in entries:
            if i.strip(): ## Remove an empty line from end file
                ent_splt = i.strip('\n').split('\t')
                chromo = ent_splt[1].split('_')[1]
                start,end = ent_splt[4].split('-')
                if target == 'T':
                    akey ='%s-%s-%s' % (chromo.strip().split("|")[0],start,end) ## Since its local analysis, transcrpt name would be 
                else:
                    akey ='%s-%s-%s' % (chromo.strip().replace('chr',''),start,end)
                aval = [sum(i) for i in zip([int(start)]*11,phaseList)]
                PHASdict[akey] = aval
                print ('Key:%s | Value:%s' % (akey,aval))
                #fh_out.write('%s\t%s\t%s\t%s\t%s\tNONE\tNONE\n' % (phase,pval,chromo.strip(),start,end))##Chromosome has space before it which later gives error while key matching
                #fh_out.write('%s\t%s\n' % (str(akey),str(aval).strip('[]')))
    else:
        print('Please enter correct PHASED result type')
        pass
    fh_in.close()
    #fh_out.close()

        
    return PHASdict

##Reads CL and MPPP file to return alist
def PAREreader(PAREres,PAREpval):
    
    fh_in = open(PAREres, 'r')
    header = fh_in.readline() ## Header
    entries = fh_in.readlines()
    acount = 0  ## Total entries
    bcount = 0  ## Filtered entries on targetScore (predicted) or p-value (validated)
    resList = [] ## miRNA,Target,cleavesite,whole entry
    # print ("\n\n\nCreating list of cleave sites\n")

    
    if PAREresType == 'S' and mode == 0: ## sPARTA server
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
                chrid = ent_splt[15]
                strand = ent_splt[16]
                resList.append((mirName,tarName,cleaveSite,chrid,i.strip('\n')))
                # print("\nList values for this entry",mirName,tarName,cleaveSite,chrid,i.strip('\n'))
                bcount += 1

    elif PAREresType == 'S' and mode == 1: ## sPARTA Local
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
                chrid = ent_splt[1] ## In local mode there is no chr_id, only transcripts
                strand = 'None'
                resList.append((mirName,tarName,cleaveSite,chrid,i.strip('\n')))
                # print(mirName,tarName,cleaveSite)
                bcount += 0

    elif PAREresType == 'T' and mode == 0: ## revFerno targets i.e. miRferno revmapped using revFerno v02 or above
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.strip('\n').split(',')
            # print("\nThis is target entry",ent_splt)
            if float(ent_splt[5]) <= tarScore:
                # print("pass")
                bindStart = ent_splt[10]
                bindEnd = ent_splt[11]
                cleaveSite1 = int(bindEnd)-10 ## 10th pos
                cleavSite2 = int(bindEnd)-11 ## 11th pos
                cleaveSite = (cleaveSite1,cleavSite2)
                # print(cleaveSite)
                mirName = ent_splt[0]
                tarName = ent_splt[1]
                chrid = ent_splt[8] ## In local mode there is no chr_id, only transcripts
                strand = ent_splt[9]
                resList.append((mirName,tarName,cleaveSite,chrid,i.strip('\n')))
                # print("List values for this entry",mirName,tarName,cleaveSite,chrid,i.strip('\n'))
                bcount += 1

    
    elif PAREresType == 'T' and mode == 1: ## sPARTA local targets only, no PARE
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.strip('\n').split(',')
            # print("\nThis is target entry",ent_splt)
            if float(ent_splt[5]) <= tarScore:
                bindStart,bindEnd = ent_splt[2].split("-")
                cleaveSite1 = int(bindEnd)-10 ## 10th pos
                cleavSite2 = int(bindEnd)-11 ## 11th pos
                cleaveSite = (cleaveSite1,cleavSite2)
                print(cleaveSite)
                mirName = ent_splt[0]
                tarName = ent_splt[1]
                chrid = ent_splt[1] ## In local mode there is no chr_id, only transcripts
                strand = 'None'
                resList.append((mirName,tarName,cleaveSite,chrid,i.strip('\n')))
                # print("List values for this entry",mirName,tarName,cleaveSite,chrid,i.strip('\n'))
                bcount += 1
        

    elif PAREresType == 'C': ## CL3
        for i in entries:
            #print (i)
            acount += 1  ## Total entries
            ent_splt = i.split(',')
            # print("\nThis is target entry",ent_splt)
            if float(ent_splt[9]) <= PAREpval:
                mirName = ent_splt[0]
                tarName = ent_splt[1]
                cleaveSite = int(ent_splt[5])
                print (cleaveSite)
                resList.append((mirName,tarName,cleaveSite,i.strip('\n')))
                bcount += 1
        
    
    print("Input file:%s | File entries:%s | Passed p-val/score:%s | List length:%s" % (PAREres,acount,bcount,len(resList)))
    
    return resList,header

##Identify if CleaveSite in phasedloci
def validatePHAS(resList,PHASdict,phase,con,header):
    '''
    '''
    fh_out = open('%s_PHASmatched.csv' % (PAREres),'w')
    fh_out.write('%s,PHASLoci,PHASindex\n' % (header.strip('\n')))
    totalCount = 0 ##Counter for total entries in reslist
    phasedCount =0 ##Counter for matched entries from the reslist
     
    for ent in resList: ##ent format: miRNA,Target,cleavesite,whole entry
        totalCount +=1
        #print (ent)
        mirName = ent[0]
        tarName = ent[1]
        chrid = ent[3]
        if PAREresType == 'T':
            cleaveSite1,cleaveSite2 = ent[2] ## Predicted 10th and 11th as cleave site
            print("Scanning phased loci for these cleave sites:",mirName,tarName,cleaveSite1,cleaveSite2,chrid)
        else:
            cleaveSite = int(ent[2])    ## PARE supported cleave site
            print("Scanning phased loci for this cleave site:",mirName,tarName,cleaveSite,chrid)
        
        if mode == 0: ## targets coordinates match to genome DB i.e. server target analysis was done
            
            cur= con.cursor()
            #print('Query DB:',genomeDB,'| Query Target:',tarName.split('_up')[0].split('_down')[0])
            
            ## Get strand of the predicted/validated entry
            ## In v07 - strand is included in resList - this comes from sPARTA or revferno - so this 'if' loop is really not required - i.e. no my SQL functionality required
            if tarName.rpartition("_")[2] == "all":
                ## sPARTA1.16 onwards scaffolds with no genes were also included in analysis like chr10715_c_all
                tarStrand = tarName.rsplit("_")[1]
                print("Target strand for entry with no annotated gene:%s" % (tarStrand))
            
            else:
                cur.execute("SELECT strand FROM %s.gene_master where gene like '%s'" % (genomeDB,tarName.split('_up')[0].split('_down')[0]))### Convert intergenic to gene name so as to get strand
                info = cur.fetchall()
                tarStrand = info[0][0]
                #print (tarStrand)
            
            ## Offset cleave site if on crick strand
            if tarStrand == 'c': ##Add 3 to cleave site
                #print("\nCrick strand:%s" % (tarStrand))
                if PAREresType == 'S':
                    cleaveSite += 3    
                elif PAREresType == 'T':
                    cleaveSite1 += 3
                    cleaveSite2 += 3
                else:
                    print("Check PAREresType")
                    sys.exit()

            else: ## Do nothing to cleave site
                #print("\nWatson strand:%s" % (tarStrand))
                pass
            
            ## DO validation
            if PAREresType == 'S' or PAREresType == 'C':
                ##PHAS LOCI MATCH - Fuck you -  Was not checking for cleave site in exact chromosome - Thats why you got shit load of results where these are scaffolds
                for akey in PHASdict.keys():
                    aval = PHASdict[akey]
                    # print ('This is the key:', akey, aval)
                    if cleaveSite in aval:
                            if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite,tarStrand,akey,aval,(aval.index(cleaveSite)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval.index(cleaveSite)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                phasedCount += 1
                            else:
                                # print ('Cleave site matched phase index but chromosome was different')
                                pass           

                    else:
                        ## Cleave site doesnt matches with phase site, try offsets maybe
                        if offset == 1:
                            aval1 = [x+1 for x in aval]
                            aval2 = [x-1 for x in aval]
                            if cleaveSite in aval1:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite,tarStrand,akey,aval,(aval1.index(cleaveSite)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval1.index(cleaveSite)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1

                            elif cleaveSite in aval2:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite,tarStrand,akey,aval,(aval2.index(cleaveSite)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval2.index(cleaveSite)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1     
                            
                            else:
                                # print ('Cleave site matched phase index but chromosme was different')  
                                pass

                        else:
                            print("Offset checking is OFF")
                            pass

            elif PAREresType == 'T': ## miRFerno results with two possible cleave sites both needs to be checked
                
                for akey in PHASdict.keys():
                    aval = PHASdict[akey]
                    # print ('This is the key:', akey, aval)
                    if cleaveSite1 in aval:
                            if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite1,tarStrand,akey,aval,(aval.index(cleaveSite1)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval.index(cleaveSite1)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                phasedCount += 1
                            else:
                                # print ('Cleave site matched phase index but chromosome was different')
                                pass

                    elif cleaveSite2 in aval:
                            if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite2,tarStrand,akey,aval,(aval.index(cleaveSite2)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval.index(cleaveSite2)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                phasedCount += 1
                            else:
                                # print ('Cleave site matched phase index but chromosome was different')
                                pass
                    
                    else:
                        ## Cleave site doesnt matches with phase site, try offsets maybe
                        if offset == 1:
                            aval1 = [x+1 for x in aval]
                            aval2 = [x-1 for x in aval]
                            if cleaveSite1 in aval1:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite1,tarStrand,akey,aval,(aval1.index(cleaveSite1)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval1.index(cleaveSite1)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1

                            elif cleaveSite2 in aval1:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite2,tarStrand,akey,aval,(aval1.index(cleaveSite2)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval1.index(cleaveSite2)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1 

                            elif cleaveSite1 in aval2:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite1,tarStrand,akey,aval,(aval2.index(cleaveSite1)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval2.index(cleaveSite1)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1

                            elif cleaveSite2 in aval2:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite2,tarStrand,akey,aval,(aval2.index(cleaveSite2)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval2.index(cleaveSite2)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1    
                            
                            else:
                                # print ('Cleave site matched phase index but chromosme was different')  
                                pass


            else:
                print("Input correct PARE result type in user settings")
                sys.exit()


        ##### LOCAL MODE ###########################################
        ############################################################
        
        elif mode == 1: ## Local analysis - Targets are transcripts and no chr_id or strand is required from server

            #### PARE validated cleave sites #######################
            #######################################################
            if PAREresType == 'S' or PAREresType == 'C':
                for akey in PHASdict.keys():
                    aval = PHASdict[akey]
                    #print ('This is the key:', akey, aval)
                    if cleaveSite in aval:
                            if akey.split('-')[0] == chrid: ## i.e. clevage site is on same transcript
                                print('\nmiRNA:%s | Target:%s |CleaveSite:%s | PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite,akey,aval,(aval.index(cleaveSite)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval.index(cleaveSite)+1-6)))
                                phasedCount += 1
                            else:
                                print ('Cleave site matched phase index but chromosme was different')            
                    else:
                        ## Cleave site doesnt matches with phase site, try offsets maybe
                        if offset == 1:
                            aval1 = [x+1 for x in aval]
                            aval2 = [x-1 for x in aval]
                            if cleaveSite in aval1:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite,tarStrand,akey,aval,(aval1.index(cleaveSite)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval1.index(cleaveSite)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1


                            elif cleaveSite in aval2:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite,tarStrand,akey,aval,(aval2.index(cleaveSite)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval2.index(cleaveSite)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1


                            else:
                                # print ('Cleave site matched phase index but chromosme was different')  
                                pass

            #### Cleave sites from miRferno #####################
            #####################################################
            elif PAREresType == 'T':
                for akey in PHASdict.keys():
                    aval = PHASdict[akey]
                    print ('This is the key:', akey, aval)
                    # print(cleaveSite1,cleaveSite2)
                    if cleaveSite1 in aval:
                            if akey.split('-')[0] == chrid: ## i.e. clevage site is on same transcript
                                print('\nmiRNA:%s | Target:%s |CleaveSite:%s | PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite1,akey,aval,(aval.index(cleaveSite1)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval.index(cleaveSite1)+1-6)))
                                phasedCount += 1
                            else:
                                print('\nmiRNA:%s | Target:%s |CleaveSite1:%s | PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite1,akey,aval,(aval.index(cleaveSite1)+1-6)))
                                print ('Cleave site matched phase index but chromosme was different') 

                    elif cleaveSite2 in aval:
                            if akey.split('-')[0] == chrid: ## i.e. clevage site is on same transcript
                                print('\nmiRNA:%s | Target:%s |CleaveSite:%s | PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite2,akey,aval,(aval.index(cleaveSite2)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval.index(cleaveSite2)+1-6)))
                                phasedCount += 1
                            else:
                                print('\nmiRNA:%s | Target:%s |CleaveSite2:%s | PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite2,akey,aval,(aval.index(cleaveSite2)+1-6)))
                                print ('Cleave site matched phase index but chromosme was different')            
                    else:
                        ## Cleave site doesnt matches with phase site, try offsets maybe
                        if offset == 1:
                            aval1 = [x+1 for x in aval]
                            aval2 = [x-1 for x in aval]
                            if cleaveSite1 in aval1:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite1,tarStrand,akey,aval,(aval1.index(cleaveSite1)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval1.index(cleaveSite1)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1
                            
                            elif cleaveSite2 in aval1:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite2,tarStrand,akey,aval,(aval1.index(cleaveSite2)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval1.index(cleaveSite2)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1

                            elif cleaveSite1 in aval2:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite1,tarStrand,akey,aval,(aval2.index(cleaveSite1)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval2.index(cleaveSite1)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1

                            elif cleaveSite2 in aval2:
                                if akey.split('-')[0] == chrid: ## i.e. clevage site is on same chromosme
                                    print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite2,tarStrand,akey,aval,(aval2.index(cleaveSite2)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                    fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval2.index(cleaveSite2)+1-6))) ## +1 to convert to human readable and -6 to get in reference to phase site in phase file
                                    phasedCount += 1
                            else:
                                # print ('Cleave site matched phase index but chromosme was different')  
                                pass
            else:
                print("Input correct PARE result type in user settings")
                pass

        
    print('Total entries in PARE list:%s and matched:%s' % (totalCount,phasedCount))
    
    pass

##Connect to DB
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

##MAIN
def main():

    print("Step 1/3: Creating PHAS Dictionary")
    PHASdict = CreatePHASdict(PHASres,phase,pVal)
    print("Step 1/3 DONE!!\n\n")
    time.sleep(1) ## To improve user experience

    print("Step 2/3: Creating cleave Site list")
    resList,header = PAREreader(PAREres,PAREpval)
    print("Step 2/3 DONE!!\n\n")
    time.sleep(1) ## To improve user experience
    
    con = ConnectToDB(dataserver,0)
    print("Step 3/3: Matching cleave sites to phase sites")
    validatePHAS(resList,PHASdict,phase,con,header)
    print("Step 3/3:DONE!!\n\n")
    time.sleep(2)

if __name__ == '__main__':
    main()
    sys.exit()


################# LOG ######################################################
############################################################################

## v01 -> v02 ##Feb-9th
##Addition of category to MPPP results let to change in PAREreader indexes and output - Should be used on new result format
##Fixed if '_down' is encountered while making mysql querty to get the strand
## header for CL result from PAREredaer should be removed as added to latest version of CL itself

## v02 -> 03
##Comaptible to MPPPv095 output i.e cleave site index in PARE read was modified

## v03 -> v04 ## Major bug fix
## Cleavasite and chromsome are matched instead of just cleave site

## v04 -> v05
## Added functionality to vaildate phased loci on local trasncripts rather genome

## v05 -> v06
## Aded phase match capability for predicted targets
## Added two switches - target, and tarScore
## If transcript is the input for phasing analysis, make sure target prediction transcript name matches with phsing analysis.
## As of now, for PACBio - we use trans name for all analysis. So, full name isi split to match transname

## v06 -> v07
## Updated to accept sPARTA1.16 validated results that have targets for scaffodls with no annotated gene
## Added functionalty to match cleave sites to +1/-1 phase offset
## Added funtionality to use miRferno revmapped results for mode 0



### TO DO ###################################################################
## 1. The validation modes 0/1 with predicted/validated file type is redundant - Can be simplified like same loop for mode 0 or 1 - Needs free time 

################### README ##################################################
#############################################################################

## MODE 0
## Target file could be a combined "validated" file from genic and intergenic runs
## predicted targets could be used but need to be revmapped using revFerno.v02.py or above - Use Mode:0 , targetRes = 'C' and pareResType = 'T'
## Phased file should have genomic co-oridinates i.e. PHAS prediction done at whole genome level

## MODE 1
## Target file should correspond to transcripts generated from pacBIo or whatever
## Phase file must be done on same transcripts so in both files chr_id must be transcript name