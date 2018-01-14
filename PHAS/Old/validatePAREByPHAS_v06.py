#!/usr/local/bin/python3
##V02 has modification for new MPPP format i.e with category column

import sys,time,os,glob
import mysql.connector as sql

#### Settings
genomeDB = 'RICE_MSU7_genome' ## Used to identify strand of target
dataserver = 'raichu.dbi.udel.edu'

mode = 1                        ## 0: Target coordinates match genome DB 1: transcripts or non-genome-DB targets (local analysis)
target = 'T'                    ## 'T' is transcripts 'C' if genomic i.e. what is in the chromosome column of phaseing results - transcript name or chromosome

PAREpval = 0.25                 ## Cutoff of corrected p-value to use to filter out non-relevant predictions
tarScore = 7                    ## Only if PAREresType == 'T'
PAREresType = 'T'               ## 'S' for sPARTA PARE validated, 'C' for CL3 PARE Validated and 'T' for just sPARTA target prediction file
PAREres = 'All.targs.parsed.csv' ## Mode = 0 revmapped file with 14th column as cleave site | Mode: 1 -local validation file with cleave site at column 9

phase = 21                      ## Phase length that you are using for double validation - Give corresponding file below
PHASresType = 'PP'              ## J: jixian results format | P:Ping list file | PP: Phas redundant script generated *.ping format
PHASres = 'Final_PHASLociID_1e-07_21ALL.ping'
pVal = '1e-07'                  ## P-value of phased loci to test for


##Read phased loci file and make dictionary
def CreatePHASdict(PHASres,phase,pVal):
    ''' Reads the phased result file
    and creates a dictionary of -5/+7 phased locations as value
    '''
    
    fh_in = open(PHASres,'r')
    #fh_out = open('%s_Dict' % (PHASres),'w')
    PHASdict = {}
    entries = fh_in.readlines()
    
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
                print(phase,pval,chromo,start,end,trash1,trash2)
                if target == 'T':
                    akey ='%s-%s-%s' % (chromo.strip().split("|")[0],start,end) ## Since its local analysis, transcrpt name would be 
                else:
                    akey ='%s-%s-%s' % (chromo.strip().replace('chr',''),start,end)
                aval = [sum(i) for i in zip([int(start)]*11,phaseList)]
                if pval == pVal:
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
    
    resList = [] ## miRNA,Target,cleavesite,whole entry
    totalCount = 0
    if PAREresType == 'M' and mode == 0: ## MPPP server
        for i in entries:
            #print (i)
            ent_splt = i.split(',')
            #print(ent_splt[13])
            if float(ent_splt[13]) <= PAREpval:
                cleaveSite = int(ent_splt[18])
                #print(cleaveSite)
                tarName = ent_splt[1]
                mirName = ent_splt[0]
                chrid = ent_splt[15]
                strand = ent_splt[16]
                resList.append((mirName,tarName,cleaveSite,chrid,i.strip('\n')))
                #print(mirName,tarName,cleaveSite)
                totalCount +=1

    elif PAREresType == 'M' and mode == 1: ## MPPP Local
        for i in entries:
            #print (i)
            ent_splt = i.split(',')
            #print(ent_splt[13])
            if float(ent_splt[14]) <= PAREpval:
                cleaveSite = int(ent_splt[8])
                #print(cleaveSite)
                tarName = ent_splt[1]
                mirName = ent_splt[0]
                chrid = ent_splt[1] ## In local mode there is no chr_id, only transcripts
                strand = 'None'
                resList.append((mirName,tarName,cleaveSite,chrid,i.strip('\n')))
                #print(mirName,tarName,cleaveSite)
                totalCount +=1

    elif PAREresType == 'T' and mode == 1: ## sPARTA local targets only, no PARE
        for i in entries:
            #print (i)
            ent_splt = i.strip('\n').split(',')
            print(ent_splt)
            if float(ent_splt[5]) <= tarScore:
                bindStart,bindEnd = ent_splt[2].split("-")
                cleaveSite1 = int(bindEnd)-10 ## 10th pos
                cleavSite2 = int(bindEnd)-11 ## 11th pos
                cleaveSite = (cleaveSite1,cleavSite2)
                print(cleaveSite)
                tarName = ent_splt[1]
                mirName = ent_splt[0]
                chrid = ent_splt[1] ## In local mode there is no chr_id, only transcripts
                strand = 'None'
                resList.append((mirName,tarName,cleaveSite,chrid,i.strip('\n')))
                #print(mirName,tarName,cleaveSite)
                totalCount +=1

    
    elif PAREresType == 'C': ## CL3
        for i in entries:
            #print (i)
            ent_splt = i.split(',')
            if float(ent_splt[9]) <= PAREpval:
                mirName = ent_splt[0]
                tarName = ent_splt[1]
                cleaveSite = int(ent_splt[5])
                print (cleaveSite)
                resList.append((mirName,tarName,cleaveSite,i.strip('\n')))
                totalCount +=1
        header = ('miR,TarName,Chr,Strand,bindSite,cleaveSite,Score,mirSeq,tarSeq,p-val,Small,Large,Ratio' % ()) ## as missing from sco_inp_ext
    
    print('\n\nTotal entries passed p-val i.e in list:%s made from file:%s' % (totalCount,PAREres))
    
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
        if PAREresType == 'T':
            cleaveSite1,cleaveSite2 = ent[2] ## Predicted 10th and 11th as cleave site
        else:
            cleaveSite = int(ent[2])    ## PARE supported cleave site
        chrid = ent[3]
        #print(mirName,tarName,cleaveSite)
        
        if mode == 0: ## targets coordinates match to genome DB i.e. server target analysis was done
            if PAREresType == 'S' or PAREresType == 'C':
                cur= con.cursor()
                #print('Query DB:',genomeDB,'| Query Target:',tarName.split('_up')[0].split('_down')[0])
                cur.execute("SELECT strand FROM %s.gene_master where gene like '%s'" % (genomeDB,tarName.split('_up')[0].split('_down')[0]))### Convert intergenic to gene name so as to get strand
                tarStrand = cur.fetchall()
                #print (tarStrand)
        
                if tarStrand[0][0] == 'c': ##Add 3 to cleave site
                    #print("\nCrick strand:%s" % (tarStrand[0][0]))
                    cleaveSite += 3    
                else: ## Do nothing to cleave site
                    #print("\nWatson strand:%s" % (tarStrand[0][0]))
                    pass
                
                ##PHAS LOCI MATCH - Fuck you -  Was not checking for cleave site in exact chromosome - Thats why you got shit load of results where these are scaffolds
                for akey in PHASdict.keys():
                    aval = PHASdict[akey]
                    #print ('This is the key:', akey, aval)
                    if cleaveSite in aval:
                            if akey.split('-')[0] == chrid: ## I.e. clevage site is on same chromosme
                                print('\nmiRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d' % (mirName,tarName,cleaveSite,tarStrand[0][0],akey,aval,(aval.index(cleaveSite1)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                                fh_out.write('%s,%s,%s\n' % (ent[4],akey,(aval.index(cleaveSite)+1-6)))
                                phasedCount += 1
                            else:
                                print ('Cleave site matched phase index but chromosme was different')            
                    else:
                        pass

            elif PAREresType == 'T':
                print("Server mode is not supported with target prediction results, simply put this mode is not ready yet")
                sys.exit()

            else:
                print("Input correct PARE result type in user settings")



        elif mode == 1: ## Local analysis - Targets are transcripts and no chr_id or strand is required from server

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
                        pass

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
                        pass

            else:
                print("Input correct PARE result type in user settings")

        
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
    PHASdict = CreatePHASdict(PHASres,phase,pVal)
    resList,header = PAREreader(PAREres,PAREpval)
    
    con = ConnectToDB(dataserver,0)
    validatePHAS(resList,PHASdict,phase,con,header)

if __name__ == '__main__':
    main()
    sys.exit()

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