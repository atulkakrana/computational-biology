#!/usr/local/bin/python3
## Feb16-2014
##Script written by karana@udel.edu for Lachke's Lab
##Script for pvalue (qValues) calculation of isyTE ranks

################# LIBRARIES #######################

#import mysql.connector as sql

import re,sys,os,subprocess,csv,glob,time,random,math,shutil
import subprocess, multiprocessing
from operator import itemgetter
import sqlite3
from multiprocessing import Process, Queue, Pool
import numpy
import rpy2.robjects as R
from rpy2.robjects.packages import importr
import datetime

#################### INPUT #######################
allProbes = 'N' ## ENrichment scores to be generated for complte file or selective probes (requires testProbes)
testProbes = 'test.csv' ## tab seprated file with probe ID as first column and gene symbol as second column

chip = 2 ## AF = 1 and IL =2
annoFile = 'Illumina_WG6v2_annot_NonRed.txt' ### Clean Header lines before using here
resFile = 'IL_RESULTS_GLOB_P60_Batch_wo_P30C2_v2.txt'
resSep = '\t'
keyCol = 2 ## Depends upon key to be used and its column in result file

rankStat = 'TS' ## Possible options: 'FC' and 'TS'
#cols = [29,3] ## Fold Change -  Column number and number of tracks (total 10 tracks in current file)
#cols=[49,1] ## T-stats -  Col number and Number of tracks (total 10 tracks in current file)

#cols = [20,3] ## T-stats -Arota
#cols = [18,3] ## Tooth
cols = [39,8] ## Illumina T-stat 8 cols

nproc = '10' ## Number of processors - Possible options 'Y': 85% cores | 'N': 1 core | 'integer': Custom number of cores
randNum = 100 ## Random sampling to calculate qValue

#AA###TT############ FUNCTIONS ########UU#######LL##

def annoDictMaker(annoFile):
    
    '''Makes a dictionary of annotation file with probe_id as key and genomic co-ordinates
    as value. This dictionary is queried later for genome co-ordinates'''
    
    ##PARSE THE ANNOTATION FILE
    fh_in = open(annoFile, 'r')
    header = fh_in.readline()
    
    
    global annoDict ## Global for "relativeRankGenerator"
    annoDict = {}
    
    global annoProbeSET
    annoProbeSET = set() ## This will be used to randomly pick a probe for co-ordinates, if co-ords not found for a probe
    
    print('\nMaking a dictionary of probes names and co-ordinates from Annotation file')
    
    if chip == 1: ## Affy
        annoRead = csv.reader(fh_in)
        notFound = 0 ## a counter
        for line in annoRead:
            #print(line[0],line[12])
            if line[12] != '---': ## Not empty
                akey = line[0]
                chromosome,coords = line[12].split()[0].split(':')
                start,end = coords.split('-')
                
                ## Convert chromosomes to Numeric. Required for relativeRank carryover to 20MB region
                chromo = chromosome.replace('chr','')
                if re.match('^[0-99]{1,2}$',chromo):
                    #print('Numeric',chromo)
                    avalue = (int(chromo),int(start),int(end)) ##'chr' stripped from chromo
                    annoDict[akey] = avalue
                    annoProbeSET.add(akey)
                
                elif re.match('X',chromo):
                    #print ('X',chromo)
                    newchromo = 0 ## Turned to a universal numeric
                    avalue = (int(newchromo),int(start),int(end)) ##'chr' stripped from chromo
                    annoDict[akey] = avalue
                    annoProbeSET.add(akey)
                    
                else:
                    #print('Something else: %s' % (chromo)) ## FIlter out
                    pass
            else:
                notFound +=1
    
    elif chip == 2: ## chip is illumina
        annoRead = csv.reader(fh_in, dialect="excel-tab")
        annoList = [] ## To itnerate over file multiple times
        for i in annoRead:
            annoList.append(i)
        
        notFound = 0 ## a counter
        
        # Total chromosomes ## 
        chromoSet = set() ## Keeps all chromosmes
        for ent in annoList:
            if ent[4] != "NA": ## Gene start coord not empty
                chromo = ent[2]
                if re.match('^[0-99]{1,2}$',chromo):
                    #print(chromo)
                    chromoSet.add(int(chromo))   
        print(chromoSet)
        minChr = min(chromoSet)
        maxChr = max(chromoSet)
        print ('Min Chromo:%s | Max Chromo:%s' % (minChr,maxChr))
        ################
        
        for line in annoList:
            #print (line)
            #print(line[10],line[20])
            if line[4] != "NA": ## Not empty
                akey = line[0] ## Gene Symbol
                chromosome = line[2]
                start = line[4] ## THe co-ordinates are negative for negative strand
                end = line[5]
                
                ## Convert chromosomes to Numeric. Required for relativeRank carryover to 20MB region
                chromo = chromosome.replace('chr','')
                if re.match('^[0-99]{1,2}$',chromo):
                    #print('Numeric',chromo)
                    avalue = (int(chromo),abs(int(start)),abs(int(end))) ##'chr' stripped from chromo
                    annoDict[akey] = avalue
                    annoProbeSET.add(akey)
                
                elif re.match('X',chromo):
                    #print ('X',chromo)
                    newchromo = 0 ## Turned to a universal numeric
                    avalue = (int(newchromo),abs(int(start)),abs(int(end))) ##'chr' stripped from chromo
                    annoDict[akey] = avalue
                    annoProbeSET.add(akey)
                
                elif re.match('Y',chromo):
                    #print ('X',chromo)
                    newchromo = maxChr+1 ## Turned to a universal numeric
                    avalue = (int(newchromo),abs(int(start)),abs(int(end))) ##'chr' stripped from chromo
                    annoDict[akey] = avalue
                    annoProbeSET.add(akey)
                    
                else:
                    #print('Something else: %s' % (chromo)) ## FIlter out
                    pass
            else:
                notFound +=1
        pass
    
    else: ##
        print ('Please enter correct chip index in user settings')
        
    ## Test file
    fhOut = open('testAnno.tsv', 'w')
    for i in annoDict.keys():
        #print(i,','.join(str(x) for x in annoDict[i]))
        fhOut.write('%s,%s\n' % (i,','.join(str(x) for x in annoDict[i])))
    
    fhOut.close()
    fh_in.close()
    
    print('Length of annoDict:%s | probe with no coordinates: %s\n' % (len(annoDict),notFound))  
    
    return annoDict

def tableMaker(trackTableName,annoDict,resList,tableCols,conn):
    '''
    Make a track specific sqlite DB for probe ID and coords that will be used
    to query probes on 20MB interval. Each probe entry has following info:
    probe_id,FC,pval,chr,start,end
    
    '''
    #conn = sqlite3.connect(coordsDB)
    cur = conn.cursor()
    cur.execute('''DROP TABLE IF EXISTS %s''' % (trackTableName)) ### Drop Old table - while testing    
    #conn.commit()
    if rankStat == 'FC':
        cur.execute('''CREATE TABLE %s (probe varchar(255),fc decimal, pval decimal, chr integer, start integer, end integer)''' % (trackTableName))
        #conn.commit()
        for ent in resList:
            #print (ent[0])
            probe = ent[keyCol-1] ## Python format
            fc_col,pval_col = tableCols
            FC = ent[fc_col]
            pval = ent[pval_col]
            
            ## Get Coords - How are these used while ranking
            try:
                chromo,start,end = annoDict[probe]
            except KeyError:
                chromo = '99'
                start = 9999999999 ## 10 9s - Unobtainable in mouse
                end = 9999999999 ## 10 9s - Unobtainable in Mouse
            
            #print ('Probe:',probe,'| FC:',round(float(FC),2),'| PVAL:',round(float(pval),5),'| Chr:',chromo,'| Start:',start,'| End:',end)
            cur.execute("INSERT INTO %s VALUES ('%s',%f,%f,%d,%d,%d)" % (trackTableName,str(probe),round(float(FC),2),round(float(pval),5),int(chromo),int(start),int(end)))
            #conn.commit()
    
    elif rankStat == 'TS':
        cur.execute('''CREATE TABLE %s (probe varchar(255),tstat decimal, chr integer, start integer, end integer)''' % (trackTableName))
        #conn.commit()
        for ent in resList:
            #print (ent[0])
            probe = ent[keyCol-1] ## Python format
            ts_col,bogus = tableCols
            tstat = ent[ts_col] ##cols[0]
            
            ## Get Coords
            try:
                chromo,start,end = annoDict[probe]
            except KeyError:
                chromo = '99'
                start = 9999999999 ## 10 9s - Unobtainable in mouse
                end = 9999999999 ## 10 9s - Unobtainable in Mouse
            
            #print ('Probe:',probe,'| tstat:',round(float(FC),2),'| Chr:',chromo,'| Start:',start,'| End:',end)
            cur.execute("INSERT INTO %s VALUES ('%s',%f,%d,%d,%d)" % (trackTableName,str(probe),round(float(tstat),2),int(chromo),int(start),int(end)))
            #conn.commit()
    
    else:
        print ('''Please choose correct 'testStat'
               System will exit now''')
        sys.exit()
    conn.commit()

    ####Test
    #cur = conn.cursor()
    #testProbe = '1459888_at'
    #cur.execute("SELECT * FROM %s WHERE probe = '%s'" % (trackTable,testProbe))
    #print ('Here is your result for sqlite query:',cur.fetchall())
    
    #conn.close() ## Don't close as 'conn' is global variable and required in 'relativeRankGenerator'
    
    return trackTableName ## Return name as table will be acessed directly

def relativeRankGenerator(probe,trackTableName,chrLimit,conn):
    '''
    It will generate relative rank of provided probe, by finding flnaking genes in 20MB interval
    and ranking all genes on FC+pval or T-stats. Finally will report back the rank of the gene. Will
    use global variables - pval_col,fc_col,annoDict and resList
    '''
    relativeRank = 0 ##Initialize
    #conn = sqlite3.connect(coordsDB,check_same_thread=False) ## Open a connection
    
    ##Get input probe position
    #print ('\nRelative rank is being calculated for probe:%s' % (probe))
    try:
        probeCoord = annoDict[probe]
        #print ('Coordinates found for:%s - %s' % (probe,','.join(str(x) for x in probeCoord)))
    except KeyError:
        print ('''Coordinates not found for:%s
        #Random co-ordinates will be selected for relativeRank''' % (probe))
        subsProbes = randomSelect(resProbeSET,20)
        for i in subsProbes:
            if i in annoDict.keys():
                probeCoord = annoDict[i]
                probe = i
                print ('Randomly selected probe:',probe,'co-ordinate:%s\n\n' % (','.join(str(x) for x in probeCoord)))
                break
            else:
                continue         
        

    ## Get probes in 20MB interval
    cur = conn.cursor()
    #testProbe = '1426791_at'
    chromo = probeCoord[0]
    start = probeCoord[1] - 10000000
    end = probeCoord[2] + 10000000
    ##Test
    #cur.execute("SELECT count(*) FROM coords") ## Tested OK
    
    #print('Interval details - Probe:%s | Chr:%s | Start:%s | End:%s' % (probe,chromo,start,end))
    #cur.execute("SELECT * FROM coords WHERE chr = '%s'" % (chromo,))
    
    intervalProbes = [] ## List to hold all the probes in an interval
    
    if end >= chrLimit[chromo]:
        #print ('Probe class 2')
        
        carryOver = end-chrLimit[chromo]
        
        if chromo == max(list(chrLimit.keys())): ## Max available chromosome  than we cant add 1 to chromosme 
            carryChromo = min(list(chrLimit.keys())) ## Move around the chromosomes i.e 0,1,2,3,4....15 if O than move to 15
            #print ('Probe class 2.2')
        else:
            carryChromo = chromo+1
            #print ('Probe class 2.1')
        
        if rankStat == 'FC':
            cur.execute("SELECT probe,fc,pval FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY pval asc, fc desc" % (trackTableName,chromo,start,end))
            intProbes = cur.fetchall()
            intervalProbes.extend(intProbes)
            cur.execute("SELECT probe,fc,pval FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY pval asc, fc desc" % (trackTableName,carryChromo,0,carryOver))
            intProbes2 = cur.fetchall()
            intervalProbes.extend(intProbes2)
        
        else:
            cur.execute("SELECT probe,tstat FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY tstat desc" % (trackTableName,chromo,start,end))
            intProbes = cur.fetchall()
            intervalProbes.extend(intProbes)
            cur.execute("SELECT probe,tstat FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY tstat desc" % (trackTableName,carryChromo,0,carryOver))
            intProbes2 = cur.fetchall()
            intervalProbes.extend(intProbes2)
    
    elif start < 0 : ## i.e in negative
        #print ('Probe class 1')
        
        carryOver = abs(start) ## find co-ordinate after minus this length
        
        if chromo == min(list(chrLimit.keys())): ## That means you cant go to earlier chromosme, move to last chromosome i.e last chromosme
            carryChromo = max(list(chrLimit.keys()))
            #print ('Probe class 1.2')
        else:          
            carryChromo = chromo-1
            #print ('Probe class 1.1')
            
        carryLimit = chrLimit[carryChromo] ## Length of chromosome i.e last probe co-ordinates
        carryStart = carryLimit-carryOver ## carry start is 5' coordinate to the end to slice out region of carryOver
        

        if rankStat == 'FC':
            cur.execute("SELECT probe,fc,pval FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY pval asc, fc desc" % (trackTableName,chromo,0,end))
            intProbes = cur.fetchall()
            intervalProbes.extend(intProbes)
            cur.execute("SELECT probe,fc,pval FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY pval asc, fc desc" % (trackTableName,carryChromo,carryStart,carryLimit))
            intProbes2 = cur.fetchall()
            intervalProbes.extend(intProbes2)
        
        else:
            cur.execute("SELECT probe,tstat FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY tstat desc" % (trackTableName,chromo,0,end))
            intProbes = cur.fetchall()
            intervalProbes.extend(intProbes)
            cur.execute("SELECT probe,tstat FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY tstat desc" % (trackTableName,carryChromo,carryStart,carryLimit))
            intProbes2 = cur.fetchall()
            intervalProbes.extend(intProbes2)
    
    else: ## All is well
        #print ('Probe class 0')
        
        if rankStat == 'FC':
            cur.execute("SELECT probe,fc,pval,chr,start FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY pval asc, fc desc" % (trackTableName,chromo,start,end))
            intProbes = cur.fetchall()
            #print(intProbes)
            intervalProbes.extend(intProbes)
        else:
            cur.execute("SELECT probe,tstat,chr,start FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY tstat desc" % (trackTableName,chromo,start,end))
            intProbes = cur.fetchall()
            #print(intProbes)
            intervalProbes.extend(intProbes)
            
    #print ('Interval ranks for query:',probe,': ',intervalProbes)
    
    
    ## Find rank of probe in sorted list
    aRank = 1
    for i in intervalProbes:
        if i[0] == probe:
            relativeRank = aRank
        else:
            aRank += 1
    #print('Relative Rank of %s:%d' % (probe,relativeRank))
    
    weigRelativeRank = relativeRank/len(intervalProbes)
    #print('Weighted Relative Rank for %s: %s'% (probe,str(weigRelativeRank)))
    
    ## Error Check
    if weigRelativeRank == 0:
        #print('weighted Relative rank is 0 - Why?')
        cur.execute("SELECT probe,tstat,chr,start FROM %s WHERE chr = '%s' AND start between %s and %s ORDER BY chr asc, start asc" % (trackTableName,chromo,start,end))
        intProbes = cur.fetchall()
        #print ('Interval members:',probe,':',intProbes)
        sys.exit()
    else:
        pass

    #conn.close()
    #return relativeRank ## Integer
    
    return weigRelativeRank
   
def pvalCalcAll(col):
    '''
    Calculates P-value for every track
    '''
    DB = '%s_%s.db' % (coordsDB.split('.')[0],col)
    conn = sqlite3.connect(DB)
    
    if rankStat == 'FC':
        
        trackPvals = [] ## Temporary list that will hold results for single track
        
        ### FIND pval COLUMN ###
        if head[col].split('.')[0] == 'FC':  ## Confirms FC as column               
            
            pval_title = ('p.value.%s' % head[col].split('.')[-1])   ## Find Matching p_value column
            print('\n*****************TRACK: %s***********************'  % (head[col]))
            print('\nTitle to match: %s and FC column: %s' % (pval_title, head[col]))
            pval_col = head.index(pval_title)   ## pval column to sort
            tableCols = [col,pval_col]
    
            ## Make sqlite3 table for current track -  probe,FC,pval,chr,start,end
            trackTableName = 'track_%s' % (col)
            trackTableName = tableMaker(trackTableName,annoDict,resList,tableCols,conn) ##Makes the table and resturn the name itself - Does't matter
            chrLimit = chrLimits(trackTableName,conn)
            
            ## Itnerate every gene of current track - Do all processing
            for ent in resList: ##
                #print(ent[0])
                ## Generate Relative Rank for gene of interest
                mainProbe = ent[0] ## Probe to rank
                geneSym = ent[1]
                mainRank = relativeRankGenerator(mainProbe,trackTableName,chrLimit,conn)
                #print('Gene Symbol:',geneSym,'| mainProbe:',mainProbe,'| mainRank:',mainRank)
                
                randomProbes = randomSelect(resProbeSET,randNum) ## randomProbes has to be alist
                
                ##Generate relative ranks for Random s

                randomRanks = []
                for randProbe in randomProbes:
                    randomRank = relativeRankGenerator(randProbe,trackTableName,chrLimit,conn)
                    randomRanks.append(randomRank)
                
                ### How many times mainProbeRank was more than randomProbeRanks - Calculate p-value
                randRankArray = numpy.array(randomRanks)
                probRR = ((randRankArray<=mainRank).sum())/randNum
                #print('\nprobRR:%f\n' % (probRR))
                trackPvals.append(round(probRR,5))
                
        else:
            print(''' Please check the column specified: %s
                  Script will exit now''' % (head[fc_col]))
            sys.exit()
 
    
    elif rankStat == 'TS':

        trackPvals = [] ## Temporary list that will hold results for single run
        
        ### Confirm Column ###
        if head[col].split('.')[0] == 't':  ## Confirms FC as column
            tableCols = [col,'NA'] ## Column to name trackTable and use for relativeRank calculation
            print('\n*****************TRACK: %s***********************'  % (head[col]))
            
            ## Make sqlite3 table for current track -  probe,FC,pval,chr,start,end
            trackTableName = 'track_%s' % (col)
            trackTableName = tableMaker(trackTableName,annoDict,resList,tableCols,conn) ##Makes the table and resturn the name itself - Does't matter
            chrLimit = chrLimits(trackTableName,conn)
            
            ## Itnerate every gene of current track - Do all processing
            for ent in resList: ###
                
                ## Generate Relative Rank for gene of interest
                mainProbe = ent[0] ## Probe to rank
                geneSym = ent[1]
                mainRank = relativeRankGenerator(mainProbe,trackTableName,chrLimit,conn)
                #print('Gene Symbol:',geneSym,'| mainProbe:',mainProbe,'| mainRank:',mainRank)
                
                randomProbes = randomSelect(resProbeSET,randNum) ## randomProbes has to be alist
                
                ##Generate relative ranks for Random set

                randomRanks = []
                for randProbe in randomProbes:
                    randomRank = relativeRankGenerator(randProbe,trackTableName,chrLimit,conn)
                    randomRanks.append(randomRank)

                #randomRanks = PPResults(relativeRankGenerator,randomProbes)
                #print('mainRank:',mainRank,'| randomRanks',randomRanks)
                
                ### How many times mainProbeRank was more than randomProbeRanks - Calculate p-value
                randRankArray = numpy.array(randomRanks)
                probRR = ((randRankArray<=mainRank).sum())/randNum
                #print('\nprobRR:%f\n' % (probRR))
                trackPvals.append(round(probRR,5))
                
        else:
            print(''' Please check the column specified: %s
                  Script will exit now''' % (head[col]))
            sys.exit()


    else: ## USer inputs wrong testStat
        print(''' Please check the 'rankStat' specified: %s
              Script will exit now''' % (rankStat))
        sys.exit()
    
    conn.close()
    
    trackDict = {}
    trackDict[col] = trackPvals ##Dictionary for relatiev ranks of all probes within a track/col
    
    return trackDict

def pValCalcTest(col):
    
    '''
    Calculates P-value for every track
    Input: column number to use for calculation of p-values and testProbeList (global variable)
    '''
    DB = '%s_%s.db' % (coordsDB.split('.')[0],col)
    conn = sqlite3.connect(DB)
    
    if rankStat == 'FC':
        
        trackPvals = [] ## Temporary list that will hold results for single track
        
        ### FIND pval COLUMN ###
        if head[col].split('.')[0] == 'FC':  ## Confirms FC as column               
            
            pval_title = ('p.value.%s' % head[col].split('.')[-1])   ## Find Matching p_value column
            print('\n*****************TRACK: %s***********************'  % (head[col]))
            print('\nTitle to match: %s and FC column: %s' % (pval_title, head[col]))
            pval_col = head.index(pval_title)   ## pval column to sort
            tableCols = [col,pval_col]
    
            ## Make sqlite3 table for current track -  probe,FC,pval,chr,start,end
            trackTableName = 'track_%s' % (col)
            trackTableName = tableMaker(trackTableName,annoDict,resList,tableCols,conn) ##Makes the table and return the name itself - Does't matter
            chrLimit = chrLimits(trackTableName,conn)
            
            ## Itnerate every gene of current track - Do all processing
            for ent in testProbeList: ##
                print(ent)
                
                ## Generate Relative Rank for gene of interest
                mainProbe = ent[0] ## Probe to rank
                geneSym = ent[1] ## Gene Symbol
                
                ##Check if probe exists in the microrray results file - as studies used for cross validation might not have some probes
                if mainProbe not in resDict:
                    print('\n\n main probe misisng from microarray result file \n\n')
                    trackPvals.append(round(0.999999,5)) ## main probe is missing in this microaary experiment - so not enriched
                    continue ## Test other main probe 
                elif mainProbe not in annoDict:
                    trackPvals.append(round(0.999999,5)) ## main probe is missing in this microaary experiment - so not enriched
                    print('\n\n main probe missing annotations \n\n')
                    continue ## Test other main probes
                else:
                    pass
                
                #conn = sqlite3.connect(coordsDB) ## Open a connection
                
                mainRank = relativeRankGenerator(mainProbe,trackTableName,chrLimit,conn)
                print('Gene Symbol:',geneSym,'| mainProbe:',mainProbe,'| mainRank:',mainRank)
                
                randomProbes = randomSelect(resProbeSET,randNum) ## randomProbes has to be alist
                
                ##Generate relative ranks for Random s

                randomRanks = []
                for randProbe in randomProbes:
                    randomRank = relativeRankGenerator(randProbe,trackTableName,chrLimit,conn)
                    randomRanks.append(randomRank)
                
                
                ### How many times mainProbeRank was more than randomProbeRanks - Calculate p-value
                randRankArray = numpy.array(randomRanks)
                probRR = ((randRankArray<=mainRank).sum())/randNum
                #print('\nprobRR:%f\n' % (probRR))
                trackPvals.append(round(probRR,5))
                
        else:
            print(''' Please check the column specified: %s
                  Script will exit now''' % (head[fc_col]))
            sys.exit()
 
    
    elif rankStat == 'TS':

        trackPvals = [] ## Temporary list that will hold results for single run
        
        ### Confirm Column ###
        if head[col].split('.')[0] == 't':  ## Confirms FC as column
            tableCols = [col,'NA'] ## Column to name trackTable and use for relativeRank calculation
            print('\n*****************TRACK: %s***********************'  % (head[col]))
            
            ## Make sqlite3 table for current track -  probe,FC,pval,chr,start,end
            trackTableName = 'track_%s' % (col)
            trackTableName = tableMaker(trackTableName,annoDict,resList,tableCols,conn) ##Makes the table and resturn the name itself - Does't matter
            chrLimit = chrLimits(trackTableName,conn)
            
            ## Itnerate every gene of current track - Do all processing
            #DB = '%s_%s.db' % (coordsDB.split('.')[0],i)
            #conn = sqlite3.connect(DB) ## Open a connection
            
            for ent in testProbeList: ###
                
                print('\nMain probe for enrichment p-value: %s' % (ent))
                
                ## Generate Relative Rank for gene of interest
                mainProbe = ent[0] ## Probe to rank
                geneSym = ent[1]
                
                if mainProbe not in resDict:
                    print('\n\n main probe misisng from microarray result file \n\n')
                    trackPvals.append(round(0.999999,5)) ## main probe is missing in this microaary experiment - so not enriched
                    continue ## Test other main probe 
                elif mainProbe not in annoDict:
                    trackPvals.append(round(0.999999,5)) ## main probe is missing in this microaary experiment - so not enriched
                    print('\n\n main probe missing annotations \n\n')
                    continue ## Test other main probes
                else:
                    pass
                
                
                mainRank = relativeRankGenerator(mainProbe,trackTableName,chrLimit,conn)
                print('Gene Symbol:',geneSym,'| mainProbe:',mainProbe,'| mainRank:',mainRank)
                
                randomProbes = randomSelect(resProbeSET,randNum) ## randomProbes has to be alist
                
                ##Generate relative ranks for Random set

                randomRanks = []
                for randProbe in randomProbes:
                    randomRank = relativeRankGenerator(randProbe,trackTableName,chrLimit,conn)
                    randomRanks.append(randomRank)

                #randomRanks = PPResults(relativeRankGenerator,randomProbes)
                #print('mainRank:',mainRank,'| randomRanks',randomRanks)
                
                
                ### How many times mainProbeRank was more than randomProbeRanks - Calculate p-value
                randRankArray = numpy.array(randomRanks)
                probRR = ((randRankArray<=mainRank).sum())/randNum
                #print('\nprobRR:%f\n' % (probRR))
                trackPvals.append(round(probRR,5))
            
            
            
        else:
            print(''' Please check the column specified: %s
                  Script will exit now''' % (head[col]))
            sys.exit()


    else: ## USer inputs wrong testStat
        print(''' Please check the 'rankStat' specified: %s
              Script will exit now''' % (rankStat))
        sys.exit()
    
    conn.close()
    
    trackDict = {}
    trackDict[col] = trackPvals ##Dictionary for relative ranks of all probes within a track/col
    
    return trackDict

def trackManager(resFile,cols,nproc,annoDict):
    
    ### Prepare ###############################################################################
    '''
    This will go through every track one by one, and each gene of a track.
    relative rank of gene will be computed folllowed by selecting a random set of
    1000 genes. Relative rank of these 100 genes will be calculated. Finally p(rr)
    computed  
    '''
    
    ## Read Microarray results file once 
    fh_in = open(resFile,'r')
    global resList,head ## Global for "relativeRankGenerator" - Need better method - So many globals
    res = [line.strip('\n').split(resSep) for line in fh_in] ## Line converted to list and Read
    resList = list(res) ## List to hold results while pre-processing
    head = resList.pop(0) ## Remove header before sorting and store here
    
    ## To fetch details of probes of interest or genes in the end
    global resDict
    resDict = {}
    for i in resList:
        key = i[keyCol-1] ## Python format
        resDict[key] = i

    ## If enrichment q-val for only few probes need to calcualated
    if allProbes == 'N':
        global testProbeList ## for p-value calculation
        fh_in2 = open(testProbes, 'r')
        fileRead = [line.strip('\n').split('\t') for line in fh_in2]
        testProbeList = list(fileRead)
        #print('\nInput probes list:',testProbeList)

    ## Create probeSET - Used for random selection of probes
    global resProbeSET ## for p-val calculation
    resProbeSET = set()
    for ent in resList:
        probe = ent[keyCol-1] ## Python format
        resProbeSET.add(probe)

    print('\nTotal probes in microarray results file that were added to the SET: %s' % (len(resProbeSET)))

    #### ANALYZE ###################################################################
    
    ### All global variables to be used for relativeRanking and other stuff
    global trackTableName,coordsDB,chrLimit  ## Global for "relativeRankGenerator"
    
    ## Collect track columns #############################################
    startCol = int(cols[0])-1 ## Coverted to python format
    trackList = []
    for i in range(cols[1]):
        trackList.append(startCol)
        startCol += 1
    #print(trackList)

    ## Make DB ###################
    coordsDB = 'tracksInfo.db' ## Name to use while making DBs
    ###############################
    
    ### Get p-vals ###############
    
    if allProbes == 'N':
        allTrackPvals = []
        ## Serial
        for i in trackList:
            trackPvals = pValCalcTest(i)
            allTrackPvals.append((trackPvals))
        
        ### Parallel
        #allTrackPvals = PPResults(pValCalcTest,trackList) ## How to sort results on basis of tracks - Perhaps generate results in form of dictionary entry with track name as key
        #print('allTrackPvals:',allTrackPvals)
    
    else:
        allTrackPvals = []
        
        ### Serial
        #for i in trackList:
        #    trackPvals = pValCalcAll(i)
        #    allTrackPvals.append((trackPvals))
        
        ## Parallel
        allTrackPvals = PPResults(pValCalcAll,trackList) ## How to sort results on basis of tracks - Perhaps generate results in form of dictionary entry with track name as key
        #print('allTrackPvals:',allTrackPvals)
        
    
    ## Create one main dict with tracks cols as key and track results as value - Will be used for resorting results trackwise before merging
    allTrackPvalsDict = {}
    for i in allTrackPvals:
        allTrackPvalsDict.update(i)
    #print('allTRackPvalsDict:',allTrackPvalsDict)
    
    ###Sort the results track wise - Just to be sure
    allTrackPvalsSorted = [] ## Have track sorted list for merging
    for tracks in trackList:
        trackPvals = allTrackPvalsDict[tracks]
        allTrackPvalsSorted.append(trackPvals)
    #print('\nTrack sorted results:',allTrackPvalsSorted)

    ##Transpose pvals i.e. pvals for a probeset in a tuple
    mergedAllPvals = zip(*allTrackPvalsSorted)
##----------------------------------------------------------------------------------------------------
    #### OUTPUT ##########################################################################################
    
    ## Get header to write results
    pvalsHeader = [] ## Grab header for iSyTE2 p-vals
    for track in trackList:
        pvalsHeader.append((  'pval.%s' % (head[track].partition('.')[-1])   ))
    
    ##Open Result File
    resultFile = 'Testpvals_rand%s_track%s_%s.tsv' % (randNum,cols[1],datetime.datetime.now().strftime("%m_%d_%H_%M"))
    fh_out = open(resultFile,'w')
    fh_out.write('%s\t%s\tcomb.qval\n' % ('\t'.join(head),'\t'.join(pvalsHeader)))
    
    ### Calculate qValues
    if allProbes == 'N':
        nColumns = len(list(resDict.values())[0]) ## Number of columns in microaaray file, to use when probe is not present in microaaray file
        #print('nColumns',nColumns)
        
        for ent,resEnt in zip(mergedAllPvals,testProbeList):
            
            ## Check if probe exists in microarray results
            if resEnt[0] not in resDict:
                geneInfo = ['Absent']*nColumns
            elif resEnt[0] not in annoDict:
                geneInfo = ['No Coords']*nColumns
            else:
                geneInfo = resDict[resEnt[0]]
            #print('pvals:',ent,'\nResult Entry:',geneInfo)
            ## Pval Comb
            y = numpy.array(ent) ## Make numpy array
            stat = (-2*numpy.log(y+0.000001)).sum() ### Added '0.000001' to avoid division by zero errors
            npvals = int(len(ent))
            #print('Calculating qval\n')
            combqval = R.r['pchisq'](stat,df=2*npvals,lower = False) ## Add lower.tail = FALSE
            #print('Qval',combqval)
            
            ## Just for the sake of it
            if combqval == 0:
                combqval = 0.00001
            elif combqval == 1:
                combqval = 0.99999
            else:
                pass
            
            fh_out.write(  '%s\t%s\t%s\t%s\t%f\n' % ( resEnt[0],resEnt[1],'\t'.join(geneInfo[2:]), '\t'.join(format(x, "1.4f") for x in ent), round(float(combqval[0]),5) ) )
            
    else:
        for ent,resEnt in zip(mergedAllPvals,resList):
            #print('pvals:',ent,'\nResult Entry:',resEnt)
            ## Pval Comb
            y = numpy.array(ent) ## Make numpy array
            stat = (-2*numpy.log(y+0.000001)).sum() ### Added '0.000001' to avoid division by zero errors
            npvals = int(len(ent))
            #print('Calculating qval\n')
            combqval = R.r['pchisq'](stat,df=2*npvals,lower = True) ## Add lower.tail = FALSE
            
            ## Just for the sake of it
            if combqval == 0:
                combqval = 0.00001
            elif combqval == 1:
                combqval = 0.99999
            else:
                pass
            
            fh_out.write(  '%s\t%s\t%f\n' % ( '\t'.join(resEnt), '\t'.join(format(x, "1.4f") for x in ent), round(float(combqval[0]),5) ) )
    
    fh_out.close()
    
    return resultFile

def chrLimits(trackTable,conn):
    ''' Computes the end limits for every chromosome
    and gives back as a dictionary'''
    
    chrLimits = {} ## Chr:Limit
    
    #conn = sqlite3.connect(coordsDB)
    cur = conn.cursor()
    cur.execute('''SELECT distinct(chr) FROM %s''' % (trackTable)) ### Drop Old table - while testing
    allChrs = cur.fetchall()
    #print ('Here is your result for chromosomes:',allChrs)
    
    for achr in allChrs:
        #print ('This is the chromosme: %s' % achr[0])
        cur.execute('''SELECT chr,MAX(end) FROM %s WHERE chr = %s''' % (trackTable,int(achr[0]))) ### Drop Old table - while testing
        chrRes = cur.fetchall()
        print ('Here is your result for chromosome limits:',chrRes[0][0],chrRes[0][1])
        chrLimits[int(chrRes[0][0])] = int(chrRes[0][1]) ## Add to dictionary
    #conn.close()
    
    return chrLimits

def randomSelect(probeSET,size):
    '''
    A specifed number of genes selected randomly and reported back
    '''
    randomProbes = random.sample(probeSET,int(size))
    #print(randomProbes)
    
    return randomProbes
    
def PPResults(module,alist):
    npool = Pool(int(nproc))    
    res = npool.map_async(module, alist)
    results = (res.get())
    return results

def main():
    annoDict = annoDictMaker(annoFile)
    resultFile = trackManager(resFile,cols,nproc,annoDict)
    
    ##Clean up
    print ("**Cleaning temp files last time**")
    garbage = [file for file in os.listdir('./') if file.endswith (('.db','.xyz'))] ## Excluded 'chopped.trimmed.fastq' as used by RNA Runner
    for file in garbage:
        print("Deleting %s" % (file))
        os.remove(file)

if __name__ == '__main__':

    if nproc == 'Y':
        nproc = int(multiprocessing.cpu_count()*0.85)
    elif nproc == 'N':
        nproc = int(1)
    else:
        nproc = int(nproc)
    ##
    start = time.time()
    main()
    end = time.time()
    print ("\nThe complete run time is %s" % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting..\n')
    sys.exit()
    
###Modification make a merged DB for every track with probe,coords,FC and pval - Use this to rank

##v01- Info
## Connection if opened and closed ine ach module so as to maintain separate copied of connection in parrallel computing
##v01 -> v02
##Added fishers method for p-value combination

##v03 -> v04
##Add functionality to genrate enricment-qvals for selctive probes -Done
## Improve Ranking mechanism by weigted Ranks - Done

##v04 -> v05
## Fixed: Error when probeCoords can't be found, a probe from annotation was file was being picked, but there is gaurantee that this probe will be present in microarray result,
## So, if the new probe is absent in results, it will not be in rankInterval and its rank will be 0 always - Fixed by picking random probe from resFile
## Fixed: Paralle processing failure due to connections to same sqlite DB by multiple cores - Replicated DB number of times track - SO every track has its own DB and connection

##v05 -> v06
## Added Illumina support
## Fixed a few things like 'Y' chr is the max chromosome

