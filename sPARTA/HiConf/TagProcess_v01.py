#!/usr/local/bin/python3

import itertools as it
import sys
from operator import itemgetter

#### Settings ######
deFile = 'ExactResCPM1.tsv'
mirbase = 'SoymiRNAAll.fa' ### file from miRBASe

dePARE = 'ExactResPARE_CPM05.tsv'### Diffrentially expressed PARE tags
countPARE = 'PAREcountCPM05.tsv' ### COunts file for all the tags in diffrentially expressd
####################

###This function is to firter out known miRNAs from diffrentially regulated miRNAs file and counts file
def deMirProcess(deFile,mirbase):
    
    ###PART I Make dictionaries for tags - up,down or unregulated ############
    upDict = {}
    downDict = {}
    unDict = {}
    ### Open file and sort on basis of logFC
    fh_in = open(deFile,'r')
    fh_in.readline()
    #deRead = fh_in.readlines()
    
    parsed_in = [line.strip('\n').split('\t') for line in fh_in]
    ##One line is ['>osa-miR812a', 'LOC_Os07g14310_up', '33472-33493', '33484', '1.5', '0.440755', 'CAGGUUGCAAAUUGGCAGGCAG', 'GUCUAACGUUUGACCGUCCGUC']
    parsed_in.sort(key=itemgetter(1))##Sorted on p-value beofre removing redundant so that first entry stored is best among other same interactions across library
    for i in parsed_in:
        #print(i[1])
        key = i[0].strip('"') ## Key is unique tag from tag count file
        #print (key,i)
        #break
        logFC = float(i[1])
        #print(repr(i))
        if logFC <= -1:
            #print('Downregulated: %s' % (i))
            #downList.append(i)
            downDict[key] = (i[1],i[2],i[3])
        elif logFC >= +1:
            #print('Upregulated: %s' % (i))
            #upList.append(i)
            upDict[key] = (i[1],i[2],i[3])
        else:
            #print('Un regulated: %s' % (i))
            #unList.append(i)
            unDict[key] = (i[1],i[2],i[3])
        
    #print(len(upList),len(downList),len(unList))
    print('Up-Regulated:',len(upDict),'Down-Regulated:',len(downDict),'Un-Regulated:',len(unDict))    

    
    ###PART II Match with miRBASE miRNAs#####################################
    ##Match mirbase miRNAs with dictionaries:
    fh_in2 = open(mirbase,'r')
    miRead = fh_in2.read() ### Read altogether and not lines so that you can split
    miRs = miRead.split('>')
    acount = 0 ## Total miRNA count
    upcount = 0 ## Found miRNA count
    downcount = 0
    uncount = 0
    upFile = open('upmiR.fa','w')
    downFile = open('downmiR.fa','w') 
    unFile = open('unmiR.fa','w')
    for i in miRs[1:]:
        ent = i.split('\n')
        mirname = ent[0].split()[0]
        mirseq = ent[1].strip('\n').replace('U','T')
        acount+= 1
        
        ##Test with dicts
        if mirseq in upDict.keys():
            print ('Found in up: %s' % (mirname))
            upFile.write('>%s_up\n%s\n' % (mirname,mirseq))
            #print(upDict[mirseq])
            upcount += 1
        elif mirseq in downDict.keys():
            print ('Found in down: %s' % (mirname))
            downFile.write('>%s_down\n%s\n' % (mirname,mirseq))
            #print(downDict[mirseq])
            downcount += 1
        elif mirseq in unDict.keys():
            unFile.write('>%s_un\n%s\n' % (mirname,mirseq))
            print ('Found in un-regulated: %s' % (mirname))
            #print(unDict[mirseq])
            uncount += 1
        else:
            pass
    
    upFile.close()
    downFile.close()
    unFile.close()
    fh_in.close()
        
    print('Total entries analyzed in sRNA file: %s and total found in sRNA library: %s ' % (acount,upcount+downcount+uncount))
    print('miRNA up: %s | miRNA down: %s | miRNAs unregulated: %s' % (upcount,downcount,uncount))

    #for i in upDict.keys():
    #    print(i)
    
    #return  upList,downList,unList
    return upFile, downFile, unFile

##This function separates tags in to three files 
def PAREProcess(dePARE,countPARE):
    
    ###PART I Make dictionaries for tags - up,down or unregulated ############
    ### Open file and sort on basis of logFC
    fh_in = open(dePARE,'r')
    fh_in.readline()
    parsed_in = [line.strip('\n').split('\t') for line in fh_in]
    parsed_in.sort(key=itemgetter(1))## Sorting on FC
    
    upDict = {} ## Dictionary to hold tag as key and other expression based figures as values for DE PARE tags
    downDict = {}
    unDict = {}
    
    for i in parsed_in:
        key = i[0].strip('"') ## Key is unique tag from tag count file
        #print(key)
        logFC = float(i[1])
        if logFC <= -1:
            #print('Downregulated: %s' % (i))
            downDict[key] = (i[1],i[2],i[3])
        elif logFC >= +1:
            #print('Upregulated: %s' % (i))
            upDict[key] = (i[1],i[2],i[3])
        else:
            #print('Un regulated: %s' % (i))
            unDict[key] = (i[1],i[2],i[3])
    print('Up-Regulated:',len(upDict),'Down-Regulated:',len(downDict),'Un-Regulated:',len(unDict))
    
    ####PART II Generate category specific (PARE) tag count files - up,down and un
    
    fh_in2 = open(countPARE,'r') ## Normalized counts per million file
    fh_in2.readline()  
    acount = 0 ## Total miRNA count
    upcount = 0 ## Found miRNA count
    downcount = 0
    uncount = 0
    upPARE = open('upPARE.txt','w')
    downPARE = open('downPARE.txt','w') 
    unPARE = open('unPARE.txt','w')
    
    for i in fh_in2:
        ent = i.strip('\n').split('\t')
        tag = ent[0].strip('"')## tag count file from R has "" as text qualifier
        #print(tag,i)
        acount+= 1
        
        if tag in upDict.keys():
            #print ('Found in up: %s' % (tag))
            upPARE.write('%s\t%d\n' % (tag,round(float(ent[2])))) ## Can also use larger value among ent[1] and ent [2] as that will be the original count
            upcount += 1
        elif tag in downDict.keys():
            #print ('Found in down: %s' % (tag))
            downPARE.write('%s\t%d\n' % (tag,round(float(ent[1])))) ##Can also use larger value among ent[1] and ent [2] as that will be the original count
            downcount += 1
        elif tag in unDict.keys():
            #print ('Found in unregulated: %s' % (tag))
            unPARE.write('%s\t%d\n' % (tag,round(float(max([ent[1],ent[2]]))))) ## Larger PARE tag value in whicever sampl
            uncount += 1
        else:
            print ('Tag from cout file not found in DE PARE tags - Problem')
            print('script exciting')
            sys.exit()
    
    print('Total entries analyzed in PARE DE file: %s and Matched with Count file: %s ' % (acount,upcount+downcount+uncount))
    print('PARE up: %s | PARE down: %s | PARE unregulated: %s' % (upcount,downcount,uncount))
    
    fh_in.close()
    fh_in2.close()
    upPARE.close()
    downPARE.close()
    unPARE.close()
    
    return upPARE,downPARE,unPARE
    


def main(deFile,mirbase):
    
    upFile, downFile, unFile = deMirProcess(deFile,mirbase) ## miRNA FASTA file divided on basis of DE
    upPARE, downPARE, unPARE = PAREProcess(dePARE,countPARE) ## PARE count file divided on basis of DE
    

if __name__ == '__main__':
    main(deFile,mirbase)
    sys.exit()
    


###LOG
