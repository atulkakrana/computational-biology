#!/usr/local/bin/python3
##Script to annotate sRNA results fine and counts file with known miRs from mirbase
###Written by Atul Kakrana kakrana@udel.edu

import itertools as it
import sys
from operator import itemgetter

#### Settings ######
deResFile = 'ExactResCPM1.tsv' ## Diffrentially expressed results file
deCountsFile = 'miRCount.txt' 
mirbase = 'miRBASE_ZM.fa' ### file from miRBASE to be used as annotation

####################

###This function is to filter out known miRNAs from diffrentially regulated miRNAs file and counts file
def deMirProcess(deResFile,mirbase):
    
    ResFile = open('demiRsStats.csv','w')
    ResList = [] ## List to hold matched miRNAs results to fed in next function
    
    ###PART I Make dictionaries for tags - up,down or unregulated ############
    deTagsDict = {}
    ### Open file and sort on basis of logFC
    fh_in = open(deResFile,'r')
    deTagsHeader = fh_in.readline().strip('\n').replace('\t',',')
    detagsHeader_new = 'miRNA,Read,%s\n' % (deTagsHeader.rstrip(',')) ## One comma was present at end as there was a tab at end of header which got replaced to comma
    ResFile.write(detagsHeader_new) ##HEader added to results

    parsed_in = [line.strip('\n').split('\t') for line in fh_in]
    ##One line is ['>osa-miR812a', 'LOC_Os07g14310_up', '33472-33493', '33484', '1.5', '0.440755', 'CAGGUUGCAAAUUGGCAGGCAG', 'GUCUAACGUUUGACCGUCCGUC']
    parsed_in.sort(key=itemgetter(1))##Sorted on FC - Just in case required that way
    for i in parsed_in:
        #print(i[1])
        key = i[0].strip('"') ## Key is unique tag from tag count file
        #print (key,i)
        #break
        #logFC = float(i[1])
        deTagsDict[key] = (i[1],i[2],i[3],i[4])
    
    #print(len(upList),len(downList),len(unList))
    #print('Up-Regulated:',len(upDict),'Down-Regulated:',len(downDict),'Un-Regulated:',len(unDict))
    print('Entries in Dictionary:',len(deTagsDict))

    ###PART II Match with miRBASE miRNAs#####################################
    ##Match mirbase miRNAs with dictionaries:
    fh_in2 = open(mirbase,'r')
    miRead = fh_in2.read() ### Read altogether and not lines so that you can split
    miRs = miRead.split('>')
    acount = 0 ## Total miRNA count
    matchCount = 0
    for i in miRs[1:]:
        ent = i.split('\n')
        mirname = ent[0].split()[0]
        mirseq = ent[1].strip('\n').replace('U','T')
        acount+= 1
        
        ##Test with dicts
        if mirseq in deTagsDict.keys():
            print ('Found in Results: %s' % (mirname))
            results = deTagsDict[mirseq]
            #print (results)
            ResFile.write('%s,%s,%s\n' % (mirname,mirseq,','.join(results)))
            ResList.append((mirname,mirseq,results[0],results[1],results[2],results[3]))
            matchCount += 1
        else:
            pass
    print('Total entries analyzed in sRNA file: %s and total found in sRNA library: %s ' % (acount,matchCount))

    fh_in.close()
    ResFile.close()

    #for i in upDict.keys():
    #    print(i)

    return ResFile,ResList,detagsHeader_new

####This function goes through counts file and finds miRNAs from above to write a combined final file
def deMirCounts(ResList,deCountsFile,detagsHeader):
    
    ##main result file that combines stats and counts info
    combResult = open('demiRsCombined.csv','w')

    
    ##PART I: create dictionary of Counts file
    deCountsDict = {}
    ### Open file and sort on basis of logFC
    fh_in = open(deCountsFile,'r')
    deCountsHeader = fh_in.readline().strip('\n').replace('\t',',')
    ##Write combined header
    combResult.write('%s,%s\n' % (detagsHeader.strip('\n'),deCountsHeader.replace('"','')))

    parsed_in = [line.strip('\n').split('\t') for line in fh_in]
    ##One line is ['>osa-miR812a', 'LOC_Os07g14310_up', '33472-33493', '33484', '1.5', '0.440755', 'CAGGUUGCAAAUUGGCAGGCAG', 'GUCUAACGUUUGACCGUCCGUC']
    #parsed_in.sort(key=itemgetter(1))##Sorted on FC - Just in case required that way
    for i in parsed_in:
        #print(i[1])
        key = i[0].strip('"') ## Key is unique tag from tag count file
        #print (key,i)
        deCountsDict[key] = (i[1:]) ### Assign counts from all libraries as values
    
    
    ##PART II: Find counts for miRNAs identifed ealier (results file)
    for i in ResList:
        #print (i)
        mirseq = i[1]
        mirname = i[0]
        if mirseq in deCountsDict.keys():
            print ('Found in Counts: %s' % (mirname))
            #print ('These are the stats:' % ((i)))
            counts = deCountsDict[mirseq]
            #print ('These are the counts: %s' % (counts))
            combResult.write('%s,%s,%s,%s,%s,%s,%s\n' % (mirname,mirseq,i[2],i[3],i[4],i[5],','.join(counts)))
            
    
    combResult.close()
    
    return combResult



def main(deResFile,mirbase):
    
    ResFile,ResList,detagsHeader = deMirProcess(deResFile,mirbase) ## miRNA FASTA file divided on basis of DE
    combResult = deMirCounts(ResList,deCountsFile,detagsHeader)

if __name__ == '__main__':
    main(deResFile,mirbase)
    sys.exit()