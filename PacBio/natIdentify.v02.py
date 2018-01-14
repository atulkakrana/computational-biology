#!/usr/local/bin/python3
## Written by Atul Kakrana - Kakrana@udel.edu

## This script requires summary file from PACBIO as generated including results from geneOverlap.py
## This will check for type of overlap and identify NAT and spit out the configuration

## IMPORTS
import os,sys

## USER SETTINGS

summary = "summary.agp2.v3.txt"               ### Should be tab seprated file, because overlap columns have comma seprated values                        
head = 'Y'                                   ## Y: Header present in summary file | N: No header present in sumamry file
genepos =          1                         ## PacBio gene ID, generated after collapse - PB.1 i.e. without transcript part
transpos =         2                         ## Coulmn number for transcript ID - PB.1.1, not for gene ID
codStatus =         27                       ## Column number for PLEK status i.e. "Coding" and "Non-coding"
overlapGenes =      30 
overlapFlags =      31
overlapConf =       32                      ## Column number for overlap configuration
overlapStatus =     33                      ## Column number for overlap status
matOrientPos    =   18
matScorePos    =    19
matgenePos =        8                       ## Used to check if Transcript is overlapping the same gene as from matchAnnot


gtf = "../Zea_mays.AGPv3.27.gtf"


def NATdetect(summary,summDict,geneDict):
    ''' Identify NAT, type and relation between gene-trans'''

    print("\nFunction: NATdetect")

    outFile = "NATReport.txt"
    fh_out = open(outFile,'w')
    fh_out.write("Gene\tTranscript\tPartner\tTransStaus\tpartnerStatus\toverFlag\tNAT\tClass\n")

    fh_in = open(summary,'r')
    if head == 'Y':
        fh_in.readline()
    summRead = fh_in.readlines()

    ## Combine coding status dictionaries - So that if transcript overlaps with either gene or other trans, same dict could be used
    statusDict = dict(list(summDict.items()) + list(geneDict.items()))

    acount = 0          ## All transcripts
    bcount = 0          ## All overlapping
    ccount = 0          ## All natural antisense
    dcount = 0          ## All natural sense

    for i in summRead:
        acount+=1
        ent = i.strip("\n").split("\t")
        # print(ent)
        gene        = ent[genepos-1]
        transcript  = ent[transpos-1]
        transinfo   = statusDict[transcript]
        transStatus = transinfo[1]
        matScore    = ent[matScorePos-1]
        matGene     = ent[matgenePos-1]
        matOrient   = ent[matOrientPos-1]

        if ent[overlapGenes-1] != "NA":
            print("\n")
            print("Gene:%s | Trans:%s | Matgene:%s | MatScore:%s" % (gene,transcript,matGene,matScore))
            bcount+=1
            ## This trans has overlapping genes
            overGenes = ent[overlapGenes-1].replace('"','').split(",")
            overFlags = ent[overlapFlags-1].replace('"','').split(",")
            overConf = ent[overlapConf-1].replace('"','').split(",")
            overStatus = ent[overlapStatus-1].replace('"','').split(",")

            for x,y,z,w in zip(overGenes,overFlags,overConf,overStatus):
                # print(transcript,x,y,z,w)

                ## Check if entry matches with gene provided by match annot with good score
                # if x != matGene or (x==matGene and int(matScore) < 2):
                
                ## Avoid a transcript mapping to it's annotated counterpart, found by matchAnnot
                if x == matGene and int(matScore) > 2:
                    ## A final annotation score (not match annot score, gives 0 score to transcripts on opposite strand inlike matchannot)
                    print("Transcript matched to same genes as match Annot")
                    pass
                
                else:
                    ## Find overlap                    
                    overGene = x
                    overFlag = y
                    overConf = z
                    partStatus = w
                    # partinfo = statusDict[x]
                    # partStatus = partinfo[1]
                    
                    if z == 'R':
                        ## NAT config
                        print("NAT")
                        print(overGene,overFlag,partStatus)

                        if partStatus == "co" and transStatus != "co":
                            aclass = 0

                        elif partStatus != "co" and transStatus == "co":
                            aclass = 0

                        else:
                            aclass = 1

                        ccount+=1
                        fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene,transcript,overGene,transStatus,partStatus,overFlag,"nat",str(aclass)))

                    elif z == "F":
                        ## Same strand overlap
                        print("NS")
                        print(overGene,overFlag,partStatus)

                        if partStatus == "co" and transStatus != "co":
                            aclass = 0

                        elif partStatus != "co" and transStatus == "co":
                            aclass = 0

                        else:
                            aclass = 1

                        dcount+=1
                        fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene,transcript,overGene,transStatus,partStatus,overFlag,"ns",str(aclass)))
                    
                    else:
                        print("Check the entry for configuration")

        else:
            ## This transcript doesn't have any overlapping genes
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene,transcript,"NA",transStatus,"NA","NA","NA","NA"))

        
    fh_in.close()
    fh_out.close()

    print("\nTotal entries scanned: %s | Total overlaps:%s | NAT:%s | NS: %s" % (acount,bcount,ccount,dcount))
    print("Exiting function - summaryDict\n")


    return outFile

def summaryDict(summary):
    '''Create a dictionary of transcript and protein-coding status'''
    print("\nFunction: summaryDict")
    
    fh_in = open(summary,'r')
    if head == 'Y':
        fh_in.readline()

    summRead = fh_in.readlines()

    summDict = {} 
    acount = 0 ## To count the entries
    for i in summRead:
        ent = i.split("\t")  ## It has to be tab seprate dfile always
        transid =  ent[transpos-1]
        transStatus = ent[codStatus-1]

        if transStatus == "protein_coding" or transStatus == "coding" or transStatus == "Coding" or transStatus == "Cod" or transStatus == "NA" or transStatus == "co":
            transStatus2 = "co"
        elif transStatus == "non-coding" or transStatus == "Non-coding" or transStatus =="Non-Cod" or transStatus == "nc" or transStatus == "nc":
            transStatus2 = "nc"
        elif transStatus == "up":
            transStatus2 == "up"
        else:
            print("Unknown coding status '%s'" % transStatus)
            sys.exit()
        
        # print("Key:%s | Value:%s" % (transid,transStatus))
        summDict[transid] = ("pred",transStatus2)
        acount+=1

    print("Total entries scanned: %s | Length of GTF Dictionary %s" % (acount,len(summDict)))
    print("Exiting function - summaryDict\n")
    return summDict

def gtfDict(parsedGTF):
    ''' Create a dictionary of genes and their protein-coding status'''
    print("\nFunction: gtfDict")
    
    tempGene = "geneList.txt"
    fh_out = open(tempGene,'w')

    gtfDict = {}
    acount = 0 ## Count the entries 
    
    for i in parsedGTF:
        # print(i)
        achr,gStart,gEnd,gStrand,gid,gVer,gSource,gType = i
        
        if gType == "protein_coding" or gType == "coding" or gType == "Coding":
            gType2 = "cod"
        else:
            gType2 = gType
        
        val = (gSource,gType2)
        key = gid
        # print("Key:%s | Value:%s" % (key,val))
        gtfDict[key] = val
        fh_out.write("%s\t%s\n" % (key,val))
        acount+=1

    fh_out.close()

    print("Total entries scanned: %s | Length of GTF Dictionary %s" % (acount,len(gtfDict)))
    print("Exiting function - gtfDict\n")

    return gtfDict

def gtfParser(gtf):

    print("\nFunction: gtfParser")

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

def main():
    parsedGTF  = gtfParser(gtf)
    geneDict = gtfDict(parsedGTF)
    summDict = summaryDict(summary)
    NATres = NATdetect(summary,summDict,geneDict)

if __name__ == '__main__':
    main()
    sys.exit()


## Put them in classes
## Fix geneOverlap to find overlap within PacBio too and then use this script to get final results

## v01 ->v02
## Added functionality to deal with PacBio as overlapping transcript


