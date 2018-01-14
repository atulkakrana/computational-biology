#!/usr/local/bin/python3
## Written by Atul Kakrana - Kakrana@udel.edu

## This script requires summary file from PACBIO as generated including results matchAnnot, CPC and CPAT
## This will perform final classification of ncRNA and codingRNAs
## IMPORTS
import os,sys

## USER SETTINGS ##################################################################
summary         = "summary.agp2.v4.txt"       ## Should be tab seprated file, because overlap columns have comma seprated values  
sep             = '\t'                      
head            = 'T'                 ## T: Header present in summary file | F: No header present in sumamry file
genepos         = 1                   ## PacBio gene ID, generated after collapse - PB.1 i.e. without transcript part
transpos        = 2                   ## Coulmn number for transcript ID - PB.1.1, not for gene ID                
matScorePos     = 19
orfPos          = 9                   ## Length of ORF, from CPC analysis or CPAT analysis
cpatScorePos    = 26
cpatThres       = 0.22                ## Determined by ROCR analysis
cpcScorePos     = 28                  ## Using PLEK score for testing
ORFThres        = 75 * 3              ## Amino acid * 3
 
def classify(summary,cpatThres):

    ''' Classify PACBIO transcripts to coding and non-coding based on -
    matchAnnot scores, CPC score, longest ORF and CPAT score'''

    print("\nFunction:Classify")

    fh_in = open(summary,'r')
    if head == 'T':
        fh_in.readline()
    summRead = fh_in.readlines()

    acount = 0         ## Count coding
    bcount = 0         ## Count non-coding
    ccount = 0         ## Count TUCP
    dcount = 0         ## Count entries
    ecount = 0         ## Undetected

    resList = [] ## Stores final result
    for i in summRead:
        ent = i.strip('\n').split(sep)
        print(ent[1:5])
        gene = ent[genepos-1]
        trans = ent[transpos-1]
        matScore = (ent[matScorePos-1])
        cpcScore = float(ent[cpcScorePos-1])
        orf = int(ent[orfPos-1])
        cpatScore = float(ent[cpatScorePos-1])
        print(gene,trans,matScore,cpcScore,orf,cpcScore)
        dcount+=1

        ## CLASSIFY ##############

        ## Transcripts with no matching gene or isoform
        if matScore == 'NA':

            if cpcScore <= (-1):
                ## ncRNA - No further testing required
                status = 'nc'
                resList.append((trans,status))
                bcount+=1

            elif cpcScore >= (1):
                ## codRNA - ncRNA overlapping coding genes might hide here - Test on ORF
                if orf >= ORFThres:
                    status = 'co'
                    resList.append((trans,status))
                    acount+=1

                else:
                    ## orf less then ORFThres:
                    ## Check CPAT status
                    if cpatScore > cpatThres:
                        ## Both cpc,cpat says its coding
                        status = 'co'
                        resList.append((trans,status))
                        acount+=1
                    else:
                        status = 'up'
                        resList.append((trans,status))
                        ccount+=1

            else:
                ## Weak evidence - Resolve with CPAT
                # if orf < ORFThres and cpcScore  <= 0 and cpatScore < cpatThres:
                if cpcScore  <= 0 and cpatScore < cpatThres:
                    status = 'nc'
                    resList.append((trans,status))
                    bcount+=1
                elif cpcScore  >= 0 and cpatScore > cpatThres:
                    status = 'co'
                    resList.append((trans,status))
                    acount+=1
                else:
                    ## Unknown potential i.e. TUCP
                    status = 'up'
                    resList.append((trans,status))
                    ccount+=1

        ## Transcripts with good match score
        elif int(matScore) >= 3:
            ## Decide on basis of CPC and CPAT score - This has been changed in version-3
            
            if cpcScore >= 1:
                status = 'co'
                resList.append((trans,status))
                acount+=1
            elif cpcScore <= -1:
                status = 'nc'
                resList.append((trans,status))
                bcount+=1
            else:
                ## Weak Evidence - Test with CPAT
                if cpcScore > 0 and cpatScore > cpatThres:
                    ## Both says its coding
                    status = 'co'
                    resList.append((trans,status))
                    acount+=1
                elif cpcScore < 0 and cpatScore < cpatThres:
                    ## Both says its nc
                    status = 'nc'
                    resList.append((trans,status))
                    bcount+=1
                else:
                    status = 'up'
                    resList.append((trans,status))
                    ccount+=1

        ## Transcripts with lower match annot score i.e. new isoforms - This may have ncRNAs overlapping protein-coding transcripts
        elif int(matScore) <= 2 and int(matScore) >= 0:
            ## These are new isoforms - equal probability of ncRNA and codRNA

            if cpcScore >= 1:
                ## Coding unless short orfs, ncRNA overlapping coding gene
                if orf >= ORFThres:
                    status = 'co'
                    resList.append((trans,status))
                    acount+=1
                else:
                    ## Short ORF - Check with CPAT
                    if cpatScore > cpatThres:
                        ## Both cpc,cpat says its coding
                        status = 'co'
                        resList.append((trans,status))
                        acount+=1
                    else:
                        ## Short ORF and conflicting cpc,cpat scores
                        status = 'up'
                        resList.append((trans,status))
                        ccount+=1

            elif cpcScore <= -1:
                ## Non-coding 
                status = 'nc'
                resList.append((trans,status))
                bcount+=1

            else:
                ## Weak evidence - Test with cpc,orf, and CPAT
                if cpcScore <= 0 and cpatScore < cpatThres:
                    status = 'nc'
                    resList.append((trans,status))
                    bcount+=1
                elif cpcScore > 0 and cpatScore > cpatThres:
                    status = 'co'
                    resList.append((trans,status))
                    acount+=1
                else:
                    status ='up'
                    resList.append((trans,status))
                    ccount+=1

        else:
            print("Please check the matchAnnot score for entry:%s" % (matScore))
            ecount+=1

    fh_in.close()
    print("5 entries from resList:",resList[1:5])
    print("Total entries in file:%s | Coding:%s | Non-coding:%s | TUCP:%s | undetected: %s" % (dcount,acount,bcount,ccount,ecount))
    print("Exiting - Classify\n")
    return resList

def writer(resList):
    '''Output tab seprated file'''

    print("\nFunction: Writer")

    outFile = "%s.ncClass.txt" % (summary)
    fh_out = open(outFile,'w')

    acount = 0 ## Entries written
    for i in resList:
        trans,aclass = i
        fh_out.write("%s\t%s\n"% (trans,aclass))
        acount +=1

    print("Entries in list %s | entries written:%s" % (len(resList),acount))

    fh_out.close()
    print("Exiting - writer")
    return outFile


def main():
    resList=classify(summary,cpatThres)
    resFile = writer(resList)

if __name__ == '__main__':
    main()


## v.01 [Jul29]

## Changes:
## 1. DO we need to use ORF filter in weak evidences?

## v01 - > v03
## Dependency on match annot score, if greater then 3 is removed