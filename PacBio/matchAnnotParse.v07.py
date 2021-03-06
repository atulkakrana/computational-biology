#!/usr/local/bin/python3


import os,sys

## USER SETTINGS ###
afile           = "trinity.hybird.collapsed.matchedAnnot.agpv2.out" ## Output file from matchAnnot
collapsetype    = 0 ## 0: collapsed using genome | 1: collapsed using cd-hit 

def matchAnnoParse(afile):

    outFile = "%s.parsed.txt" % (afile)
    fh_out = open(outFile,'w')
    fh_out.write("PacBioGene\tPacBioTrans\tmapStart\tmapEnd\tstrand\tChr\tmapLen\tmatGene\tmatGeneStart\tmatGeneEnd\tnatGeneStrand\tstartDiff\tendDiff\tmatTrans\t matTransExons\ttransExon\tmatScore\tmatOrient\tannotScore\tannoStatus\n")

    fh_in = open(afile,'r')
    fileRead = fh_in.read().split("\n\n")

    transList = [] ## Temporary list of annotated transcripts that matched to PacBio transcripts
    acount = 0 ## All PacBio transcripts in file
    bcount = 0 ## All those with annoStatus == N, i.e. new isoforms
    ccount = 0 ## PacBio transcripts with no overlap with annotated ones
    for block in fileRead:
        acount+=1
        # print("##############Block:\n",block)
        info = block.split("\n")
        # print("################Info:\n",info)
        
        geneList    = [] ## Temporary list of annotated genes that matched pacBio transcripts
        transList   = [] ## Temporary list of annotated transcripts that matched to PacBio transcripts

        ### Check if there is an alingment ###
        ######################################
        isoforminfo = info[0]
        isoforment  =  isoforminfo.split()
        
        if len(isoforment) == 2:
            print("No alingments found - writing res and skipping") ## If PacBio transcipts are collapsed w/o genome, some transcripts from sam file (input of matchannot) might not have hits to genome
            annoStatus = 'N' ## The transcript didn't mapped to genome i.e. a novel transcript missing from genome
            if collapsetype     == 0:
                pacTrans,loc,trans  = isoforment[1].split("|") ## PB.1.1|chromosome1:4773-9701(-)|ac44850/f1p1/2043
                pacGene             = pacTrans.rpartition(".")[0]
                fh_out.write("%s\t%s\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\t%s\n" % (pacGene, pacTrans,annoStatus))
                continue
            elif collapsetype   == 1:
                pacTrans    = isoforment[1].strip() ## bc28888/f1p1/1446
                pacGene     = pacTrans
                fh_out.write("%s\t%s\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\tnoa\t%s\n" % (pacGene, pacTrans,annoStatus))
                continue
            else:
                print("Please input correct 'collapsetype' mode")
                continue

        else:
            # print("There is an alingment - Parse the results")
            pass

        ### Proceed for parsing ########
        ################################
        for i in info:
            if i.startswith("isoform"):
                ent =  i.split()
                # print("\n\n",ent)
                if collapsetype     == 0:
                    pacTrans,loc,trans  = ent[1].split("|") ## PB.1.1|chromosome1:4773-9701(-)|ac44850/f1p1/2043
                    pacGene             = pacTrans.rpartition(".")[0]
                elif collapsetype   == 1:
                    pacTrans            = ent[1].strip() ## bc28888/f1p1/1446
                    pacGene             = pacTrans
                else:
                    print("Please input correct 'collapsetype' mode")
                    pass

                print('\nThis is the transcript ID: %s' % (pacTrans))
                
                mapStart    = ent[2]
                mapEnd      = ent[3]
                achr        = ent[4]
                strand      = ent[5]
                transLen    = ent[6]

            elif i.startswith("cigar"):
                ent =  i.strip().split()
                # print(ent[1:])
                alen = len(ent[1:]) ## Remove the first entry 'cigar:'
                cigarExons = int((alen+1)/2) ## This is the number of exons in PacBio transcripts from CIGAR
                print('Exons from CIGAR:%s' % (cigarExons))

            ## Annotated gene details
            elif i.startswith("gene"):
                ent2 = i.split()
                # print(ent2)
                matGene = ent2[1]
                geneStart = ent2[2]
                startDiff = ent2[3]
                geneEnd = ent2[4]
                endDiff = ent2[5]
                geneStrand = ent2[6]
                geneList.append((matGene,geneStart,startDiff,geneEnd,endDiff,geneStrand))

            ## Annotated transcripts details
            elif i.startswith("tr:"):
                ent3 = i.split()
                # print("This is ent:",ent3)

                if ent3[1] != '(none)':
                    matTrans,matTransExons = ent3[1],ent3[5]
                    # print("This is matTrans %s and its Exons %s " % (matTrans,matTransExons))
                    transList.append((matTrans,matTransExons))
                else:
                    ## None of the overlapping transcripts had match > score 0, so none matched to PacBio transcript
                    ## Skip this none line, because best matched annotated trans with 0 score is already added in 
                    ## transList, and will be fetched later while preparing final results
                    # print(ent3)
                    # print('Loose end-1 - Atul')
                    # sys.exit()
                    pass

            ## Mapping Result details
            elif i.startswith("result"):
                ent4 = i.split()
                print(ent4)
                if ent4[2] != "no_genes_found":
                    matTransIndex   = ent4.index("ex:")-1 ## -1 because transcript is one index before "ex"
                    # matGeneBlock    = ent4[2:-5]
                    # print (matGeneBlock,len(matGeneBlock))
                    
                    ### Get matched gene names - bug-6/7/16
                    if matTransIndex > 3:
                        ## The gene name has space between full name like "ACT 1"
                        matGeneBlock    = ent4[2:matTransIndex]
                        matGeneFinal        = '_'.join(matGeneBlock)
                        expectedBlockLen    = len(matGeneBlock)+7
                    else:
                        matGeneFinal    = ent4[2]
                        matTransFinal   = ent4[3] ## index = 3
                        transExon       = ent4[5] ## index = 5
                        matScore        = ent4[7] ## index = 7
                        expectedBlockLen = 8 ## Normally there are 8 entries in list

                    print("Match Gene:%s | Match Trans:%s | Trans Exons:%s | Match Score:%s" % (matGeneFinal,matTransFinal,transExon,matScore) )
                    print("Exons in PacBio transcript from Results:%s" % (transExon))
                    if len(ent4) > expectedBlockLen:
                        ## Additional information exists i.e. UTR diff if score greater then 2 or match strand if on reverse
                        if ent4[8] == 'rev': ## 'rev' in end, transcript and matched gene on different strands, flag it Reverse
                            matOrient = 'R'
                            fiScore = 0
                        else:
                            ## This result has a score of 3 or above and ent4[8] consists of 5' and 3' results
                            matOrient = 'F'
                            fiScore =  matScore
                    else:
                        matOrient = 'F'
                        fiScore = matScore

                    ### Assign annotation status based on score and orientation
                    if int(matScore) <= 2 and matOrient == 'F' :
                        ## It's a novel isoform
                        annoStatus = 'N'
                        bcount+=1
                    elif matOrient == 'R' :
                        ## Annotated gene on other strand, definately PacBio transcript is novel gene
                        annoStatus = 'N'
                        ccount+=1
                    else:
                        annoStatus = 'Y'
                
                elif ent4[2] == "no_genes_found": ###############
                    matGeneFinal    = "None"
                    annoStatus      = "N"
                    ccount+=1
                    # matTransFinal = "None"
                    transExon = cigarExons ## Tested OK
                    # matScore = "None"
                    # print("No match found - exons in PacBio transcripts:%s" % (transExon))
                    # sys.exit()

                else:
                    print("This point can't be reached normally - You screwed up")
                    sys.exit()
            
            else:
                ## Trash entry 
                pass


        if matGeneFinal != "None":
            ## Get exons of best matched transcripts, whose name is mentioned in result
            for i in transList:
                # print(i)
                matTrans1,matTransExons1 = i
                if matTrans1 == matTransFinal:
                    # print("Matched")
                    matFinalTransExons = matTransExons1
                    continue
                else:
                    # print("Not Matched")
                    pass
            
            ## Get gene that best matched to PacBio transcript, whose name is mentioned in result line
            for i in geneList:
                # print(i)
                matGene1 = i[0]
                if matGene1 == matGeneFinal:
                    # print("Matched")
                    geneStartFi = i[1]
                    startDiffFi = i[2]
                    geneEndFi = i[3]
                    endDiffFi = i[4]
                    geneStrandFi = i[5]
                    continue

                else:
                    # print("Not Matched")
                    pass

            # print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pacGene, pacTrans, mapStart, mapEnd, strand, achr, transLen, matGeneFinal, geneStartFi, geneEndFi, geneStrandFi, startDiffFi, endDiffFi, matTransFinal, matFinalTransExons, transExon, matScore))
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pacGene, pacTrans, mapStart, mapEnd, strand, achr, transLen, matGeneFinal, geneStartFi, geneEndFi, geneStrandFi, startDiffFi, endDiffFi, matTransFinal, matFinalTransExons, transExon, matScore, matOrient, fiScore, annoStatus))

        else:
            # print("%s\t%s\t%s\t%s\t%s\t%s\t%s\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t%s\tNone\n" % (pacGene, pacTrans, mapStart, mapEnd, strand, achr, transLen, transExon))
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%s\tNA\tNA\tNA\t%s\n" % (pacGene, pacTrans, mapStart, mapEnd, strand, achr, transLen, transExon, annoStatus))

    print("Total transcripts in file: %s | Putative new isoforms:%s | New Loci:%s" % (acount,bcount,ccount))



    fh_out.close()
    fh_in.close()

    return outFile

def main():
    outFile = matchAnnoParse(afile)


if __name__ == '__main__':
    main()
    sys.exit()


## v01 -> v02
## Fixed a minor error in provifing matched transcripts exons

## v02 -> v03
## Fixed header and variable names for chromosome and transLen

## v02 -> v03
## Catched transcriprts match from different strands and flagged them with annotation score 0, annoScore is the final score

## v04 -> v05
## Fixed wrong exons number issue for genes with no match, now their exons are fetched from CIGAR line
## match on other strand now counted as new gene

## v05 -> v06
## Added functionality to work with cd-hit collapsed transcripts

## v06 -> v07
## Fixed a bug related to gene name in "result" line when it's seprated by a space like "ACT 1". In such case all the list indexes become different. Such cases are caught, and expected number of entries computed. If list if longer then expected entries then presence of extra info is detcted 

