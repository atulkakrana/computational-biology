#!/usr/local/bin/python3


import os,sys

## USER SETTINGS ###
afile = "matchedAnnot.agpv2.out" ## Output file from matchAnnot

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
        # print(block)
        info = block.split("\n")
        # print(info)
        
        geneList = [] ## Temporary list of annotated genes that matched pacBio transcripts
        transList = [] ## Temporary list of annotated transcripts that matched to PacBio transcripts
        for i in info:
            if i.startswith("isoform"):
                ent =  i.split()
                # print(ent)
                pacTrans,loc,trans = ent[1].split("|")
                pacGene = pacTrans.rpartition(".")[0]
                mapStart = ent[2]
                mapEnd = ent[3]
                achr = ent[4]
                strand = ent[5]
                transLen = ent[6]

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

            elif i.startswith("tr:"):
                ent3 = i.split()
                # print("This is ent:",ent3)

                if ent3[1] != '(none)':
                    matTrans,matTransExons = ent3[1],ent3[5]
                    # print("This is matTrans %s and its Exons %s " % (matTrans,matTransExons))
                    transList.append((matTrans,matTransExons))


            elif i.startswith("result"):
                ent4 = i.split()
                # print(ent4)
                if ent4[2] != "no_genes_found":
                    matGeneFinal = ent4[2]
                    matTransFinal = ent4[3]
                    transExon = ent4[5]
                    matScore = ent4[7]
                    if len(ent4) > 8:
                        ## Additional information exits i.e. UTR diff if score greater then 2 or match strand if on reverse
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

                    if int(matScore) <= 2 or matOrient == 'R' :
                        ## See matchAnno Scoring - Score <= 2 suggest that exons did not mapped with overlapping gene and it's a novel isoform
                        annoStatus = 'N'
                        bcount+=1
                    else:
                        annoStatus = 'Y'
                else:
                    matGeneFinal = "None"
                    annoStatus = "N"
                    ccount+=1
                    # matTransFinal = "None"
                    # transExon = "None"
                    # matScore = "None"
            else:
                pass


        if matGeneFinal != "None":

            ## Get exons of best matched transcripts, whose name is mentioned in result
            for i in transList:
                # print(i)
                matTrans1,matTransExons1 = i
                if matTrans1 == matTransFinal:
                    print("Matched")
                    matFinalTransExons = matTransExons1
                    continue

                else:
                    print("Not Matched")
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
                    print

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
## Catched transcriprts match from different starands and flagged them with annotation score 0, annoScore is the final score