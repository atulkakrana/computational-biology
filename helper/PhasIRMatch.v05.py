#!/usr/local/bin/python3
#### Script to find Inverted repeats that enclave PHAS loci
#### Requires coordinates form phasing analysis and whole genome IR analysis 
######## IMPORT ################
import sys,os,subprocess,multiprocessing,datetime,time

################################

genome      = 'MAIZE_AGPv2'   ## For coordinate ID
phasFile    = 'Final_PHAS_Loci_1e-05_24ALL.csv'
head1       = 'T'
phasSep     = ','
phasName    = 1
phasChr     = 3
phasStart   = 4 
phasEnd     = 5

IRfile      = 'RES_all.fa_11_28_20_23.repLen5k.csv'
head2       = 'T'
IRsep       = ','
IRchr       = 1              ## If no chr column then give first column as its basically chrxx
IR5start    = 7
IR5end      = 8
IR3start    = 9
IR3end      = 10
IRloop      = 11
IRgap       = 5
IRscore     = 2
IRMatch     = 3

def reader(phasFile,IRfile):
    ''' Prepares list of phased loci and 
    Inverted repeats from respective file'''

    print ("Preparing list of phased loci and inverted repeats\n")
    phasList    = [] ## chr, start and end
    IRlist      = [] ## IR chr, start and end

    print ("Creating list of phased loci")
    fh_in = open(phasFile,'r')
    if head1 == 'T':
        fh_in.readline()
    

    for i in fh_in:
        ent = i.split(phasSep)
        # print(ent[phasName-1],ent[phasChr-1],ent[phasStart-1],ent[phasEnd-1])
        phasList.append((ent[phasName-1],ent[phasChr-1],ent[phasStart-1],ent[phasEnd-1])) ## Coords converted to python format

    fh_in.close()

    print ("Creating list of inverted repeats")
    fh_in2 = open(IRfile,'r')
    if head2 == 'T':
        fh_in2.readline()

    for j in fh_in2:
        ent2 = j.split(IRsep)
        iChr = ent2[IRchr-1].replace("Chr","").replace("chr","")
        # print(ent2[IRchr-1].replace("chr",""),ent[IRgap-1],ent2[IRscore-1],ent2[IRMatch-1],ent2[IRstart-1],ent2[IRend-1],ent[IRloop-1]) ## Chr, score, match,Gap, Start, End , Loop
        IRlist.append((iChr,ent2[IRscore-1],ent2[IRMatch-1],ent2[IRgap-1],ent2[IR5start-1],ent2[IR5end-1],ent2[IR3start-1],ent2[IR3end-1],ent2[IRloop-1].strip("\n"))) ## Coords converted to python format
    fh_in2.close()

    return phasList,IRlist

def Match(phasList,IRlist):
    '''Finds the inverted repats for phased loci'''

    print ("Searching the inverted repeats for phased loci")

    fileOut = 'MatchedPHAS.tsv'
    fh_out = open(fileOut, 'w')
    fh_out.write("PhasCoord\tIRcoord\tChr\tScore\tMatch\tGap\tStart5\tEnd5\tStart3\tEnd3\tLoop\tPHAS\tphas.start\tphas.end\t5flank\t3flank\n")

    pcount = 0
    icount = 0
    acount = 0
    for phas in phasList:
        pName,pChr,pStart,pEnd = phas
        pCluster = "%s.%s.%s.%s" % (genome,pChr,pStart,pEnd)
        print("\n",pName,pChr,pStart,pEnd)

        pcount +=1
        for inv in IRlist:
            # print(inv)
            iChr,iScore,iMatch,iGap,i5Start,i5End,i3Start,i3End,iLoop = inv
            # print(iChr,iScore,iMatch,iGap,iStart,iEnd,iLoop)
            icount += 1

            if str(pChr) == str(iChr): ## Same chromosome - don't use int here because some chromosomes are from mitochondria and all
                if int(pStart) >= int(i5Start): ## Phas Loci is witin start of inverted
                    if int(pEnd) <= int(i3End): ## Phas loci is before end of inverted
                        print("Phased locus:",phas,"is enclaved by inverted repeat:",inv,"\n")
                        iName = "%s.%s.%s.%s" % (genome,iChr,i5Start,i3End)

                        ## Compute overlap with 5'ARM
                        if int(pEnd) <= int(i5End):
                            ## Phas is contained in 3 ARM
                            flank5 = int(pEnd)-int(pStart)
                        else:
                            flank5 = int(i5End)-int(pStart)
                        
                        ## Compute overlap with 3' ARM
                        if int(pStart) >= int(i3Start):
                            flank3 = int(pEnd)-int(pStart)
                        else:
                            flank3 = int(pEnd)-int(i3Start)

                        fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pName,iName,pChr,iScore,iMatch,iGap,i5Start,i5End,i3Start,i3End,iLoop,pCluster,pStart,pEnd,flank5,flank3))
                        acount+=1 
    fh_out.close()
    print("\nPhased Locus:%s | IRs: %s | Matched: %s" % (pcount,icount,acount))
    return fileOut 


def main():

    phasList,IRlist = reader(phasFile,IRfile)
    fileOut = Match(phasList,IRlist)


    pass


if __name__ == '__main__':
    main()
    sys.exit()


## v0.2 -> 0.3
## Fixed bug when genomic file is used for IR prediction and has "Chr"/"chr" in chromosome
## Also mitrochondrial chromosome will pass w/o any eror

## v04 -> v05
## Added 5end and 3start to results

