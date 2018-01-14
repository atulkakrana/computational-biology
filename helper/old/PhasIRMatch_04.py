#!/usr/local/bin/python3
#### Script to find Inverted repeats that enclave PHAS loci
#### Requires coordinates form phasing analysis and whole genome IR analysis 
######## IMPORT ################
import sys,os,subprocess,multiprocessing,datetime,time

################################

genome      = 'RICE_MSU7'   ## For coordinate ID
phasFile    = 'Final_21PHAS_Loci_ALL_v6.txt'
head1       = 'T'
phasSep     = '\t'
phasName    = 4
phasChr     = 9
phasStart   = 10 
phasEnd     = 11

IRfile      = 'RES_Asparagus.V2.0.genome.stable_modifed.fa_10_11_13_33.csv'
head2       = 'T'
IRsep       = ','
IRchr       = 1              ## If no chr column then give first column as its basically chrxx
IRstart     = 7
IRend       = 10
IRloop      = 11
IRgap       = 5
IRscore     = 2
IRMatch     = 3

def reader(phasFile,IRfile):
    ''' Prepares list of phased loci and 
    Inverted repeats from respective file'''

    print ("Preparing list of phased loci and inverted repeats\n")
    phasList = [] ## chr, start and end
    IRlist = [] ## IR chr, start and end

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
        IRlist.append((iChr,ent2[IRscore-1],ent2[IRMatch-1],ent2[IRgap-1],ent2[IRstart-1],ent2[IRend-1],ent2[IRloop-1].strip("\n"))) ## Coords converted to python format
    fh_in2.close()

    return phasList,IRlist

def Match(phasList,IRlist):
    '''Finds the inverted repats for phased loci'''

    print ("Searching the inverted repeats for phased loci")

    fileOut = 'MatchedPHAS.tsv'
    fh_out = open(fileOut, 'w')
    fh_out.write("PhasCoord\tIRcoord\tChr\tScore\tMatch\tGap\tStart\tEnd\tLoop\n")

    pcount = 0
    icount = 0
    acount = 0
    for phas in phasList:
        pName,pChr,pStart,pEnd = phas
        print("\n",pName,pChr,pStart,pEnd)

        pcount +=1
        for inv in IRlist:
            # print(inv)
            iChr,iScore,iMatch,iGap,iStart,iEnd,iLoop = inv
            # print(iChr,iScore,iMatch,iGap,iStart,iEnd,iLoop)
            icount += 1

            if str(pChr) == str(iChr): ## Same chromosome - don't use int here because some chromosomes are from mitochondria and all
                if int(pStart) >= int(iStart): ## Phas Loci is witin start of inverted
                    if int(pEnd) <= int(iEnd): ## Phas loci is before end of inverted
                        print("Phased locus:",phas,"is enclaved by inverted repeat:",inv,"\n")
                        iName = "%s.%s.%s.%s" % (genome,iChr,iStart,iEnd)
                        fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (pName,iName,pChr,iScore,iMatch,iGap,iStart,iEnd,iLoop))
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



