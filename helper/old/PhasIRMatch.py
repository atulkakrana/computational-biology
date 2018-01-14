#!/usr/local/bin/python3
#### Script to find Inverted repeats that enclave PHAS loci
#### Requires coordinates form phasing analysis and whole genome IR analysis 
######## IMPORT ################
import sys,os,subprocess,multiprocessing,datetime,time

################################

genome = 'ASPARAGUS_UGA1'   ## Too add to coordinates
phasFile = 'Final_24PHAS_Loci_ALL_v4.csv'
phasSep = ','
phasName = 2
phasChr = 5
phasStart = 6 
phasEnd = 7

IRfile = 'RES_Asparagus_edited.fa_03_31_16_58.csv'
IRsep = ','
IRchr = 1                   ## If no chr column then give first column as its basically chrxx
IRstart = 7
IRend = 10

def reader(phasFile,IRfile):
    ''' Prepares list of phased loci and 
    Inverted repeats from respective file'''

    print ("Preparing list of phased loci and inverted repeats\n")
    phasList = [] ## chr, start and end
    IRlist = [] ## IR chr, start and end

    print ("Creating list of phased loci")
    fh_in = open(phasFile,'r')
    for i in fh_in:
        ent = i.split(phasSep)
        print(ent[phasName-1],ent[phasChr-1],ent[phasStart-1],ent[phasEnd-1])
        phasList.append((ent[phasName-1],ent[phasChr-1],ent[phasStart-1],ent[phasEnd-1])) ## Coords converted to python format

    fh_in.close()

    print ("Creating list of inverted repeats")
    fh_in2 = open(IRfile,'r')
    for j in fh_in2:
        ent2 = j.split(IRsep)
        print(ent2[IRchr-1].replace("chr",""),ent2[IRstart-1],ent2[IRend-1])
        IRlist.append((ent2[IRchr-1].replace("chr",""),ent2[IRstart-1],ent2[IRend-1])) ## Coords converted to python format
    fh_in2.close()

    return phasList,IRlist

def Match(phasList,IRlist):
    '''Finds the inverted repats for phased loci'''

    print ("Searching the inverted repeats for phased loci")

    fileOut = 'MatchedPHAS.tsv'
    fh_out = open(fileOut, 'w')
    fh_out.write("PhasCoord\tIRcoord\tchr\tstart\tend\n")

    pcount = 0
    icount = 0
    acount = 0
    for phas in phasList:
        pName,pChr,pStart,pEnd = phas
        pcount +=1
        for inv in IRlist:
            iChr,iStart,iEnd = inv
            icount += 1

            if pChr == iChr: ## Same chromosome
                if pStart >= iStart: ## Phas Loci is witin start of inverted
                    if pEnd <= iEnd: ## Phas loci is before end of inverted
                        print("Phased locus:",phas,"is enclaved by inverted repeat:",inv,"\n")
                        iName = "%s.%s.%s.%s" % (genome,iChr,iStart,iEnd)
                        fh_out.write("%s\t%s\t%s\t%s\t%s\n" % (pName,iName,pChr,pStart,pEnd))
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



