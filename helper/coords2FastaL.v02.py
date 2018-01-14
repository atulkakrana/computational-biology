#!/usr/local/bin/python3


### This script extracts seqeunce from local genome file, given a ta-sperated or csv file with co-ordinates

import os,sys,time,sqlite3,operator
from operator import itemgetter

genomeFasta     = "/home/kakrana/99.genomes/maize.agp.v2/Zea_mays.AGPv2.17.dna.allchromosome"

coordsFile      = "FL.target-based.longprm.txt"
delim           = "\t"                                                      ## Delimiter for summary file
head            = 'Y'                                                       ## Header is summary file: 'Y' else: 'N'
namepos         = 1
chrpos          = 3
strandpos       = 
startpos        = 4
endpos          = 5


def cacheGenome(genomeFasta):
    '''
    Read genome file and retruns a dict
    '''
    print("\nFunction: cacheGenome")

    genomeD     = {} ## A dictionary to store fasta files
    fh_in       = open(genomeFasta)
    fasta       = fh_in.read()
    fasta_splt  = fasta.split('>')
    acount      = 0 ## count the number of entries
    for i in fasta_splt[1:]:
        ent     = i.split('\n')
        name    = ent[0].split()[0].strip()
        seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        alen    = len(seq)
        genomeD[name] = seq
        acount+=1
        print("seqeunce %s of length %s cached" % (name,alen))

    print("Sequences cached:%s" % (acount))

    print ("Exiting function - cacheGenome\n")

    return genomeD

def parseCoords(coordsFile):
    '''
    Reads coords file
    '''
    print("\nFunction: parseCoords")
    
    fh_in = open(coordsFile,'r')
    if head == 'Y':
        fh_in.readline()
    coordsRead = fh_in.readlines()
    fh_in.close()

    coordsL = [] ## Store all coords
    for i in coordsRead:
        ent     = i.split(delim)
        aname   = ent[namepos-1].strip()
        achr    = int(ent[chrpos-1].replace("chr",""))
        astart  = int(ent[startpos-1])
        aend    = int(ent[endpos-1])
        astrand = ent[strandpos-1]

        print("Phas Name %s | chr:%s | start:%s | end:%s | astrand:%s" % (aname,achr,astart,aend,astrand))
        coordsL.append((aname,achr,astart,aend,astrand))

    print("Coords list prepared:%s" % len(coordsL))
    print("Exiting function - parseCoords\n")
    
    return coordsL

def fetchSequences(genomeD,coordsL):
    '''
    Fetches and writes seqeunces
    '''
    print("\nFunction: fetchSequences")
    
    fastaFile   = "%s.fa" % coordsFile.rpartition(".")[0]
    fh_out      = open(fastaFile,'w')

    ## Sorts coords on chr to speed up
    coordsL_sorted = sorted(coordsL, key=itemgetter(1))
    for i in coordsL_sorted:
        aname,achr,astart,aend,astrand = i
        print(aname,achr,astart,aend,astrand)
        chrseq = genomeD[str(achr)]

        aseq = chrseq[astart:aend+1]
        if astrand == "w":
            fh_out.write(">%s\n%s\n" % (aname,aseq))
        elif astrand == "c":
            fh_out.write(">%s\n%s\n" % (aname,aseq[::-1].translate(str.maketrans("tagcTAGC","atcgATCG"))))
        else:
            print("Unknown strand encountered")
            sys.exit()

    print("Exiting function - fetchSequences\n")

    return fastaFile

def main():
    genomeD     = cacheGenome(genomeFasta)
    coordsL     = parseCoords(coordsFile)
    fastaFile   = fetchSequences(genomeD,coordsL)


if __name__ == "__main__":
    main()
    print("Script finished sucessfully")
    sys.exit()