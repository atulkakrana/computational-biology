#!/usr/local/bin/python3
## Scipt for processing miRNA annotation results from BLAST to give final results

import os,sys

blastTab =  str(sys.argv[1])
# blastTab = "BLASTAsparagusNewLibsRelaxed21_22notInFInalv1.0.txt"

def process(blastTab):

    '''process the tabular results from BLAST to compute overhangs, mismatches, matches and ttotalVariance.
    The unalign includes - overhangs because these are the regions which did not aligned
    '''

    fh_in = open(blastTab,'r')
    fileRead = fh_in.readlines()


    resList = [] ## To hold results
    for i in fileRead:
        ent = i.strip('\n').split('\t')
        print("query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore")
        print(ent)

        ## Format - query sseqid pident length qlen qstart qend slen sstart send gaps mismatch positive evalue bitscore
        qlen = int(ent[4])
        qstart = int(ent[5])
        qend = int(ent[6])
        slen = int(ent[7])
        sstart = int(ent[8])
        send = int(ent[9])
        align = int(ent[3])
        pid =  int(round(float(ent[2]),0))

        hang5 = 0
        hang3 = 0
        
        # ## 5' overhang - Tested OK --- OLD - Deprectaed to compute overhang in subject too
        # if qstart > sstart:
        #     hang5 = qstart-sstart
        # else:
        #     # print("No 5' overhand")
        #     pass

        ## 5' overhang
        qReminder5 = 1 - qstart ## Remaining query nucleotides at 5' before alignment starts
        sReminder5 =  1 - sstart ## Remaining subject nucleotides at 5' before alignment starts
        hang5 = abs(qReminder5) - abs(sReminder5) ## If negative that means subject had overhang
        
        ## 3' overhang
        qReminder3 = qlen - qend ## Remaining query ucleotides after alignment ends
        sReminder3 =  slen - send ## Remaining subject nucleotides after alignment ends
        hang3 = qReminder3 - sReminder3 ## If negative that means subject had overhang

        ## Compute mismatch in aligned, unaligned and total variance
        print ("align:",align)
        match = round(align*(pid/100),0)
        mismatch = align-int(match)
        unalign = (qlen-align) ## This included sequnce that was not aligned, including 3' and 5' overhang
        totalVar = unalign+mismatch

        ## Test of fitness
        
        ## Catch overhang in 5' under cutoff, check for positive overhang i.e. in query
        if hang5 <= 2 and hang3 == 0:
            if totalVar <= 5:
                status = "pass"
            else:
                status ='fail'
        
        ## atch overhang in 3' under cutoff, check for positive overhang i.e. in query
        elif hang3 <= 2 and hang5 == 0:
            if totalVar <= 5:
                status = "pass"
            else:
                status ='fail'
        
        ## If query is long i.e. 24 nt and subject is small 21nt - 1-3nt overhang expected
        ## so overhangs + 3 mismatches or unalign is cutoff 6
        elif hang5 >= 1 and hang3 >= 1: ## There is overhang at both sides
            if totalVar <= 6:
                status = "pass"
            else:
                status ='fail'
        
        ## No overhang - perfect situation
        elif hang3 == 0 and hang5 == 0 and totalVar <= 4:
            status = "pass"
        
        ## Cases with more then 2nt overhang and mismatches+unalign - fail
        else:
            status ='fail'

        ent.extend((str(hang5),str(hang3),str(match),str(mismatch),str(unalign),str(totalVar),status))
        resList.append((ent))

        # print (resList)

    fh_in.close()
    return resList

def writer(resList,blastTab):

    print ("\nWriting results to file")

    outFile = ("Final_%s" % (blastTab))
    fh_out = open(outFile,'w')
    fh_out.write("query\tsseqid\tpident\tlength\tqlen\tqstart\tqend\tslen\tsstart\tsend\tgaps\tmismatch\tpositive\tevalue\tbitscore\thang5\thang3\tmatch\tmismatch\tunalign\ttotalVariance\tstatus\n")

    for res in resList:
        print(res)
        fh_out.write("%s\n" % (str('\t'.join(res))) )

    fh_out.close()
    return outFile

def main():
    resList = process(blastTab)
    writer(resList,blastTab)

if __name__ == '__main__':
    main()
    print("\nScript finished succesfully")
    sys.exit()

## Companion script of miRanno.sh
## v1 -> v2
## Overhangs computed in both reference and subject are reflected as positive and negative integers
    ## Overhangs in qury i.e. positive integers aare the one those matters