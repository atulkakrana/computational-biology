#!/usr/local/bin/python3
## Scipt for processing miRNA annotation results from BLAST to give final results

import os,sys

# blastTab =  str(sys.argv[1])
blastTab = "BLASTAsparagusNewLibsRelaxed21_22notInFInalv1.0.txt"

def process(blastTab):

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
        ## 5' overhang
        if qstart > sstart:
            hang5 = qstart-sstart

        ## 3' overhang
        elif send > qend: 
            # hang3 = send - qend
            if send == slen:  ## Reached end of subject, everything is 3' overhang
                hang3 = qlen-qend
            elif slen > send: ## More subject remaining but not alinged
                hang3 = qlen - ((slen-send)+qend)
                
                if hang3 < 0: ## Means subject was longer then query - set to zero
                    print("caught")
                    hang3 = 0
                else:
                    pass
            else:
                ## Not sure, if any other condition exits
                pass

        ## No overhang    
        else:
            pass

        print ("align:",align)
        match = round(align*(pid/100),0)
        mismatch = align-int(match)
        unalign = (qlen-align)-hang3-hang5 ## This included sequnce that was not aligned, 3' and 5' overhang
        totalMis = unalign+mismatch

        ## Test of fitness
        if hang3 <= 2 and hang5 ==0:
            if hang3+totalMis <= 4:
                status = "pass"
            else:
                status ='fail'
        elif hang5 <= 2 and hang3 == 0:
            if hang5 + totalMis <= 4:
                status = "pass"
            else:
                status ='fail'
        elif hang3 == 0 and hang5 == 0 and totalMis <= 4:
            status = "pass"
        
        else:
            status ='fail'

        ent.extend((str(hang5),str(hang3),str(match),str(mismatch),str(unalign),str(unalign+mismatch+hang3+hang5),status))
        resList.append((ent))

        # print (resList)

    fh_in.close()
    return resList


def writer(resList,blastTab):

    print ("\nWriting results to file")

    outFile = ("Final_%s" % (blastTab))
    fh_out = open(outFile,'w')
    fh_out.write("query\tsseqid\tpident\tlength\tqlen\tqstart\tqend\tslen\tsstart\tsend\tgaps\tmismatch\tpositive\tevalue\tbitscore\thang5\thang3\tmatch\tmismatch\tunalign\ttotalMis\tstatus\n")

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
    print("Script finished succesfully")
    sys.exit()

## Companion script of miRanno.sh