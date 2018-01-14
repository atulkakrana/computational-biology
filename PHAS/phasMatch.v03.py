#!/usr/local/bin/python3

### This script takes phaser summary file and phastank pred_tab file
### and matches those with each other to give a summary of overlap

import sys,os,time,datetime
from operator import itemgetter

phaster  = str(sys.argv[1])
phastank = str(sys.argv[2])
revferno = str(sys.argv[3])
flanksize= 150              ## Flaking region used to match at head and tail

def revferno_parse(referno):
    '''
    parses validated revferno results, chosses best trigger and provides a dict for phaster_parse
    '''    

    print("\nFn: revferno_parse ####################")
    revfernoL = []
    revfernoD = {}

    fh_in   = open(referno,'r')
    fh_in.readline()
    afile   = fh_in.readlines()
    fh_in.close()

    acount  = 0 ### Counts all interactions in triggers file
    bcount  = 0 ### Counts unique PHAS
    trigsL  = []
    tempset = set() ## TO keep track of unique PHAS
    for i in afile:
        # print("Ent:",i)
        ent     = i.strip('\n').split(',')
        amir    = ent[0].split('-')[1]
        aphas   = ent[1].split('_')[0]
        aloci   = ent[20]
        aindex  = int(ent[21])
        aflag   = ent[22]
        tempset.add(aphas)
        acount  +=1
        # print("Parsed:",amir,aphas,aloci,aindex,aflag)
        if abs(aindex) <= 3:
            # print("Included index",abs(aindex))
            trigsL.append((amir,aphas,aloci,aindex,aflag))
        else:
            # print("Excluded index:",abs(aindex))
            pass
    
    # print("-Total PHAS:%s | Total interactions:%s | Total included:%s" % (len(tempset),acount,len(trigsL)))

    ## Sort the list
    trigsL_s = sorted(trigsL,key=itemgetter(3),reverse = False)
    # print("-Snippet of sorted triggers",trigsL_s[:5])

    phasset     = set() ## Record PHAS that have been processed 
    ### Seggregate results and select best candidates
    for aent in trigsL_s:
        amir,aphas,aloci,aindex,aflag = aent

        if aphas not in phasset:
            if amir.startswith("miR2275") or amir.startswith("miR2118") or amir.startswith("miR482"):
                revfernoL.append(aent)
                revfernoD[aphas] = aent
                phasset.add(aphas)
        else:
            ## Trigger for PHAS already added
            pass

    ### Unique other triggers, add those w/o any ealrier identified trigger
    for aent in trigsL_s:
        amir,aphas,aloci,aindex,aflag = aent

        if aphas not in phasset:
            if not amir.startswith("miR2275") or not amir.startswith("miR2118") or not amir.startswith("miR482"):
                revfernoL.append(aent)
                revfernoD[aphas] = aent
                phasset.add(aphas)
        else:
            ## Trigger ofor PHAS already added
            pass

    print("Trigs in file:%s| Trigs cached:%s" % (acount,len(revfernoL)))

    return revfernoL,revfernoD

def phaster_parse(phaster,revfernoD):
    '''
    parses phaster file
    '''
    print("\nFn: phaster_parse ####################")
    phasterL = []
    phasterD = {}

    fh_in   = open(phaster,'r')
    fh_in.readline()
    afile   = fh_in.readlines()
    fh_in.close()

    acount = 0
    for i in afile:
        ent     = i.strip('\n').split('\t')
        aname   = ent[0]
        apval   = ent[1]
        achr    = ent[2]
        astart  = ent[3]
        aend    = ent[4]
        xtrig   = revfernoD.get(aname,"na")
        # print(xtrig)
        if xtrig    == "na":
            atrig   = "na"
            aindex  = "na"
            aloci   = "na"
        else:
            atrig   = xtrig[0]
            aindex  = xtrig[3]
            aloci   = xtrig[2]

        avalue  = (aname,achr,astart,aend,apval,atrig,aindex,aloci)
        # print("Phaster parsed:",avalue)
        phasterL.append(avalue)
        phasterD[aname] = avalue
        acount          +=1

    print("PHAS in file:%s| PHAS cached:%s" % (acount,len(phasterL)))

    return phasterL,phasterD

def phastank_parse(phastank):
    '''
    parses phasetank file
    '''

    print("\nFn: phastank_parse ###################")
    phastankL = []
    phastankD = {}

    fh_in   = open(phastank,'r')
    fh_in.readline()
    afile   = fh_in.readlines()
    fh_in.close()

    acount = 0
    for i in afile:
        # print(i)
        ent     = i.strip('\n').split('\t')
        aname   = ent[1]
        ascore  = ent[7]
        achr    = aname.rpartition("_")[0].replace("chromosome_","").replace("chr","") ## Need better handling
        astart,aend  = ent[2].strip().split(":")
        xtrig   = ent[8]
        if xtrig != "NONE":
            atrig = xtrig.split('-')[1]
        else:
            atrig = xtrig

        avalue  = (aname,achr,astart,aend,ascore,atrig)
        # print("PhaseTank parsed:",avalue)
        phastankL.append(avalue)
        phastankD[aname] = avalue
        acount          +=1

    print("PHAS in file:%s| PHAS cached:%s" % (acount,len(phastankL)))


    return phastankL,phastankD

def matchPHAS(phasterL,phastankL):
    '''
    Matches results and generated summary
    '''

    print("\nFn: matchPHAS ######################")

    ##$ Output
    outfile     = "matched_%s.txt" % (datetime.datetime.now().strftime("%m_%d_%H_%M"))
    summaryfile = "summary_%s.txt" % (datetime.datetime.now().strftime("%m_%d_%H_%M"))
    fh_out      = open(outfile,'w')
    fh_out.write("phaster-name\tphaster-chr\tphaster-start\tphaster-end\tphaster-pval\tphaster-trig\tphaster-trig\tphaster-trig-loci\tphastank-name\tphastank-chr\tphastank-start\tphastank-end\tphastank-score\tphastank-trig\tmatch-side\n")
    sm_out      = open(summaryfile,'w')

    ### Store
    amatched    = 0   ## Phaster counter for matched
    bmatched    = 0   ## Phastank counter for matched
    aset        = set()  ## Store unique matched in phaster      
    bset        = set()  ## Store unique matched in phasetank

    for aent in phasterL:
        matchflag = False
        aname,achr,astart,aend,apval,atrig,aindex,aloci = aent
        # print("Query:%s | chr:%s | start:%s | end:%s | p-val:%s" % (aname,achr,astart,aend,apval))

        for bent in phastankL:
            bname,bchr,bstart,bend,bpval,btrig = bent

            if achr == bchr:
                ## Check for overlap at 5'-end and 3-end
                start5  = int(astart) - flanksize
                start3  = int(astart) + flanksize
                end5    = int(aend)   - flanksize
                end3    = int(aend)   + flanksize

                if int(bstart) >= start5 and int(bstart)  <= start3:
                    ## Match found
                    matchflag   = True
                    amatched    +=1
                    bmatched    +=1
                    aset.add(aname)
                    bset.add(bname)
                    overlapside = 5
                    # print("-- Matched akey:%s-%s-%s | bkey:%s-%s-%s" % (achr,astart,aend,bchr,bstart,bend))
                    fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (aname,achr,astart,aend,apval,atrig,aindex,aloci,bname,bchr,bstart,bend,bpval,btrig,overlapside))


                elif int(bend) >= end5 and int(bend) <= end3:
                    ## Match found
                    matchflag   = True
                    amatched    +=1
                    bmatched    +=1
                    aset.add(aname)
                    bset.add(bname)
                    overlapside = 3
                    # print("-- Matched akey:%s-%s-%s | bkey:%s-%s-%s" % (achr,astart,aend,bchr,bstart,bend))
                    fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (aname,achr,astart,aend,apval,atrig,aindex,aloci,bname,bchr,bstart,bend,bpval,btrig,overlapside))

                elif (int(bstart) >= int(astart) and int(bend) <= int(aend)) or (int(astart) >= int(bstart) and int(aend) <= int(bend)):
                    ## Match found - Contained PHAS
                    matchflag   = True
                    amatched    +=1
                    bmatched    +=1
                    aset.add(aname)
                    bset.add(bname)
                    overlapside = 3
                    # print("-- Matched akey:%s-%s-%s | bkey:%s-%s-%s" % (achr,astart,aend,bchr,bstart,bend))
                    fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (aname,achr,astart,aend,apval,atrig,aindex,aloci,bname,bchr,bstart,bend,bpval,btrig,overlapside))

                else:
                    ## No match found 
                    pass


            else:
                ## Chromosomes do not match
                pass

        ### catch unmatched and write results
        if matchflag == False:
            ## Write unmatched phaster results to keep the number same in results file
            # print("-- Unmatched akey:%s-%s-%s | bkey:x-x-x" % (achr,astart,aend))
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tx\tx\tx\tx\tx\tx\tna\n" % (aname,achr,astart,aend,apval,atrig,aindex,aloci))

    ### Write unmatched PhaseTank entries
    for bent in phastankL:
        bname,bchr,bstart,bend,bpval,btrig = bent
        if bname not in bset:
            fh_out.write("x\tx\tx\tx\tx\tx\tx\tx\t%s\t%s\t%s\t%s\t%s\t%s\tna\n" % (bname,bchr,bstart,bend,bpval,btrig))

    print("phaster PHAS total   :%s | phastank PHAS total   :%s" % (len(phasterL),len(phastankL)))
    print("phaster matched count:%s | phastank matched count:%s" % (amatched,bmatched))
    print("phaster matched uniq :%s | phastank matched uniq :%s" % (len(aset),len(bset)))

    ### Write Summary file
    sm_out.write("phaster file:%s\n" % (phaster))
    sm_out.write("phastankfile:%s\n" % (phastank))
    sm_out.write("flanksize:%snt\n" % (flanksize))
    sm_out.write("phaster PHAS total   :%s | phastank PHAS total   :%s\n" % (len(phasterL),len(phastankL)))
    sm_out.write("phaster matched count:%s | phastank matched count:%s\n" % (amatched,bmatched))
    sm_out.write("phaster matched uniq :%s | phastank matched uniq :%s\n" % (len(aset),len(bset)))



    fh_out.close()
    sm_out.close()

    return outfile

def main():
    revfernoL,revfernoD   = revferno_parse(revferno)
    phasterL,phasterD   = phaster_parse(phaster,revfernoD)
    phastankL,phastankD = phastank_parse(phastank)
    matchfile           = matchPHAS(phasterL,phastankL)


if __name__ == '__main__':

    main()
    print("\nBye!!\n")
    sys.exit()


### Chage Log
## v01 -> v02
## removes "chromosome" from PhasTAnk results while parsing

## v02 -> v03
## Added score/p-value to output
## Added triggers to output
