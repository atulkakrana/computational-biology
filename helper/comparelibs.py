#!/usr/local/bin/python3

### Script compares two libraries in tab seprated format
### Written bu Atul for Sandra libraries

import sys,os

afile   = sys.argv[1]
bfile   = sys.argv[2]

# afile = "5993.test" 
# bfile = "5993.server.test"



def tagcount_parse(alib):
    '''
    parse tag count files
    '''

    print("\nFn: Parse Tagcount ###########################")

    alist = []
    fh_in = open(alib,'r')
    aread = fh_in.readlines()

    for i in aread:
        aseq,acount = i.strip('\n').split('\t')
        alist.append((aseq.strip(),acount.strip()))


    print("%s entries cached for %s" % (len(alist),alib))

    return alist

def compare(alist,blist):
    '''
    Compare libs
    '''
    print("\nFn: Compare files ###########################")

    auniqL      = []
    buniqL      = []
    commonL     = []
    negList     = [] ## Stores only tags for matching
    for aent in alist:
        aseq    = aent[0]
        acount  = aent[1]
        matFlag = False
        for bent in blist:
            bseq    = bent[0]
            bcount  = bent[1]

            if aseq == bseq:
                commonL.append(aent)
                negList.append(aseq)
                matFlag = True
                break
            else:
                pass

        ## Catch Uniq in A          
        if matFlag == False:
            ## Not found in B
            auniqL.append(aent)

    print(commonL)
    ## Catch Uniq to B
    for bent in blist:
        bseq    = bent[0]
        bcount  = bent[1]
        if bseq not in negList:
            buniqL.append(bent)
        else:
            print("common:%s" % (bent[0]))


    print("Common tags:%s | Uniq to %s:%s | Uniq to %s:%s" % (len(commonL),afile,len(auniqL),bfile,len(buniqL)))

    return auniqL,buniqL,commonL

def writer(auniqL,buniqL,commonL):
    '''
    writes those uniq to A, to B and common
    '''

    print("\nFn: writing files ###########################")

    afileout    = "%s_uniq.txt" % (afile) 
    bfileout    = "%s_uniq.txt" % (bfile)
    cfileout    = "%s_%s_common.txt" % (afile,bfile)

    fh_out1     = open(afileout,'w')
    fh_out2     = open(bfileout,'w')
    fh_out3     = open(cfileout,'w')
    fh_out_log  = open("log.txt",'w')

    acount = 0 
    for i in auniqL:
        fh_out1.write("%s\t%s\n" % (i[0],i[1]))
        acount+=1
    
    bcount = 0
    for j in buniqL:
        fh_out2.write("%s\t%s\n" % (j[0],j[1]))
        bcount+=1

    ccount = 0
    for k in commonL:
        fh_out3.write("%s\t%s\n" % (k[0],k[1]))
        ccount+=1

    print("%s entries common for both files" % (ccount))
    print("%s entries uniq for %s" % (acount,afile))
    print("%s entries uniq for %s" % (bcount,bfile))
    fh_out_log.write("%s entries common for both files" % (ccount))
    fh_out_log.write("%s entries uniq for %s" % (acount,afile))
    fh_out_log.write("%s entries uniq for %s" % (bcount,bfile))


    fh_out1.close()
    fh_out2.close()
    fh_out3.close()

    return None


def main():

    alist = tagcount_parse(afile)
    blist = tagcount_parse(bfile)

    auniqL,buniqL,commonL = compare(alist,blist)
    writer(auniqL,buniqL,commonL)

    pass


if __name__ == '__main__':
    main()
    print("\nScript finished sucessfully\n")
    sys.exit()



### ChangeLog