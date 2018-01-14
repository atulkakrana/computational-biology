#!/usr/local/bin/python3

import os,math,sys,random

## Script to perform random sampling from two sets and output matrix of variable one and variable two. For Ex - cutoffs from non-coding and coding RNAs, and their original class

posSet =  "cod.zma.cpat.txt"          ## Tab-separated file with variable 1 and variable 2
x1Pos =   6          ## Variable 1 like cutoffs from CPAT
x2Pos =   8          ## Variable 2 like labels - noncoding or coding
posHead = 'T'

negSet =  "ncrna.zma.cpat.txt"              ## Tab-separated file with variable 1 and variable 2
y1Pos =   6                                 ## Variable 1 like cutoffs from CPAT
y2Pos =   8                                 ## Variable 2 like labels - noncoding or coding
negHead = 'T'

size = 100                                  ## Size of sample neg+pos
samples = 100                               ## Total number of samples




def sampler(posSet,negSet,size,samples):

    '''Reads postive and negative set files, extracts values, does random sampling and 
    outputs a single list with a balanced random sample as sublist'''

    print("\FUNCTION: sampler")

    samplesList = [] ## List to hold random samples as sublist, and two variables as tuples in sublist

    ## Make list of values in postive set
    posList = []
    fh_in = open(posSet,'r')
    if posHead == 'T':
        fh_in.readline()
    posRead = fh_in.readlines()
    for i in posRead:
        ent = i.strip('\n').split("\t")
        x1 = ent[x1Pos-1] ## Cutoff, or some other value
        x2 = ent[x2Pos-1] ## Label, or some other value
        print("posList:%s,%s"% (x1,x2))
        posList.append((x1,x2))
    print("Total number of entries in positive List:%s" % (len(posList)))
    print("posList 20 elements",posList[:12])

    ## Make list of values in negative set
    negList = []
    fh_in2 = open(negSet,'r')
    if negHead == 'T':
        fh_in2.readline()
    negRead = fh_in2.readlines()
    for i in negRead:
        ent = i.strip('\n').split("\t")
        x1 = ent[y1Pos-1] ## Cutoff, or some other value
        x2 = ent[y2Pos-1] ## Label, or some other value
        # print("negList:%s,%s"% (x1,x2))
        negList.append((x1,x2))
    print("Total number of entries in negative List:%s" % (len(negList)))
    print("negList 20 elements",negList[:20])

    ## Random picking from bothlist
    tempSample = []
    acount = 0  ## Count of sampling events
    for n in range(samples):
        randomPos = random.sample(posList,int(size/2))
        randomNeg = random.sample(negList,int(size/2))
        tempSample.extend(randomPos)
        tempSample.extend(randomNeg)

        ## Add this sample to main random sample list
        samplesList.append(tempSample)
        acount+=1

        print("%s items in sample:%s" % (len(tempSample),acount))
        # print("tempSample",tempSample)
        tempSample=[] ## Empty after every sampling

    fh_in.close()
    fh_in2.close()
    print("Exiting: sampler\n")
    return samplesList

def writer(samplesList,size):

    '''Writes results for both variables in a different matrix'''
    print("\nFUNCTION: writer")

    var1File = "var1Samples.txt"
    var2File = "var2Samples.txt"

    fh_out1 = open(var1File,'w')
    fh_out2 = open(var2File,'w')

    acount = 0  ## Count of postition in sample
    var1List = [] ## List in order for matrix 
    var2List = []

    for i in range(size):
        
        tempvar1 = [] ## Temp list to store variables from samples that are at same position
        tempvar2 = []

        bcount = 0 ## Current sample
        for sam in samplesList:
            # print("Sample:%s"%bcount,sam)
            var1,var2 = sam[acount]
            print(var1,var2)
            tempvar1.append(str(var1))
            tempvar2.append(str(var2))
            bcount += 1

        # print("tempvar1",tempvar1)
        # print("\ntempvar2",tempvar2)

        fh_out1.write("%s\n" % ('\t'.join(x for x in tempvar1)))
        fh_out2.write("%s\n" % ('\t'.join(x for x in tempvar2)))
        
        var1List.append(tempvar1)
        var2List.append(tempvar2)

        acount+=1

    # print("\nFinalList1",var1List)
    # print("FinalList2",var2List)

    print("Files for %s samples with %s elements made\n" % (bcount,acount))

    fh_out1.close()
    fh_out2.close()

    return var1File,var2File

def main():
    samplesList = sampler(posSet,negSet,size,samples)
    var1File,var2File = writer(samplesList,size)

if __name__ == '__main__':
    main()
    sys.exit()

