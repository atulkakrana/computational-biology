#!/usr/local/bin/python3

import sys,os,time

###############

afile = "24AlldistMat.ClustalO.melt.tsv"
head = 1            ## 0: No and 1:Yes
sep = '\t'          ## Separator
simThres = 40       ## Threshold for similarity/distance (non-negative)
nThres = 1          ## Number of alignments above cutoff
################


def filterMat(afile,simThres,nThres):
    ''' distance/similarity matrix in melted format i.e.
    x,y and z where x and y are seqeunce names and z is the distance/similairty 
    between them'''
    fh_out = open("failedNames.txt",'w')
    fh_in = open(afile,'r')
    if head == 1:
        fh_in.readline()
    aread = fh_in.readlines()

    nameSet = set()
    dataList = [] ## Data stored in list format for later operations
    for i in aread:
        # print(i.strip().split(sep))
        x,y,z = i.strip().split(sep)
        nameSet.add(y)
        dataList.append((x,y,z))

    print("Total unique elements in input file:%s" % len(nameSet))

    finalSet = set()  ## Final names that have atleast one alingment above cutoff, other then self
    failSet = set()   ## Names that do not have even a single alignment above cutoff,these will be deleted
    
    for seqname in nameSet:
        acount = -1          ## Number of pairwise comparision greater then cutoff, 1 is expected that is with itself
        
        for i in dataList:
            x,y,z = i ## x and y are seqeunces and  value is similarity between them

            if y == seqname:
                if float(z) >= simThres:
                    acount += 1
            else:
                pass

        if acount >= nThres:
            finalSet.add(seqname)
        else:
            failSet.add(seqname)

    ## Write failes names that did not pass cutoff

    for i in failSet:
        fh_out.write("%s\n" % (i))

    fh_in.close()
    fh_out.close()


    print("Total sequences %s that have atleast %s alingment with similarity/distance >= %s" % (len(finalSet),nThres,simThres))

    return dataList,finalSet

def writer(dataList,finalSet):
    outFile = "%s_sim%s_%sn.txt" % (afile,simThres,nThres)
    fh_out = open(outFile,'w')

    acount = 0
    for i in dataList:
        x,y,z = i
        if y in finalSet: 
            if x in finalSet:
                fh_out.write("%s\t%s\t%s\n" % (x,y,z))
                acount +=1

    ## Check if correct entries are written
    expEnt = len(finalSet)*len(finalSet)
    if acount == expEnt: ## That is all final entries had alingment with all final entries
        print("Correct number of entries written: %s" % (str(expEnt)))
    else:
        print("Expected number of entries %s were not written" % (str(expEnt)))
        print("Script will exit now")
        sys.exit()

    fh_out.close()

    return outFile

def main():
    dataList,finalSet = filterMat(afile,simThres,nThres)
    outFile = writer(dataList,finalSet)

if __name__ == '__main__':
    start = time.time()
    main()
    print ('It took', round(time.time()-start,2), 'sec')
    print ('\nCheers script finished sucessfully\n')
    sys.exit()


## Jun-24-15
