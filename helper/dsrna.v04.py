#/usr/local/bin/python3

### This script was written by Atul kakrana to analyze dsRNA seq data for maize and asparagus
### To generate input files - see bwtool.sh and thresholdCompute.v01.py

import os,sys,time
import numpy as np
from scipy.stats import binned_statistic

### Uses Settings ####

mode        = 1                     ## Mode:0 PHAS with info for both strands Mode:1 miR with info for one strand - - For mode 0, file should have same number of w and c entries in repFiles 
repFile1    = "181.w.c.txt" ## File generated with bwtools extract function for regions of interest
repFile2    = "182.w.c.txt" ## File generated with bwtools extract function for regions of interest
nbins       = 100
upthres     = 2.08          ## Upper threshold, usually 97.5th percentile of bigwig signals, or signals of specific regions, -These are computed from "thresholdCompute" script
lowthres    = -2.08         ## Lower threshold, usuall lower 2.5th percentile of bigwig signals, or signals of specific regions -  These are computed from "thresholdCompute" script 

stattype    = 0                       ## 0: Max datapoint represents bin score 1: Mean 2: Sum 
# binType     = 0                     ## 0: Unoform binning i.e. last bin will take values from second last bin to keep size uniform 1: Last bin score will be based on remainder, this could be of any size
#####################
#### Functions #####
def collapseScores(repFile1,repFile2):
    '''This is master module the merges the scores from different strands, libraries and 
    bins them to certain length so that different regions could be alinged'''

    mergedScores = [] ## List of merged and binned scores

    ### Read repFiles
    afileList,awstrandD,acstrandD = fileReader(repFile1)
    bfileList,bwstrandD,bcstrandD = fileReader(repFile2)

    ### Merge locus-specific data points from different libraries
    astrandmergeD = mergeStrandScores(awstrandD,acstrandD,afileList) ## Library one
    bstrandmergeD = mergeStrandScores(bwstrandD,bcstrandD,bfileList) ## Library two

    ### Merge locus-specific data points from different libraries
    mergedscoresD = mergeLibScores(astrandmergeD,bstrandmergeD) ## Mode 0 will have keys with names of regions and mode 1 will have unique strand specific keys
    outFile       = writer(mergedscoresD)

    ### Bin scores to make all features comaprable
    binnedList=binScores2(mergedscoresD) ## [(Name,merged data points),(),...]
    outFile = writer2(binnedList)

    return mergedScores

def binScores(mergedscoresD):
    '''This function makes the regions of different length comaprable by binning them'''

    binnedList = [] ### Stores result after binnnig, remeber that metric from last bin will 
    ####            include data points from second last bin

    for aname in mergedscoresD.keys():
        print("Ent:%s" % (aname))
        datap   = mergedscoresD[aname]
        datalen = len(datap)
        binsize = int(datalen/nbins) ## Here int() will return number by ignoring point values so 1.6 will be 1, remander is caught below
        lastbinvals = int(binsize%nbins)
        print("Total datapoints:%s | binsize:%s | last bin size:%s" % (datalen,binsize,lastbinvals))

        ### Sanity check that no sequenc ewith length smaller then bin size is used

        if datalen < nbins:
            print("The length of data is smaller then required bins - Exiting")
            sys.exit()

        binstart = 0 ### This keeps tracck of bin start

        tempscoreL = [] ## List to hold means scores of the bins
        for abin in range(nbins):
            binend  = binstart+binsize
            bindata = datap[binstart:binend]
            print("bindata",bindata)
            
            binscore    = 0
            validvalues = 0
            for i in bindata:
                if i != 'NA':
                    binscore+= float(i)
                    validvalues+=1
                else:
                    print("NA skipped")

            binscore = round(binscore/validvalues,2)
            tempscoreL.append(binscore)
            print("This is binscore:%s" % (binScores))

        ### Score computation complete for all perfectly sized bins, now scores for last bin needs to be computed
        if lastbinvals > 0:
            if binType == 0:
                ### Make last bin uniform as others, take some values froms econf last bin
                extravals = binsize-lastbinvals

                x = abin[-extravals:]
                y = datap[-lastbinvals:]
                print("These are values taken from second last bin",x)
                print("These are values from last bin",y)
                print("Required binsize:%s | Vals from second last bin:%s | Vals. from last bin:%s" % (binsize,len(x),len(y)))
                bindata = x+y

                binscore    = 0
                validvalues = 0
                for i in bindata:
                    if i != 'NA':
                        binscore+= float(i)
                        validvalues+=1
                    else:
                        print("NA skipped")
                binscore = round(binscore/validvalues,2)
                tempscoreL.append(binscore)
                print("This is binscore:%s" % (binScores))


            elif binType == 1:
                ### Compute bin score just based on remainder data points, the binsize is not uniform in this case
                y = datap[-lastbinvals:]

                print("These are values from last bin",y)
                print("Required binsize:%s |Vals. from last bin:%s" % (binsize,len(y)))
                bindata = y

                binscore    = 0
                validvalues = 0
                for i in bindata:
                    if i != 'NA':
                        binscore+= float(i)
                        validvalues+=1
                    else:
                        print("NA skipped")
                binscore = round(binscore/validvalues,2)
                tempscoreL.append(binscore)
                print("This is binscore:%s" % (binScores))


            else:
                print("Please input correct bin type - exiting")
                sys.exit()

        print("score computations for all bins of %s complete" % (aname))
        binnedList.append((aname,tempscore))


    return binnedList

def binScores2(mergedscoresD):
    '''This function bins score using scipy libraries, that are capable of merging into bins to avoid maximium uniformity and cover complete range'''

    print("\nFN - binScores")

    binnedList = [] ### Stores result after binnnig, remeber that metric from last bin will 
    ####            include data points from second last bin

    misscount = 0 ## Count of those smaller then binsize
    for aname in mergedscoresD.keys():
        print("Ent:%s" % (aname))
        adatap   = mergedscoresD[aname]
        print("Data points:",adatap)

        ## Remove NAs before binning
        adata    = []
        for x in adatap:
            if x == 'NA':
                adata.append(0)
            else:
                adata.append(x)

        ### Get the positions as list to be used for binning
        adatalen = len(adata)
        avalues  = list(range(1,adatalen+1,1))
        print("-Transformed data for binning:",adata)
        print("-Positions for binning:",avalues)
        print("-Length of data:%s | Length of positions array:%s" % (len(adata),len(avalues)))

        ## Do the binning 
        binstat,binedge,binids = binned_statistic(avalues, adata,statistic='mean', bins=nbins)
        # print("-binstat",binstat)
        # print("-bin assignments",binids)

        ### Sanity check that no sequence with length smaller then bin size is used
        if (adatalen+1) < nbins:
            misscount+=1
            print("The length of data is smaller then required bins - Exiting")
            # sys.exit()
            pass

        ## Merge two data points with bin ids
        binneddata = list(zip(adata,binids))
        print("+This is snippet of binned data",binneddata[0:20])

        ### Capture data for bins into a single list 
        binwisevals = [] 
        for anid in set(binids):
            # print("Checking values for bin:%s" % anid)
            tempdata = [] ## Store data for a bin to copute stats
            for x in binneddata:
                # print("Bin",x)
                aval,abin = x
                if abin == anid:
                    # print("match")
                    tempdata.append(float(aval))
                else:
                    pass
                    # sys.exit()

            binwisevals.append(tempdata)
        
        print("+This is the snippet of binwisevalues:",binwisevals[0:20])

        ### Compute bin wise stats
        binwisestats = [] ## This will hold final stat of every bin - max, sum, mean (only if a value)
        for abin in binwisevals:
            print("\n-abin data points:",abin)
            nvals       = len(abin)                       ## All values
            validdata   = list(x for x in abin if x != 0)   ## All non-zero values
            nvalid      = len(validdata)
            print("Total values:%s | data points:%s" % (nvals,nvalid))
            if validdata:
                maxstat     = max(validdata)
                meanstat    = round(sum(validdata)/nvalid,2)
                sumstat     = round(sum(validdata),2)
            else:
                print("No data point recorded- All are 'NA's'")
                maxstat     = 0
                meanstat    = 0
                sumstat     = 0
            print("Stats - max:%s | mean:%s | sum:%s " % (maxstat,meanstat,sumstat))

            if stattype == 0:
                binwisestats.append(maxstat)
            elif stattype == 1:
                binwisestats.append(meanstat)
            elif stattype == 2:
                binwisestats.append(sumstat)
            else:
                print("Please choose correct 'stattype'")
                pass

        ### Add ninned stats and loci name to results
        print("Snippet of %s binned stats:%s" % (aname,binwisestats[0:20]))
        binnedList.append((aname,binwisestats))
    print("+Features skipped due length samller then bin size:%s" % (misscount))
    time.sleep(1)
    print("FN Exiting - binScores")

    return binnedList

def mergeLibScores(adict,bdict):
    '''This function merges the scores from multiple libraries - If data point recorded by both at a position then average of data points is taken other wise recoded data point is taken'''

    print("\nFN - mergeLibScores")

    ### Getting list of names, assuming that both replicate libraries have same features, it' just values are diffrerent
    features = adict.keys()

    print("Merging Library scores")
    mergedscoresD = {} ## Store name as keys and lib merged data points as values
    for aname in features:
        # print("Ent:%s" % (aname))
        adatap = adict[aname]
        bdatap = bdict[aname]

        acount      = 0
        tempscores  = [] ## List to store merged scores for a loci
        
        ### Length of both lists should be same because smae feature will have same length
        for a,b in zip(adatap,bdatap):
            acount+=1
            # print("pos %s values:%s, %s" % (acount,a,b))
            
            if a == 'NA' or b == 'NA':
                ## Take the available data point as is
                if a == 'NA':
                    ascore = b
                else:
                    ascore = a
                # print("- %s single data point:%s,%s" % (aname,a,b))
                tempscores.append(ascore)
            else:
                # ascore = round((float(a)+float(b))/2,3)
                ascore = max([float(a),float(b)])
                # print("- %s two data points found - average: %s" % (aname,str(ascore)))
                tempscores.append(ascore)
                # sys.exit()
                pass
        
        mergedscoresD[aname] = (tempscores)
        # sys.exit()
    print("FN Exiting - mergeLibScores")


    return mergedscoresD

def mergeStrandScores(wdict,cdict,afileList):
    '''This function merges strand scores from different strands, by taking average of data points, else taking max i.e. datapoint'''

    print("\nFN - mergeStrandScores")

    strandmergeD = {} 
    print("Merging 'w' and 'c' strands")

    if mode == 0:
        ### Collect a unique set of names
        alociS = set() ## Set of unique loci irrespective of strand to use as key for w and c dictionaries
        acount = 0 ## To count the numbe rof entries in the w/c file list
        for i in afileList:
            acount  +=1
            aname   = i[3]
            alociS.add(aname)
        print("Total entries for w and c:%s | unique entries:%s" % (acount,len(alociS)))

        print("Both 'w' and 'c' data points exists and will be merged")
        for i in alociS:
            # print("Ent:%s" % (i))
            aname   = i
            wkey    = "%s_w" % (aname)
            ckey    = "%s_c" % (aname)
            wlist   = wdict[wkey]
            clist   = cdict[ckey]

            acount  = 0
            tempscores = [] ## List to store merged scores for a loci
            for w,c in zip(wlist,clist):
                acount+=1
                # print("pos %s values:%s, %s" % (acount,w,c))
                if w == 'NA' or c == 'NA':
                    ## Take the available data point as is
                    if w == 'NA':
                        ascore = c
                    else:
                        ascore = w
                    tempscores.append(ascore)
                else:
                    ascore = round((float(w)+float(c))/2,3)
                    # print("Two data points found - average: %s" % (str(ascore)))
                    tempscores.append(ascore)
                    # sys.exit()
                    pass
            strandmergeD[aname] = (tempscores)

    if mode == 1:
        # print("Either 'w' or 'c' data point exists - No merging required")
        ## Combining dictionaries into one dict
        strandmergeD = wdict.copy()
        strandmergeD.update(cdict)


        # print("Scores:",strandmergedL)
        # sys.exit()
    print("FN Exiting - mergeStrandScores")


    return strandmergeD

def fileReader(afile):
    '''This reads the output of bwtools, generated with bed6 file. First six columns
    are from bedfile and 7th column is length of region, 8 th column has comma seprated 
    data points. Data for both strands (if exists) is stored in same file'''
    print("\nFN - fileReader")

    wstrandD    = {} ### Holding entries of positive strand
    cstrandD    = {} ### Holding strands for negative strand
    fileList    = [] ### All records in alist

    fh_in = open(afile,'r')
    aread = fh_in.readlines()

    print("Caching file:%s" % (afile))
    for i in aread:
        ent = i.strip('\n').split('\t')
        fileList.append((ent))
        aname   = ent[3]
        astrand = ent[5].translate(str.maketrans("+-","wc"))
        datap   = ent[-1].replace('"','').split(',') ### Its a list
        # print("datapoints:",datap)
        if astrand == 'w':
            akey           = "%s_%s" % (aname,astrand) ## Strand specific key to help in keeping these as seprate
            wstrandD[akey] = (datap)
            # wstrand.append((aname,astrand,datap))
        elif astrand == 'c':
            # print(datap)
            # print(datap[::-1])
            akey            = "%s_%s" % (aname,astrand) ## Strand specific key to help in keeping these as seprate
            cstrandD[akey]  = (datap[::-1]) ## For crick strand bwtools reverses the positions, setting them back
            # cstrand.append((aname,astrand,datap))
        else:
            print("Unknown stand encountered - Debug")
            sys.exit()


    print("Total entries:%s | Strand-w:%s | strand-c:%s" % (len(fileList),len(wstrandD),len(cstrandD)))
    print("FN Exiting - fileReader")
    # sys.exit()

    return fileList,wstrandD,cstrandD

def writer(mergedscoresD):

    '''Write merged scores'''

    outFile = "merged.txt"
    fh_out = open(outFile,'w')

    for i in mergedscoresD.keys():
        aname = i
        datap = mergedscoresD[aname]
        fh_out.write("%s\t%s\n" % (aname,'\t'.join(str(x) for x in datap)))

    fh_out.close()

    return outFile

def writer2(binnedList):

    '''Write binned scores'''
    
    print("\nFN - writer2")

    outFile = "merged.bin.txt"
    fh_out  = open(outFile,'w')
    binnames= ['bin']*nbins
    if mode == 0:
        fh_out.write("Name\tmaxstat\tsumstat\tmeanstat\t%s\n" % ('\t'.join(x for x in binnames)))
    elif mode == 1:
        fh_out.write("Name\tstrand\tmaxstat\tsumstat\tmeanstat\t%s\n" % ('\t'.join(x for x in binnames)))
    else:
        print("Please input correct 'mode'")
        pass

    for i in binnedList:
        print("Ent:",i)
        aname,datap = i
        maxstat     = max(datap)
        sumstat     = round(sum(datap),2)
        validdata   = list(x for x in datap if x != 0)   ## All non-zero values
        if validdata:
            meanstat = round(sumstat/len(validdata),2)
        else:
            meanstat = 0
        
        if mode == 0: 
            fh_out.write("%s\t%s\t%s\t%s\t%s\n" % (aname,maxstat,sumstat,meanstat,'\t'.join(str(x) for x in datap)))
        elif mode ==1:
            bname,trash,bstrand   = aname.rpartition("_")
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (bname,bstrand,maxstat,sumstat,meanstat,'\t'.join(str(x) for x in datap)))
        else:
            print("Please input correct 'mode'")
            pass

    fh_out.close()

    return outFile

##### MAIN ############
#######################

def main():
    mergedScores = collapseScores(repFile1,repFile2)


if __name__ == '__main__':
    main()
    print("\n\nI_I` Free beer ---> ")


## v01 -> v02
## Added binning capability of dsRNA data to make be able to compare between different features - No all fearures are fo same bins, rather then different lenghs

## v02 -> v03
## Added functionality to could those cases where the feature length is smaller then bin size, and skip them

## v03 -> v04
## Added functionality to keep scores of PHAS from different strands separate and generate strand specific bins
## Fixed a bug in mergeStrand function, where list for input files, holding same PHAS for both strands was itnerated on PHAS twice
##
