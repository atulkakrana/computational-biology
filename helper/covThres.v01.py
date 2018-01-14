#!/usr/local/bin/python3

### This script does random sampling to compute pthreshold percentiles

### Imports ###############
import os,sys,random,subprocess,multiprocessing
import numpy as np
# from statistics import median

##### User settings ########
coordsFile  = "Aspa.v2.scaffolds.txt" ## Tab separated file with chr, start and end - Make sure chr has suffix matching feature file
bigwigw    = "asparagus_bm14_181_structure.plus.bw"
bigwigc    = "asparagus_bm14_181_structure.minus.bw"

featureFile = "24PHAS.v9.collapse.w.bed.txt" ## File with features of variable length whose threshold needs to be computed, strand doesnt matter as all  we need is name and length of feature- random samples of same length will be drawn
delimiter   = "\t"
namepos     = 4 ## Name position in excel format
achrpos     = 1
startpos    = 2 ## Start position in excel format
endpos      = 3 ## End position in excel format

niter       = 10
nsamples    = 100
stats       = 1     ## 0: Mean | 1: Median - This will be used compute single percentiles from list of percentiles from diffrent itnerations


##### Functions ############
def thresCompute(coordsFile,featureFile):
    '''Computes length-specific coverage threshold'''
    ## parse coords file
    coordsD     = coordsparser(coordsFile)
    ## parse input features file
    featuresD   = featureparser(featureFile)
    ## sample random regions
    randcoordsD = sampler(coordsD,niter,nsamples,featuresD)

    ### Output file
    if stats == 0:
        fh_out = open("coveragethres.mean.txt","w")
    elif stats == 1:
        fh_out = open("coveragethres.median.txt","w")
    else:
        print("Check 'stat' user setting - Exiting")
        sys.exit()

    fh_out.write("Name\tPerc95\tperc5\n")

    ### Get data points for this random sampled data and compute 95th and 5th percentiles
    for afeature in randcoordsD.keys():
        print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print("+++++++Computing percentile for feature:%s" % (afeature))
        randdata = randcoordsD[afeature]

        perc95L   = [] ## Percentile from all iters for this feature/PHAS 
        perc5L    = [] ## Percentile from all iters for this feature/PHAS
        acount    = 0 ## Keep track of itneration 
        for aniter in randdata:
            acount+=1
            # print("\nSamples for itneration-",acount,aniter)
            bedw,bedc   = bedwriter(aniter,afeature)
            resw,resc   = bwtool(bedw,bedc,bigwigw,bigwigc)
            alldatap    = bwparser(resw,resc)
            aperc95,aperc5= perccompute(alldatap) ## Get percentile for this iter, append to perclist
            perc95L.append(aperc95)
            perc5L.append(aperc5)

        ## 95th and 5th percentile of coverages for this feature recorded for all itnerations, get final percentiles
        if stats == 0:
            perc95 = round(sum(perc95L)/len(perc95L),2)
            perc5  = round(sum(perc5L)/len(perc5L),2)
        elif stats == 1:
            perc95 = round(np.median(np.array(perc95L)),2)
            perc5  = round(np.median(np.array(perc5L)),2)
        else:
            print("Input correct 'stat' to combine percentiles from multiplel itnerations - Exiting")
            sys.exit()

        achr,alen  = featuresD[afeature] ## Just to get an idea of len while looking at threholds in next print statment
        print("Coverage percentiles for %s of length %s are %s and %s" % (afeature,alen,perc95,perc5))
        fh_out.write("%s\t%s\t%s\n" % (afeature,perc95,perc5))
    
    fh_out.close()

    return None

def perccompute(reslist):
    '''
    This function takes a list of datapoints, and return top and bottom percentile of coverages found. Remember that these percentile corresponds to feature length and will be different for different lengths
    '''
    print("\nFn: perccompute #########")

    ### Get the number of data points for randomly sampled regions of this itneration

    foundcovL = [] ## Alist that stores number of datapoints found from list of randomly slected regions - This is later used to compute 95th and 5th percentile of covereages for this feature of specific length.
    ###
    for ares in reslist:
        # print("Sampled region",ares)
        ndatap = 0
        for x in ares:
            if x != "NA":
                ndatap+=1
            else:
                # print("No data point recorded at this position")
                pass
        foundcovL.append(ndatap)
    print("Coverages for this itneration",foundcovL) ## datapoints found in each sample, length of list is equal to nsamples, data points per samples

    ### Compute percentiles
    aperc95 = np.percentile(foundcovL, 95) ##return 50th percentile, e.g median
    aperc5  = np.percentile(foundcovL, 5)
    print("95th perc:%s | 5th perc:%s" % (aperc95,aperc5))

    print("Fn: perccompute exit #########")
    return aperc95,aperc5 

def bwparser(resw,resc):
    '''
    This function reads w and c results from bwtool and return a single list with data points
    '''
    print("\nFn: bwparser #########")
    alldatap = [] ## List that will hold data points from both w and c strands

    ## Parse data points from 'w' results file
    fh_inw = open(resw,'r')
    wread = fh_inw.readlines()
    fh_inw.close()
    acount = 0 ## COunt entries in w file
    for i in wread:
        ent = i.strip("\n").split("\t")
        datap = ent[-1]
        datapL = datap.split(",")
        alldatap.append(datapL)
        acount +=1

    ## Parse data points from 'c' results file
    fh_inc = open(resc,'r')
    cread = fh_inc.readlines()
    fh_inc.close()
    bcount = 0
    for i in cread:
        ent = i.strip("\n").split("\t")
        datap = ent[-1]
        datapL = datap.split(",")
        alldatap.append(datapL)
        bcount+=1

    print("Entries is w file:%s | Entries in c file:%s | Entries in combined result list:%s" % (acount,bcount,len(alldatap)))

    print("Fn: bwparser exit #########")

    return alldatap

def bwtool(bedw,bedc,bigwigw,bigwigc):
    '''
    given a bed file runs bwtool extract function
    '''
    print("\nFn: bwtool #########")

    ## Fetch data for w strand
    # print("Input bed:%s | bigwig:%s" % (bedw,bigwigw))
    resw    = "tempresw.txt"
    retcode = subprocess.call(["bwtool","extract","bed",bedw,bigwigw,resw])
    if retcode == 0:
        # print("-Data extracted sucessfully for %s 'w' strand" % (bedw) )
        pass
    else:
        print("Something wrong happened while extracting from 'w' bigwig - Debug for reason")
        print("Make sure bed file has first column with 'chr' or without it - as expected by bigwig")## The feature file and coords file should have chr_id as expected by bigwig file
        sys.exit()
    
    ## Fetch data for c strand
    # print("Input bed:%s | bigwig:%s" % (bedc,bigwigc))
    resc = "tempresc.txt"
    retcode = subprocess.call(["bwtool","extract","bed",bedc,bigwigc,resc])
    if retcode == 0:
        # print("-Data extracted sucessfully for %s 'c' strand" % (bedc) )
        pass
    else:
        print("Something wrong happened while extracting from 'c' bigwig - Debug for reason")
        print("Make sure bed file has first column with 'chr' or without it - as expected by bigwig") ## he feature file and coords file should have chr_id as expected by bigwig fileread
        sys.exit()

    print("Fn: bwtool exit #########")

    return resw,resc

def bedwriter(alist,featurename):
    '''
    This function writes temporary bed file
    chr,start,end,name,score,strand
    '''
    
    print("\nFn: bedwriter #########")

    wtempbed    = "temp.w.bed"
    ctempbed    = "temp.c.bed"
    fh_outw     = open(wtempbed,'w')
    fh_outc     = open(ctempbed,'w')

    for i in alist:
        # print("-acoord:",(i))
        achr,astart,aend,astrand = i
        if astrand == "+":
            fh_outw.write("%s\t%s\t%s\t%s\t0\t%s\n" % (achr,astart,aend,featurename,astrand))
        elif astrand == "-":
            fh_outc.write("%s\t%s\t%s\t%s\t0\t%s\n" % (achr,astart,aend,featurename,astrand))
        else:
            print("-The strand is unrecognized- Debug",i)
            sys.exit()


    fh_outw.close()
    fh_outc.close()

    # sys.exit()
    print("Fn: bedwriter exit #########")
    return wtempbed,ctempbed

def sampler(coordsD,niter,nsamples,featuresD):
    '''
    Generates a dict of randomly picked regions from specific chromosome equal to sample size, and repeat the same for number of itnerations
    '''
    print("\nFn: sampler #########")
    strandL = ['+','-'] ## Used later to select strand
    randcoordsD = {} ## Dictionary that holds coords (eualing sample size) x nitenerations i.e.[(List from itneration1),(itneration2)...]
    for afeature in featuresD.keys():
        print("+Random sampling for %s" % (afeature))
        reqchr,reqlen = featuresD[afeature] ## Chromosome-specific sampling is not required in case of asparagus, becaus ethese have scaffolds

        featurecoords = [] ### Store random coords from all itnerations
        for aniter in range(niter):
            print("-iter-%s" % aniter)
            ### Draw random chromosomes - double the same size, just to keep some extra
            iterL       = [] ## Random coords for itneration 
            randchrL    = random.sample(coordsD.keys(),nsamples*2) ## Draw chromosomes for required sasmple size
            for arandchr in randchrL:
                ### Now pick random location withiin chromosome
                chrstart,chrend = coordsD[arandchr]
                randstartL   = random.sample(range(int(chrstart), int(chrend)), 1)
                # print(randstart) ## It is a list 
                randstart   = randstartL[0]
                randend     = randstart+reqlen
                randstrand  = random.choice(strandL)
                if randend < chrend:
                    # print("Valid region")
                    iterL.append((arandchr,randstart,randend,randstrand))
                else:
                    # print("-This random region %s-%s goes beyond chr end %s - skipping" % (randstart,randend,chrend))
                    pass
            
            print("-iter-%s complete" % (aniter))
            featurecoords.append(iterL[:nsamples])
            # print("Length of feature coords for %s is %s" % (afeature,len(featurecoords)))

        ## All itnerations complete for this feature, appending results
        randcoordsD[afeature] = (featurecoords)
        # print("random dict snippet",randcoordsD[afeature])

    print("Random sample dictionary has coords for  %s entries" % (len(randcoordsD)))
    

    print("Fn: sampler exit #########")
    return randcoordsD

def coordsparser(coordsFile):
    '''Takes the coords file and parses to give a dict with chr as key and start, end as values
    '''
    print("\nFn: coordsparser #########")
    fh_in = open(coordsFile,'r')
    fh_in.readline() ## Remove Header
    fileread = fh_in.readlines()

    coordsD = {} ## [name/chr],[start,stop]
    for i in fileread:
        ent = i.strip('\n').split('\t')
        achr,astart,aend = ent
        coordsD[achr] = (int(astart),int(aend))

    print("Coords dict with prepared with %s entries" % (len(coordsD)))
    print("Fn: coordsparser exit #########")
    
    return coordsD

def featureparser(featureFile):
    '''
    reads feature file providea a dict with name as key and length as value
    '''
    print("\nFn: featureparser #########")
    fh_in = open(featureFile,'r')
    fh_in.readline() ## Remove Header
    fileread = fh_in.readlines()

    featuresD = {} ## [name/chr],[start,stop]
    for i in fileread:
        ent     = i.strip('\n').split('\t')
        aname   = ent[namepos-1]
        achr    = ent[achrpos-1] ## The feature file should have chr_id as expected by bigwig file
        astart  = int(ent[startpos-1])
        aend    = int(ent[endpos-1])
        alen    = int(aend-astart) 

        ## Sanity check
        if alen <= 0:
            print("The length of %s is %s - Please debug" % (aname,alen))
            sys.exit()
        
        featuresD[aname] = (achr,alen)

    print("Features dict with prepared with %s entries" % (len(featuresD)))
    print("Fn: featureParser exit #########")

    return featuresD

##### MAIN ###########
def main():
    thresCompute(coordsFile,featureFile)
    pass

if __name__ == '__main__':
    main()
    print("\n\n 'A still mind can make whole world surrender'\n")


## v01 - 11th March
## We can use genecoords, or some other transcripts
## Median or MAD (mean absolute deviation) percentiles from multiple itneration can be used instead of mean