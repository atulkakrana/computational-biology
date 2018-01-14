#!/usr/local/bin/python3

'''This script is for specific scenario unlike phasiExtract, see read me below'''

import os,sys,time
import operator

fileType            = "cluster.boundary.without.PARE.validation.list"
seqType             = 0                                     ## 0: Genomic coords and normal seq files 1: PacBio - with lots of stuff in header
phasiLenFilter      = 'Y'                                   ## 'Y' then tags of phase length will be written from the cluster | 'N' - All tags will be written
minAbun             = 10                                    ## Minimum abundance of tag to be written

def getBestPhas(fileType):
    ''' phased loci seqeunces were used to re-map sRNA and identify phasiRNAs
    since fasta file had seqeunced fro best/selected phased loci, also non-redundant there
    is no need to check for best loci. In this module only the longest longest report of a phased loci
    is idenifies and then these co-ordinates are used to exatrct phasiRNAs
    '''

    ResFls = [file for file in os.listdir('./') if file.endswith (fileType)]
    resList = []

    for fls in ResFls:
        print("This is the file:",fls)
        file_id = fls.split(".")[0]
        fh_in = open(fls,'r')
        entries = fh_in.readlines()
        for i in entries:
            if i.strip(): ## Remove an empty line from end file
                ent_splt = i.strip('\n').split('=')
                #print(ent_splt[0].split('|'))
                pval,phase,trash = ent_splt[0].strip().split('|')
                chromo_start,end = ent_splt[1].strip().split('..')
                chromo,sep,start = chromo_start.rpartition(':')

                if seqType == 0:
                    resList.append((file_id,phase,float(pval),chromo.strip().replace('chr',''),start,end))
                elif seqType == 1: ## Header has lots of stuff
                    resList.append((file_id,phase,float(pval),chromo.strip(),start,end))
        
        fh_in.close()
    
    print("\nResults list ready with %s values" % (len(resList)))
    print("Example values:%s" % (resList[0:10]))

    ## Sort list on p-values
    resList_sort = sorted(resList, key=operator.itemgetter(2),reverse=False)
    print("\nExample sorted values:%s" % (resList_sort[0:10]))

    ## Make uniq
    resList_names = []     ## List to store uniq names 
    resList_uniq = []      ## Tp store uniq results
    for i in resList_sort:
        name = i[3]
        # print(name)
        if name not in resList_names:
            resList_names.append(name)
            resList_uniq.append(i)
        else:
            pass

    print("\nUniq results list ready with %s values" % (len(resList_uniq)))
    print("Example uniq values:%s" % (resList_uniq[0:10]))

    return resList_uniq

def phasiExtract(resList):

    ## Open output file
    outfile = "extractedPhasi.fasta"
    fh_out = open(outfile,"w")

    ## Sort results on library ###
    resList_sort = sorted(resList, key=operator.itemgetter(1,2),reverse=False)

    ## Get file-specific dict with PHAS name as key and phasiRNAs as values
    fileList = [] ## If file to be read
    for i in resList_sort:
        lib,phase,pval,getname,astart,aend = i
        print("\nThis is the entry:",lib,phase,pval,getname,astart,aend)
        clustfile = "%s.txt.score_p%s_sRNA_%s_out.cluster" %  (lib,pval,phase)
        
        ## Load File ######
        if clustfile not in fileList:
            fileList.append(clustfile)
            print ("File being loaded:%s" % (clustfile))
            clustDict = fileRead(clustfile,phase)


        else:
            print("File for this library-p-value already loaded")
            pass

        ## Extract Entries ###
        phasiClust = clustDict[getname] ## Get the phasiRNAs for this name

        phasiCyc = phasiClust[1]
        phasiSig = phasiClust[2]
        for i in phasiClust[0]:
            phasiname   = i[0]
            phasiseq    = i[1]
            phasilen    = i[2]
            phasiabun   = i[3]
            phasistrand = i[4]
            phasihits    = i[5]
            print("phasiRNAs",phasiname,phasiseq,phasilen,phasiabun,phasistrand,phasihits,phasiSig,phasiCyc)

            if phasiLenFilter == "Y" and phasiabun >= minAbun:
                if phasilen == int(phase):
                    fh_out.write(">%s_%s_%s_%s_%s\t%s\n%s\t%s\n" % (getname,phasiname,phasistrand,phasihits,phasilen,phasiabun,phasiseq,phasiabun))

            elif phasiabun >= minAbun:
                fh_out.write(">%s_%s_%s_%s_%s\t%s\n%s\t%s\n" % (getname,phasiname,phasistrand,phasihits,phasilen,phasiabun,phasiseq,phasiabun))

            else:
                pass

    fh_out.close()

def fileRead(clustfile,phase):

    '''Helper module - Which reas file for first time and returns a dict of PHAS name as key
    and phasiRNA list as value '''

    fh_in = open(clustfile,'r')
    clusters = fh_in.read().split('>')
    fh_in.close()

    ## Raed File and make dictionary ######
    
    clustDict = {} ## dictionary to store phasiRNAs as value and name of PHAS loci as key
    for aclust in clusters[1:]:
        aclust_splt = aclust.split('\n')
        header      = aclust_splt[0].split()
        clust_id    = header[2]
        name        = header[6]
        start       = header[10]
        end         = int(header[12])+1 ## 1 added because when opening a range using start and end, end number is not included in range
        akey = name

        tempPhasiList = [] ## cluster-specific phasiRNA list
        phasiCyc = 0 ## count of cycles
        phasiSig = 0 ## count of abundances
        for ent in aclust_splt[1:-1]:## Because header was the first entry of block and Last entry is always empty
            # print (ent)
            phasient = ent.split('\t')
            phasistrand = phasient[2].translate(str.maketrans("+-","wc"))
            phasiname = phasient[4].replace("|","_")
            phasiseq = phasient[5]
            phasilen = int(phasient[6])
            phasiabun = int(phasient[7])
            phasihits = phasient[10].split("=")[1]
            # print(phasiname,phasiseq,phasilen,phasiabun)
            tempPhasiList.append((phasiname,phasiseq,phasilen,phasiabun,phasistrand,phasihits))

            if int(phasilen) == int(phase):
                phasiCyc +=1
                phasiSig += phasiabun

        aval = (tempPhasiList,phasiCyc,phasiSig)
        clustDict[akey] = aval

    print("Loading complete, %s clusters cached\n" % (len(clustDict)))

    return clustDict

def main():
    resList_uniq = getBestPhas(fileType)
    phasiExtract(resList_uniq)

### RUN ########################
if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print ('Complete run time is %s' % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting..\n')
    sys.exit()

###### LOG ############
## v01 [Nov-5-2015] ###

## v01 -> v02
## Added minAbun, which was not functional before


##### README ##########
### This script was written for a special scenario of Asparagus study. Since, PHAS loci
### were idenified in different studies, extracting their phasiRNAs was a problem.
### So, in this case PHAS loci seqeunces were extracted with some flanking region
### and phasing analysis was run. This was specifically done to extract the phasiRNAs.
### In this script, we choose best library for each phasiRNA locus and then extract 
### phasiRNAs from cluster file of that library
