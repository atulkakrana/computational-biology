#!/usr/local/bin/python3

## Requires fastQC, trimmomatic and tally in $PATH variable
## Script written by ATUL Kakrana: kakrana@udel.edu

## run : python3 ScriptName.py

############## IMPORTS #####################

import sys,os,re,time,timeit,csv,glob, string
import subprocess, multiprocessing
import shutil
import itertools as it
import operator
from multiprocessing import Process, Queue, Pool
import mysql.connector as sql
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import datetime

### PPRE-PROCESSING SETTINGS ##########################################
genomeDB = 'MAIZE_AGPv2_genome'

#libs = [3726,3727,3225,3226,3486,3519,3520,3521,3522,3523,3481,3482,3715,3716,3717,3718,3722,3723,3227,3228,3728,3733,3724,3725,3729,3734,3229,3230,3730,3735,3731,3736,3732,3737,]
libs = [4179,4181,4180,4182,4178,4185,4183,4186,4184,4187,4233,4234,4235]
#reps = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,18,19,19,20,20,21,21,22,4,5] ## for final step i.e. normalizing abundances

#Optional Steps
QCheckStep = 0 ## Optional -Performs preliminary QC
preProGraphsStep = 0 ## Generates before chopping graphs3228

## Required Step
trimLibsStep = 0 ## Trim fastq files
chopLibsStep = 1 ## Chops
fastQ2CountStep = 0 ## Converts chopped to tag count
mapperStep = 0 ## Maps final chopped files and generates graphs
cleanupStep = 0 ## Final cleanup
###########################################################

############# ADVANCED SETTINGS ###########################
Local = 'Y' ## Run local Analysis: Depends upon maxLen,minLen, maxReadLen, adpaters.fa and files locally
#genoIndex = 'xyz'
#maxLen = 32 ## Max length of the tag allowed
#minLen = 18 ## Min length of tag allowed
maxReadLen = 100 ## max allowed unchopped read length for graph generation

numProc = 0 ## Coarse grain PP
nthread = 10 ## Fine grain PP

masterDB = 'master'
dataServer = 'raichu.dbi.udel.edu'
##########################################################

############# RNA RUNNER ##################################
##Files
genoFile = '/home/kakrana/99.genomes/maize.AGP.v3/Zea_mays.AGPv3.22.dna_sm.genome.fa'
gffFile = '/home/kakrana/99.genomes/maize.AGP.v3/Zea_mays.AGPv3.22.gtf' ## http://plants.ensembl.org/info/website/ftp/index.html
libType = 0 ## From cuffNorm manual 1) 0 = fr-unstranded 2) 1 = fr-firststrand 3) 2 = fr-secondstrand

##Tools - updated versions
tophat = '/home/kakrana/tools/tophat-2.0.12/tophat2'
cufflinks = '/home/kakrana/tools/cufflinks-2.2.1/cufflinks'
cuffmerge = '/home/kakrana/tools/cufflinks-2.2.1/cuffmerge'
cuffquant = '/home/kakrana/tools/cufflinks-2.2.1/cuffquant'
cuffnorm = '/home/kakrana/tools/cufflinks-2.2.1/cuffnorm'

##Steps
indexBuilderStep = 0
spliceMapperStep = 1
cuffLinksStep = 1
cuffMergeStep = 0
cuffQuantStep = 0
cuffNormStep = 0

############# FUNCTIONS #######################

## Output: "libName_fastqc" folder
def QCheck(aninput):
    print(aninput)
    lib,ext,nthread,infile = aninput
    print('****Checking quality of %s library****' % (lib))    
    toolPath = "%s/svn/Tools/FastQC/fastqc" % (os.getenv('HOME'))
    retcode2 = subprocess.call([toolPath,infile])

    return None

## Output: "libName.trimmed.fastq"
def trimLibs(aninput):
    print(aninput)
    lib,ext,nthread,infile,adp_5p,adp_3p,minTagLen = aninput
    print('****Trimming %s library with min length %s****' % (lib,minTagLen))
    
    trimmedFile = '%s.trimmed.%s' % (lib,ext)
    
    adapter = open('%s_adapter.fa' % (lib), 'w')
    adapter.write('>adapter_5p\n%s\n>adapter_3p\n%s' % (adp_5p,adp_3p))
    adapter.close()
    
    ## Just Trim
    #trimlog = '%s.trim.log' % (lib)
    
    toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:./%s_adapter.fa:2:30:10" % (lib), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:%s" % (minTagLen)])
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n****Trimming for %s complete****' % (infile) )
    else:
        print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
        sys.exit()
        
    
    ### Make plot
    #charts(mappedList,mappedAbunList,allTagsList,allAbunList,mode)
    
    return None

## Output: "libName.chopped.fastq" 
def chopLibs(aninput):
    print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('****Chopping %s library to length %s****' % (lib,maxTagLen))
    
    #choplog = '%s_chop.log' % (lib)
    trimmedInFile = '%s.%s' % (lib,ext) ####
    choppedOutFile = '%s.chopped.%s' % (lib,ext)
    
    print("\n")
    toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (maxTagLen)])
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Chopping for %s complete****' % (trimmedInFile) )
    else:
        print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
        sys.exit()
    
    return None

## Output: "libName.processed.txt"
def fastQ2Count(aninput):
    print(aninput)
    lib,ext,nthread = aninput
    print('\n****Converting %s_%s file to tag count****\n' % (lib,ext))
    infile = '%s.%s' % (lib,ext)
    outfile = '%s.processed.txt' % (lib)
    retcode = subprocess.call(["tally", "-i", infile, "-o", outfile, "--nozip", "-format","%R%t%X%n"])
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Conversion to tag count format for %s complete****' % (infile) )
    else:
        print("Something wrong happened while converting to %s library tag count - Debug for reason" % (lib))
        sys.exit()

## Output: "libName.map"
def mapper(con,rawInputs,mode):
    
    ## For all libs - One by one
    for aninput in rawInputs:
        print('\nInput:',(aninput))
        lib,ext,nthread,maxTagLen = aninput
        
        # Resolve index path #################################
        cur = con.cursor()
        cur.execute("SELECT bowtie_index_path FROM master.genome_db WHERE genome_db like '%s'" % (genomeDB)) ##bowtie_index_path
        genoIndexPath = cur.fetchall()
        
        genoIndex = genoIndexPath[0][0].replace('$ALLDATA', '/alldata') ### Index file
        print ('Genomic index being used for mapping: %s\n'% (genoIndex))
        #genoIndex = 'ASPARAGUS_UGA1_genome' ## Test
        
        ### Prepare ###########################################
        inFile = '%s.%s' % (lib,ext)
        print ('Processing %s for mapping to genome' % (inFile))
        fastaFile = tagCount2FASTA(inFile,'N') ## Unique reads to FASTA format 

        mapFile = ('./%s.map' % (lib))
        print(genoIndex,inFile,fastaFile,mapFile)
        
        ## Map to index ##########################################
        print ('Mapping %s processed file to genome' % (lib))
        nproc2 = str(nproc)
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])
        
        ## Bowtie 1 - So as to be compatible with current indexes
        retcode = subprocess.call(["bowtie","-f","-n 0","-p", nproc2,"-t" ,genoIndex, fastaFile, mapFile])
        
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\nBowtie mapping for %s complete' % (inFile) )
        else:
            print ("There is some problem with mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()
        
        ### Prepare lists for plotting
        mappedList,mappedAbunList = mappedStats(aninput,mode)
        print('Mapped Reads:',mappedList)
        print('Abundance of mapped:',mappedAbunList)
        allTagsList,allAbunList = tagCountStats(aninput,mode)
        print('\nAll Reads:',allTagsList)
        print('Abundance of all sizes:',allAbunList)
        
        
        #### Test
        ###mappedList  =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95075, 166790, 278740, 869086, 735439, 1515217, 7389751, 694494, 122211, 60005, 46023, 39329, 33565, 26818, 19973, 15328, 11599, 842, 648, 579, 653, 1280, 1217, 1219, 1277, 955, 856, 749, 1268, 960, 766, 708, 1983, 28293, 0, 0, 0, 0, 0, 0, 0, 0]
        ###mappedAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 594218, 805020, 1025890, 5581017, 4444132, 4992476, 20590608, 1714861, 805331, 732898, 595526, 476446, 392119, 299055, 216764, 151625, 91236, 1205, 851, 862, 1039, 3765, 3022, 2628, 3144, 1791, 1727, 1300, 2696, 1905, 2014, 1783, 9453, 856855, 0, 0, 0, 0, 0, 0, 0, 0]
        ###
        ###allTagsList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126163, 220695, 370421, 1103866, 954861, 1886032, 9010585, 1012559, 274245, 140174, 105363, 91338, 82506, 83528, 54283, 56415, 56744, 16843, 20320, 25321, 21814, 41079, 29515, 27635, 23628, 26212, 17507, 13588, 18378, 10826, 8296, 10611, 28215, 483608, 0, 0, 0, 0, 0, 0, 0, 0]
        ###allAbunList =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 660160, 944285, 1217495, 6338895, 5015388, 5567509, 23419384, 2145615, 1029584, 858822, 709178, 672526, 658077, 416777, 348543, 248074, 173785, 21838, 23572, 28526, 26472, 77881, 53331, 41566, 36627, 33736, 22249, 17419, 24912, 13704, 10567, 14170, 42449, 1689522, 0, 0, 0, 0, 0, 0, 0, 0]
        ###
        ###mappedList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95075, 166790, 278740, 869086, 735439, 1515217, 7389751, 694494, 122211, 60005, 46023, 39329, 33565, 26818, 86042]
        ###mappedAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 594218, 805020, 1025890, 5581017, 4444132, 4992476, 20590608, 1714861, 805331, 732898, 595526, 476446, 392119, 299055, 1803027]
        ###
        ###allTagsList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126163, 220695, 370421, 1103866, 954861, 1886032, 9010585, 1012559, 274245, 140174, 105363, 91338, 82506, 83528, 644905]
        ###allAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 660160, 944285, 1217495, 6338895, 5015388, 5567509, 23419384, 2145615, 1029584, 858822, 709178, 672526, 658077, 416777, 2948943]
        ###
        ###mappedList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95, 166, 278, 869, 735, 1515, 7389, 694, 122, 600, 460, 39, 33, 26, 86]
        ###mappedAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 594, 805, 1025, 5581, 4444, 4992, 20590, 1714, 805, 732, 595, 476, 392, 299, 180]
        ###
        ###allTagsList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126, 220, 370, 1103, 954, 1886, 9010, 1012, 274, 140, 105, 913, 825, 835, 644]
        ###allAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 660, 944, 1217, 6338, 5015, 5567, 23419, 2145, 1029, 858, 709, 672, 658, 416, 294]

        ## Plot
        charts(lib,mappedList,mappedAbunList,allTagsList,allAbunList,mode) ## Mode 1 - Preprocess graphs 2: processed files graphs
        
    
    return None

def indexBuilder(genoFile):
    print ("\n**Deleting old index 'folder' !!!!!!!!!!!***\n")
    print('If its a mistake cancel now by pressing ctrl+D and continue from index step by turning off earlier steps- You have 30 seconds')
    time.sleep(30)
    
    shutil.rmtree('./index', ignore_errors=True)

    os.mkdir('./index')
    genoIndex = './index/%s' % (genoFile.rpartition('/')[-1].rpartition('.')[0]) ## Can be merged with genoIndex from earlier part if we use bowtie2 earlier
    print('**Creating index of cDNA/genomic sequences:%s**\n' % (genoIndex))
    retcode = subprocess.call(["bowtie2-build", genoFile, genoIndex])
    return genoIndex

def splicedMapper(rawInputs,genoIndex,gffFile):
    
    for aninput in rawInputs:
        print('\nInput:',(aninput))
        lib,ext,nthread = aninput
        
        ### Prepare ###########################################
        inFile = '%s.%s' % (lib,ext)
        print ('Mapping %s for to genome for RNA-Seq analysis' % (inFile))
        
        mapFile = ('./%s.tophat.map' % (lib))
        print(genoIndex,inFile,mapFile)
        
        ## Map to index ##########################################
        print ('Mapping %s processed file to genome' % (lib))
        nproc2 = str(nproc)
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])
        
        ## Bowtie 1 - So as to be compatible with current indexes
        retcode = subprocess.call([tophat,"-p", nproc2,"-G", gffFile, "-o", mapFile, genoIndex, inFile])
        
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\nTophat mapping for %s complete' % (inFile) )
        else:
            print ("There is some problem with Tophat mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()

def cuffLinks(rawInputs):
    
    print ('\n**Cufflink anlysis in progress**\n')
    
    assemblyOut = open('assemblies.txt','w')
    for aninput in rawInputs:
        
        print('\nInput:',(aninput))
        lib,ext,nthread = aninput
        
        ### Prepare ###########################################
        inFile = '%s.%s/accepted_hits.bam' % (lib,ext)
        print ('Cufflinking %s for to genome for RNA-Seq analysis' % (inFile))
        
        outDir = ('./%s.cufflink.out' % (lib))
        print(inFile,outDir)
        
        ## Map to index ##########################################
        print ('Mapping %s processed file to genome' % (lib))
        nproc2 = str(nproc)
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])
        
        ## Bowtie 1 - So as to be compatible with current indexes
        print ('****Binaries specified in user settings will be used: %s****' % (cufflinks))
        retcode = subprocess.call([cufflinks,"-p", nproc2,"-o", outDir, inFile])
        
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\nCuffLinks for %s complete' % (inFile) )
        else:
            print ("There is some problem with Cufflinks of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()
            
    return None

def cuffMergeQuant2(rawInputs,gffFile,genoFile):
    
    ### Cuffmerge
    print ('Creating a merged assembly\n')
    assemblyFile = 'assemblies.path'
    assemblyOut = open(assemblyFile,'w')
    for aninput in rawInputs:  
        lib,ext,nthreads = aninput
        assemblyPath = './%s.%s/transcripts.gtf\n' % (lib,ext)
        print(assemblyPath)
        assemblyOut.write(assemblyPath)
    
    assemblyOut.close()
    
    #print('\n**Running cuffMerge**\n')
    #print ('****Binaries specified in user settings will be used: %s****' % (cuffmerge))
    #nproc2 = str(nproc)
    ##print('%s,%s,%s,%s' % (gffFile,genoFile,nproc2,assemblyOut))
    #retcode = subprocess.call([cuffmerge, "-g", gffFile, "-s", genoFile,"-p", nproc2, assemblyFile])
    
    
    ### Cuffquants
    print('\n**Running cuffQuants**\n')
    
    mappedFiles = []

    for aninput in rawInputs[:2]:
        lib,ext,nthreads = aninput
        mappedPath = './%s.tophat.map/accepted_hits.bam ' % (lib)
        mappedFiles.append(mappedPath)
    
    mappedFilesPath = "%s" % (' '.join(mappedFiles))
    print ('%s' % mappedFilesPath)

    print ('\n**** Binaries specified in user settings will be used: %s ****\n' % (cuffquant))
    nproc2 = str(nproc)
    mergedAssembly = './merged_asm/merged.gtf' ## Merged annotatiions from cufflinks and original genome gtf file
    
    '''
    -M/--mask-file <mask.(gtf/gff)> Tells Cuffquant to ignore all reads that could have come from transcripts in this GTF file.
                                    We recommend including any annotated rRNA, mitochondrial transcripts other abundant transcripts you wish to ignore in your analysis in this file. 
    -u/--multi-read-correct	 Tells Cuffquant to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome.
    '''
    
    command = "%s -p %s -u %s %s" % (cuffquant,nproc2,mergedAssembly,mappedFilesPath) ## Had problem supplying multiple files seprated with space using normal switch based approach approach
    print(command)
    retcode =subprocess.call(command,shell=True)
        
    return None

def cuffMerge(rawInputs,gffFile,genoFile):
    
    '''
    Create a merged transcript
    '''
    
    print ('\n**Creating a merged assembly**\n')
    assemblyFile = 'assemblies.path'
    assemblyOut = open(assemblyFile,'w')
    for aninput in rawInputs:  
        lib,ext,nthreads = aninput
        assemblyPath = './%s.%s/transcripts.gtf\n' % (lib,ext)
        print(assemblyPath)
        assemblyOut.write(assemblyPath)
    
    assemblyOut.close()
    
    print('\n**Running cuffMerge**\n')
    print ('****Binaries specified in user settings will be used: %s****' % (cuffmerge))
    nproc2 = str(nproc)
    #print('%s,%s,%s,%s' % (gffFile,genoFile,nproc2,assemblyOut))
    print ('System:', cuffmerge, "-g", gffFile, "-s", genoFile,"-p", nproc2, assemblyFile)
    retcode = subprocess.call([cuffmerge, "-g", gffFile, "-s", genoFile,"-p", nproc2, assemblyFile])
    
    return None
    
def cuffQuant(aninput):
    '''
    Quatify the reads to merged assembly for every sample
    '''
    print('\n**Running cuffQuants**\n')   
    print ('\n**** Binaries specified in user settings will be used: %s ****\n' % (cuffquant))
    
    lib,ext,nthread = aninput
    inFile = '%s.%s/accepted_hits.bam' % (lib,ext)
    
    ## From cuffMerge
    mergedAssembly = './merged_asm/merged.gtf' ## Merged annotatiions from cufflinks and original genome gtf file

    outDir = '%s.cuffquant.out' % (lib)
    print(inFile,outDir)        
    retcode = subprocess.call([cuffquant,"-p", nthread,"-u", mergedAssembly, inFile, "-o", outDir])
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        print('\nCuffquants for %s complete' % (inFile) )
    else:
        print ("There is some problem with Cuffquants of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
        print ("Script exiting.......")
        sys.exit()
        
    return None

def cuffNorm(libs,reps,libType):
    
    print ('\n** Executing cuffNorms **')
        
    sampleDict = dict(zip(libs,reps)) ### Make a dictionary of samples and replicates
    print(sampleDict)
    repSet = set(reps) ### Get unique set of reps
    #repsCountDict = {x:reps.count(x) for x in reps} ## How many replicates of reach sample
    
    repList = [] ## List to hold replicates in a tuple and later join with space
    for rep in repSet:
        templist = [] ## Hold filenames for replicates
        print ('Replicate',rep)
        
        for lib in libs:
            #print ('Lib:',lib)
            aval = sampleDict[lib]
            if aval == rep:
                print ('Lib:',lib,'Rep:',rep)
                cxbFile = '%s.cuffquant.out/abundances.cxb' % (lib)
                templist.append(cxbFile)
        repList.append(templist)

    
    design = '' ## Will hold replicate wise samples for input to cuffNorm
    sampleList = []
    for samples in repList:
        replicates = ','.join(samples)
        sampleList.append(replicates)
    
    design = ' '.join(sampleList)
    print('\nThis is the sample-replicate design:',design)

    mergedAssembly = './merged_asm/merged.gtf' ## From cuffMerge
    outDir = 'Final.cuffNorm.out'
    nproc2 =  str(nproc)
    if libType == 0:
        command = '%s -p %s --library-type fr-unstranded %s %s -o %s' % (cuffnorm,nproc2,mergedAssembly,design,outDir)
    if libType == 1:
        command = '%s -p %s --library-type fr-firststrand %s %s -o %s' % (cuffnorm,nproc2,mergedAssembly,design,outDir)
    if libType == 2:
        command = '%s -p %s --library-type fr-secondstrand %s %s -o %s' % (cuffnorm,nproc2,mergedAssembly,design,outDir)
    print('System:', command)
    retcode = subprocess.call(command,shell = True)
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        print('\nCuffnorms for %s complete for replicate-sample design' % (design) )
    else:
        print ("There is some problem with Cuffnorms for replicate-sample design %s - Debug for reason" % (design))
        print ("Script exiting.......")
        sys.exit()
        
        
    return None

def compileResults():
    ## Combine 
    pass

############STANDARAD FUNCTIONS ##############

def PP(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(module, alist)

def PPBalance(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = round((nproc/int(nthread))+1) ## 1 added so as to avoid 0 processor being allocated in serial mode
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)
    
def ConnectToDB(server,infile):
    
    ##infile values are '0' when you dont want to uplaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE bioinfo_data.mir_page_results FIELDS TERMINATED BY ','")
    
    print ('\nTrying to connect to mySQL server on %s' % (server))
    # Try to connect to the database
    try:
        con=sql.connect(host= server, user='kakrana', passwd='livetheday')###local_infile = 1 not supported yet so a table has to be updated on row basis
        print ('Connection Established\n')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con

def tagCount2FASTA(inFile,Exprs):
    
    fh_in=open(inFile, 'r')
    outFile = '%s.fa' % (inFile.rpartition('.')[0])
    fh_out =open(outFile, 'w')
    tag_num = 1 ### For naming tags

    if Exprs=='Y':  ### Write as raw sequencing file with tag repeate dnumber of times it appears in tag_count 
        ##Write to file
        print('\nWriting expression file for %s tagcount file' % (inp_file_name))
        print('\n---PLEASE BE PATIENT---')
        
        for ent in fh_in:##All the entries of the library
            #if len(ent[0]) == 20:
            ent = ent.split('\t')
            tag_count = int(ent[1])
            for count in range(tag_count):##Number of times the tag_count file
                fh_out.write('>Tag%s\n%s\n' % (tag_num, ent[0]))
                tag_num += 1
                
    else: ##Convert tag count to FASTA
        for i in fh_in:
            ent = i.strip('\n').split('\t')
            #print(ent)
            fh_out.write('>Tag%s_%s\n%s\n' % (tag_num,ent[1],ent[0]))
            tag_num += 1
            
    fh_in.close()
    fh_out.close()
    
    return outFile

## Parse map file and collect statics for graph generation
def mappedStats(aninput,mode):
    
    #print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('\nCollecting statistics of matched reads for Lib:%s' % (lib))
    
    inFile = '%s.map' % (lib)
    fh_in = open(inFile,'r')
    mapFile = fh_in.read().split('\n')
    
    if mode == 1:
        mappedList = [0]*(maxReadLen) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        mappedAbunList = [0]*(maxReadLen) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    elif mode == 2:
        mappedList = [0]*(maxTagLen+1) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        mappedAbunList = [0]*(maxTagLen+1) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    else:
        print('\nThe mode selected for collecting mapped stats is not correct - Debug for reason')
    
    for anent in mapFile[:-1]: ## Last entry is empty due to split on newline
        ent = anent.split('\t')
        #print(ent)
        
        tagName,Abun = ent[0].split('_')
        #print(tagName,Abun,len(ent[4]))
        mappedList[len(ent[4])] += 1
        mappedAbunList[len(ent[4])] += int(Abun)
    
    return mappedList,mappedAbunList

## Get stats for all the reads from tagCount file
def tagCountStats(aninput,mode):
    
    #print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('\nCollecting statistics of total reads for Lib:%s' % (lib))
    
    inFile = '%s.processed.txt' % (lib)
    print(inFile)
    fh_in = open(inFile,'r')
    tagCountFile = fh_in.read().split('\n')
    
    if mode == 1:
        allTagsList = [0]*(maxReadLen) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        allAbunList = [0]*(maxReadLen) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    elif mode == 2:
        allTagsList = [0]*(maxTagLen+1) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        allAbunList = [0]*(maxTagLen+1) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    else:
        print('\nThe mode selected for collecting tag count stats is not correct - Debug for reason')

    
    for anent in tagCountFile[:-1]: ## Last entry is empty due to split on newline
        #print(anent)
        tagSeq,Abun = anent.split('\t')
        #print(tagSeq,Abun)
        
        allTagsList[len(tagSeq)] += 1
        allAbunList[len(tagSeq)] += int(Abun)
    
    #print('Total Tags',allTagsList,'\n','Total Abundance',allAbunList)
    
    return allTagsList,allAbunList

def charts(lib,mappedList,mappedAbunList,allTagsList,allAbunList,mode):
    
    ### Graphs of pre-processed files i.e trimmed files
    if mode == 1:
        
        print ("\n**Generating graphs for trimmed files**\n")
        ## Get all the tag sizes from disticnt tags list
        #print('alltagsList:', allTagsList)
        indexList = [i for i,x in enumerate(allTagsList) if x != 0]
        #print ('indexList:',indexList)
        minLen = min(indexList)
        maxLen = max(indexList)
        
        #### CHART-1: Distinct mapped vs distinct all
        plotFile = ('%s_distinct_dist_before_chop.png' % (lib))##Plot results file
        
        #bottomList = [i for i in mappedList if i > 0]
        #upList = [i for i in allTagsList if i > 0]
        bottomList = list(mappedList[minLen:maxLen+1])
        upList = list(allTagsList[minLen:maxLen+1]) ## Abundance of different sizes
        upList2 = [a - b for a, b in zip(upList, bottomList)] ## Abundance substracted from bottom List - to plot the remainder on top
        
        maxAbun = max(upList)
        #print (len(bottomList),len(upList),maxAbun)
        ybreak = 500000
        
        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1
        ind=np.arange(N)
        #print('np.arange',ind)
        width = 0.5
        
        ##plotting variables
        p1 = plt.bar(ind, bottomList, width, color = 'g',)
        p2 = plt.bar(ind, upList2, width, color = 'b', bottom=bottomList)
        
        plt.ylabel('Count of distinct tags mapped to genome (before chopping)', fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length', fontproperties=font_manager.FontProperties(size=8))
        plt.title('Distinct tags in %s library matched to %s genome before pre-processing' % (lib,genomeDB),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N),    np.arange(int(minLen),int(maxLen)+1),   rotation = 45,  fontproperties=font_manager.FontProperties(size=9))
        plt.yticks(np.arange(0,maxAbun+ybreak,ybreak),  fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))
        
        plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file
        
        
        
        ### CHART 2:  Mapped abundance vs total abundance ######
        plotFile2 = ('%s_abund_dist_before_chop.png' % (lib))##Plot results file
        #bottomListAll = [i for i in mappedAbunList if i > 0]
        #upListAll = [i for i in allAbunList if i > 0]
        bottomListAll = list(mappedAbunList[minLen:maxLen+1])
        upListAll = list(allAbunList[minLen:maxLen+1])
        upListAll2 = [a - b for a, b in zip(upListAll, bottomListAll)] ## Abundance substracted from bottom List - to plot the remainder on top
        
        maxAbun2 = max(upListAll)
        #print (upListAll,maxAbun2)
        ybreak2 = 500000
        
        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1
        #ind=np.arange(N)
        width = 0.5
        
        ##plotting variables
        p1 = plt.bar(ind, bottomListAll, width, color = 'm',)
        p2 = plt.bar(ind, upListAll2, width, color = 'c', bottom=bottomListAll)
        
        plt.ylabel('Total tags', fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length', fontproperties=font_manager.FontProperties(size=8))
        plt.title('Total tags in %s library matched to %s genome (before chopping)' % (lib,genomeDB),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N), np.arange(int(minLen),int(maxLen)+1),  rotation = 45,  fontproperties=font_manager.FontProperties(size=9) )
        plt.yticks(np.arange(0,maxAbun2+ybreak2,ybreak2),   fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))

        plt.savefig(plotFile2, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file
    
    ### Graphs of processed files
    elif mode == 2: ## Graph of final processed reads:
        
        print ("\n**Generating graphs for final processed files**\n")
        
        indexList = [i for i,x in enumerate(allTagsList) if x != 0]
        #print ('indexList:',indexList)
        minLen = min(indexList)
        maxLen = max(indexList)
        
        #### CHART-1: Distinct mapped vs distinct all
        plotFile = ('%s_distinct_dist_after_chop.png' % (lib))##Plot results file
        
        bottomList = list(mappedList[minLen:maxLen+1])
        upList = list(allTagsList[minLen:maxLen+1]) ## Abundance of different sizes
        upList2 = [a - b for a, b in zip(upList, bottomList)] ## Abundance substracted from bottom List - to plot the remainder on top
        
        
        maxAbun = max(upList)
        #print (upList,maxAbun)
        ybreak = 500000
        
        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1 
        ind=np.arange(N)
        #print('np.arange',ind)
        width = 0.5
        
        ##plotting variables
        p1 = plt.bar(ind, bottomList, width, color = 'g',)
        p2 = plt.bar(ind, upList2, width, color = 'b', bottom=bottomList)  
        
        plt.ylabel('Count of distinct tags mapped to genome', fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length', fontproperties=font_manager.FontProperties(size=8))
        plt.title('Distinct tags in %s library matched to %s genome after processing' % (lib,genomeDB),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N),    np.arange(int(minLen),int(maxLen)+1),   rotation = 45,  fontproperties=font_manager.FontProperties(size=9))
        plt.yticks(np.arange(0,maxAbun+ybreak,ybreak),  fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))
        
        plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file
        
        
        ### Chart 2:  All genomic reads
        plotFile2 = ('%s_abund_dist_after_chop.png' % (lib))##Plot results file
        bottomListAll = list(mappedAbunList[minLen:maxLen+1])
        upListAll = list(allAbunList[minLen:maxLen+1])
        upListAll2 = [a - b for a, b in zip(upListAll, bottomListAll)] ## Abundance substracted from bottom List - to plot the remainder on top
        
        maxAbun2 = max(upListAll)
        #print (upListAll,maxAbun2)
        ybreak2 = 500000
        
        ## Different sizes to be plotted
        N = int(maxLen) - int(minLen) +1 
        ind=np.arange(N)
        width = 0.5
        
        ##plotting variables
        p1 = plt.bar(ind, bottomListAll, width, color = 'm',)
        p2 = plt.bar(ind, upListAll2, width, color = 'c', bottom=bottomListAll)
        
        plt.ylabel('Total tags', fontproperties=font_manager.FontProperties(size=8))
        plt.xlabel('Tag Length', fontproperties=font_manager.FontProperties(size=8))
        plt.title('Total tags in %s library matched to %s genome after processing' % (lib,genomeDB),fontproperties=font_manager.FontProperties(size=9))
        plt.xticks(np.arange(N), np.arange(int(minLen),int(maxLen)+1),  rotation = 45,  fontproperties=font_manager.FontProperties(size=9) )
        plt.yticks(np.arange(0,maxAbun2+ybreak2,ybreak2),   fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))

        plt.savefig(plotFile2, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file

    else:
        print('\nThe mode selected for generating graph is not correct - Debug for reason')
        
    
    return None

############## MAIN #########################
def main(libs):
    
    global con
    con = ConnectToDB(dataServer,0)
    
    start = time.time() 
    runLog = open('%s_run.log' % (datetime.datetime.now().strftime("%m_%d_%H_%M")), 'w')
    
    #### 0. Initialize Input Register ###################################################
    register =[]
    print('\n\n**Total %s libraries provided for pre-processing**\n' % (len(libs)))
    for i in libs:
        cur = con.cursor()
        cur.execute("SELECT raw_path,raw_file_format,adapter_5p,adapter_3p,minimum_tag_length,maximum_tag_length FROM master.library_info WHERE lib_id like '%s'" % (i)) ##bowtie_index_path
        fileInfo = cur.fetchall()
        print ('Library:', (fileInfo[0]))
        
        filePath = fileInfo[0][0].replace('$ALLDATA', '/alldata') 
        ext = fileInfo[0][1]
        adp_5p = fileInfo[0][2]
        adp_3p = fileInfo[0][3]
        minTagLen = fileInfo[0][4]
        maxTagLen = fileInfo[0][5]
        register.append((str(i),ext,str(nthread),str(filePath),adp_5p,adp_3p,minTagLen,maxTagLen)) ## Lib, ext, nthread, rawa file path, adapter 5p, adapter 3p, min tag len, max tag len
    
    #### 1. QC Raw files ################################################################
    rawInputs = [(i[0],i[1],i[2],i[3]) for i in register]
    if QCheckStep == 1:
        print('\n**Quality check of the raw files will be performed now**\n')
        PPBalance(QCheck,rawInputs)
    else:
        print ('\n**Quality check of the raw files will be skipped as selected**\n')
    
    
    #### 2. TRIM RAW files #############################################################
    rawInputs = [(i[0],i[1],i[2],i[3],i[4],i[5],i[6]) for i in register]
    if trimLibsStep == 1:

        print('\n**Trimming of the libraries will be performed now**\n')
        #for i in rawInputs:
        #    trimLibs(i)
        PPBalance(trimLibs,rawInputs)
    else:
        print('\n**Trimming of the libraries will be skipped as selected**\n')
        pass
    
    #### 3. Prepare graphs of trimmed files #############################################
    rawInputs = [(i[0],'trimmed.fastq',i[2]) for i in register]
    if preProGraphsStep == 1:
        print('\n**Converting trimmed files to tag count format for quality graphs**\n')
        PP(fastQ2Count,rawInputs)
        con = ConnectToDB(dataServer,0)
        maps = mapper(con,rawInputs,1)
        
        ### Delete tag count, fasta files and map files to make way for real processed files
        print ("**Cleaning temp files**")
        garbage = [file for file in os.listdir('./') if file.endswith (('.map','processed.txt','processed.fa'))]
        for file in garbage:
            print("Deleting %s" % (file))
            os.remove(file)
    else:
        pass

    #### 3. Chop trimmed files ###########################################################
    rawInputs = [(i[0],'trimmed.fastq',i[2],i[7]) for i in register] ## Lib, extension,threads, max tag len
    if chopLibsStep == 1:
        print('\n**Chopping of the libraries will be performed now**\n')
        PPBalance(chopLibs,rawInputs)
    else:
        print('\n**Chopping of the libraries will be skipped as selected**\n')
        pass


    #### 4. Convert  chopped files to tagcount format ###################################
    rawInputs = [(i[0],'chopped.trimmed.fastq',i[2]) for i in register]
    if fastQ2CountStep == 1:
        print('\n**Converting processed files to tag count format**\n')
        PP(fastQ2Count,rawInputs)
    else:
        pass
    
    
    #### 5. Map tagcount to genome index ################################################
    rawInputs = [(i[0],'processed.txt',i[2],i[7]) for i in register]
    if mapperStep == 1:
        print('\n**Mapping chopped files to the genome index**\n')
        con = ConnectToDB(dataServer,0)
        mapper(con,rawInputs,2)
    else:
        print("\n**Mapping Step skipped as selected - What are you intentions?\n")
        pass

    #### 6. Clean up ###################################################################
    if cleanupStep == 1:
        print ("**Cleaning temp files last time**")
        garbage = [file for file in os.listdir('./') if file.endswith (('.map','trim.log','processed.fa','.zip'))] ## Excluded 'chopped.trimmed.fastq' as used by RNA Runner
        for file in garbage:
            print("Deleting %s" % (file))
            os.remove(file)
    else:
        pass
    
    end = time.time()
    runLog.write('Complete run time is %s' % (round(end-start,2)))
    runLog.close()
    
    print ('Complete run time is %s seconds' % (round(end-start,2)))
    
    ####7. RNA-SEQ runner ######################################
    if indexBuilderStep == 1:
        genoIndex = indexBuilder(genoFile)
    else:
        genoIndex = './index/%s' % (genoFile.rpartition('/')[-1].rpartition('.')[0])
        pass
    
    rawInputs = [(i[0],'chopped.trimmed.fastq',i[2]) for i in register]
    
    if spliceMapperStep == 1:
        print('Index being used: %s' % (genoIndex))
        splicedMapper(rawInputs,genoIndex,gffFile)
    
    if cuffLinksStep == 1:
        rawInputs = [(i[0],'tophat.map',i[2]) for i in register] ## 'tophat.map' is a folder
        assemblyOut = cuffLinks(rawInputs)
    
    if cuffMergeStep == 1:
        rawInputs = [(i[0],'cufflink.out',i[2]) for i in register]
        finalAssembly = cuffMerge(rawInputs,gffFile,genoFile)
    
    if cuffQuantStep == 1:
        rawInputs = [(i[0],'tophat.map',i[2]) for i in register]
        PPBalance(cuffQuant,rawInputs)
    
    if cuffNormStep == 1:
        rawInputs = [(i[0],'cuffquant.out',i[2]) for i in register]
        cuffNorm(libs,reps,libType)

################ Execute #######################
if __name__ == '__main__':
    
    #### Assign Cores
    if numProc == 0:
        nproc = int(multiprocessing.cpu_count()*0.90)
    else:
        nproc = int(numProc)
    
    #### Execute modules

    main(libs)
    print("\n\n-------Script finished sucessfully - CHEERS!!!! - You owe a beer !_! to Atul--------\n")
    
    sys.exit()


### CHANGELOG ---------------------------------------------------
## Version-1 preprocessing
##Screw losers i.e. Sa ta- I will write myself

##v2 -> v05
## First working copy of the script

## v5 -> v06
## Automatic resolve path for raw fasta file - so no need to copy in working directory
## Get tha adapter sequences for each library from master DB
## Get the max and min tag len from master DB
## Fixed maxTagLen issue - it was being swapped with minTagLen

## v06 -> v07
## Added RNA seq analysis till cuffmerge

## v07 -> v01 of seqAnalyze
## Completed the pipeline with most default parameters

## v07 -> Rocketv02
## Fixed toolPath for chopping module

### Future fixes -----------------------------------------------
## Which quality scores used for sequencing: phred33 or phred66
## Adpaters need to selected or specified at begning of preprocessing
## Xlabels for abundance graph needs to be fixed
## Bring back the local version for publication

