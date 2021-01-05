#!/usr/bin/python3

## version: v4.0b2

## Autor: kakrana@udel.edu
## "rocket" is a script for ab-intio assembly, it pre-process and runs Tophat/Bowtie based 
## assembly process - It can be used for assembly and abudannce quantification
## Requires fastQC, trimmomatic and tally in $PATH variable

## Run : python3 ScriptName.py

## ENVIRONMENT #################################
import os, yaml
from os.path import expanduser
HOME    = expanduser("~")
CONFIG  = yaml.load(open('rocket.yaml', 'r'), Loader=yaml.FullLoader)

## IMPORTS #####################################
import shutil
import datetime
import operator
import subprocess
import matplotlib
import numpy as np
import multiprocessing
import itertools as it
from multiprocessing import Process, Queue, Pool
import matplotlib.font_manager as font_manager
import sys,os,re,time,timeit,csv,glob,string
matplotlib.use('Agg')
import matplotlib.pyplot as plt

## REQUIRED INFO ##############################
## 1. Lib preparation protocol for value to "libType", "seqType", "adapterFileSE" or "adapterFilePE" parameters
## 2. splicesites.txt generated using HISAT2 on GTF file (Hisat oes not uses GTF files directly)

## PRE-PROCESSING SETTINGS ####################
Local           = CONFIG['dev']['local']        ## 1: Local 0: Remote | NOTE: Local Analysis: Requires - maxLen,minLen,maxReadLen and adpaters.fa and libraries
                                                ## Requires SampleInfo file with sampleNum, Filename/Lib code, reps, group
genomeDB        = CONFIG['dev']['genomeDB']     ## [Server mode]For Bowtie index path used for mapping for graphs
referenceGTF    = CONFIG['dev']['referenceGTF'] ## [optional] T: True - use reference gtf file for merging assembling and annotation | F: Do not use GTF file and report transcripts based on transcriptome. Process GTF file for hisat: hisat2_extract_splice_sites.py genes.gtf > splicesites.txt
                                                ## http://plants.ensembl.org/info/website/ftp/index.html OR USE gffread my.gff3 -T -o my.gtf (GFFREAD IS PART OF CUFFLINKS/TUXEDO)
## ROCKET #####################################
userMinTagLen   = CONFIG['seq']['userMinTagLen']       ## [server] Used if 'hardMinTagLen' is ON
userMaxTagLen   = CONFIG['seq']['userMaxTagLen']       ## [server] Used if 'hardMinTagLen' is ON
libType         = CONFIG['seq']['libType']      ## [mandatory] From HISAT2 manual choose mode for strand-specificity (https://www.biostars.org/p/262027/)
                                                ## 0 = F or reads corresponds to transcript
                                                ## 1 = R or reads correspond to reverse complemented counterpart of a transcript
                                                ## 2 = RF or fr-firststrand (dUTP method, NSR, NNSR, Illumina Tru-Seq stranded protocol)
                                                ## 3 = FR or fr-secondstrand (RNA linkers ligated)
seqType         = CONFIG['seq']['seqType']      ## [mandatory] 0: Single End; 1:Paired end (requires splitted reads - see fastq dump --split-reads for lib/or custom script)

gtfFile         = CONFIG['user']['gtfFile']     ## GTF file for stringtie
ssFile          = CONFIG['user']['ssFile']      ## Splice site file for HISAT2 (see referenceGTF setting for info)
genoFile        = CONFIG['user']['genoFile']    ## Reference genome file
genoIndexHS     = CONFIG['user']['genoIndex']   ## If index is not in $ALLDATA i.e. local analysis, then specify bowtie1 index here for pre-processing. For Seq-analysis a Bowtie2 index  will be made using 'indexBuilderStep'
genoIndexPrePro = CONFIG['user']['genoIndexPrePro'] ## Bowtie1 index used for generating mapping charts 
sampleInfo      = CONFIG['user']['sampleInfo']  ## [mandatory] Tab delimted file with three mandatory columns - num (sample numbers), id (filename,library id), rep (same number if replicates
                                                ## And one optional columns group (sample grouped for edgeR analysis). See end of code for format.
groupBy         = CONFIG['user']['groupBy']     ## [mandatory]   R: Group Samples by replicates, G: By user specified groups in sampleInfo 'group' column

## PRE_PROCESSING - OPTIONAL STEPS ###############
QCheckStep      = CONFIG['steps']['QCheckStep']  ## Optional -Performs preliminary QC

## PRE_PROCESSING - REQUIRED STEPS [Value: 0/1] ##############
trimLibsStep     = CONFIG['steps']['trimLibsStep']      ## Trim fastq files
preProGraphsStep = CONFIG['steps']['preProGraphsStep']  ## Generates before chopping graphs
chopLibsStep     = CONFIG['steps']['chopLibsStep']      ## Chops adapter trimmed files
fastQ2CountStep  = CONFIG['steps']['fastQ2CountStep']   ## Converts chopped to tag count
mapperStep       = CONFIG['steps']['mapperStep']        ## Maps final chopped files and generates graphs
summaryFileStep  = CONFIG['steps']['summaryFileStep']   ## Generates mapped summary - Never tested in Rocket, actually imported from prepro
cleanupStep      = CONFIG['steps']['cleanupStep']       ## Final cleanup

## SEQ-ANALYSIS - REQUIRED STEPS [Value: 0/1] ##############
indexBuilderStep = CONFIG['steps']['indexBuilderStep']  ## Build index for all the mappings
spliceMapperStep = CONFIG['steps']['spliceMapperStep']  ## HiSat2 Mapping
stringTieStep    = CONFIG['steps']['stringTieStep']     ## StringTie assemblies for all transcripts
stringMergeStep  = CONFIG['steps']['stringTieStep']     ## Merge GTFs to single assembly
stringCountStep  = CONFIG['steps']['stringCountStep']   ## StringTie assemblies for quantification
stringQuantStep  = CONFIG['steps']['stringQuantStep']   ## Generate counts table from all libraries

## ADVANCED SETTINGS #######################
unpairDel       = CONFIG['dev']['unpairDel']            ## [Only for paired end analysis] 0: Retain unpaired read files after trimming 1: Delete these files
#maxfrags        = CONFIG['dev']['stringQuantStep']     ##  Maximum fragments allowed in a bundle before skipping [ default: 500000 ]

ncores          = CONFIG['user']['numProc']             ## [developer]  Coarse grain PP [0: Maximize parallel processing | [1-64]: Number of Cores]
nthread         = CONFIG['dev']['nthread']              ## [developer]  Fine grain PP
maxReadLen      = CONFIG['dev']['maxReadLen']           ## [developer]  Max allowed unchopped read length for graph generation

masterDB        = 'master'              ## [server]
dataServer      = 'raichu.ddpsc.org'    ## [server]
hardMinTagLen   = CONFIG['dev']['hardMinTagLen']       ## [server] Override Chopping values from server 

## TOOL/FILE PATH ###################################
tally           = f'{HOME}/tools/tally/tally'                               ## [mandatory]
fastqc          = f'{HOME}/tools/FastQC/fastqc'                             ## [mandatory]
hisat2          = f'{HOME}/tools/hisat2-2.1.0/hisat2'                       ## [mandatory]
samtools        = f'{HOME}/tools/samtools-1.11/bin/samtools'                ## [mandatory] Needs compilation
bowtie          = f'{HOME}/tools/bowtie-1.3.0-linux-x86_64/bowtie'          ## [mandatory]
adapterFileSE   = f'{HOME}/tools/Trimmomatic-0.39/adapters/TruSeq-SE.fa'    ## [mandatory] Sequence adapter file in FASTA format - Trimmomatic has files for different kits - Copy from there
adapterFilePE   = f'{HOME}/tools/Trimmomatic-0.39/adapters/TruSeq-PE.fa'    ## [mandatory] You can try merged file with SE and PE adapter too, it gave best results to Atul, because adapters are sometimes not in pairs
stringtie       = f'{HOME}/tools/stringtie-1.3.5.Linux_x86_64/stringtie'    ## [mandatory]
prepDE          = f'{HOME}/tools/stringtie-2.1.4.Linux_x86_64/prepDE.py3'   ## [mandatory]
trimmomatic     = f'{HOME}/tools/Trimmomatic-0.39/trimmomatic-0.39.jar'     ## [mandatory]


## PREPROCESS #####################
def checkHost(allowedHost):
    print ("\n           ---> Checking HostName <----      ")
    f = subprocess.Popen("hostname", stdout=subprocess.PIPE,shell= True)
    output,err = f.communicate()
    #print (output.decode("ascii"))
    
    host = output.decode("ascii")
    print ('Current Host: ',host)
    
    ## DO not turn OFF this 'for' loop as that given an error while matching current host with allowedHost - Reason Unknown
    print ('Allowed Hosts:')
    for host in allowedHost:
        print (host)
    
    
    if str(host) in allowedHost:
        pass
    
    else:
        print ("** Violation of task allocation scheme ** ")
        print("** Data buiding and/or preprocessing is not allowed on this server **")
        print("          --> script will exit now <---   \n")
        sys.exit()

def readSet(setFile):
    print ("\n######## User Settings #############")
    
    fh_in = open("prepro.set", 'r')
    setFile = fh_in.readlines()
    
    for line in setFile:
        if line: ## Not empty
            if line.startswith('@'):
                line = line.strip().split('<')
                param,value = line[0].split('=')
                print(param,value)
                
                ##Extract values
                
                if param.strip() == '@genomeDB':
                    global genomeDB
                    genomeDB = value.replace('"','').strip()
                    #print("User input Genome: ",genomeDB)
                
                elif param.strip() == '@libs':
                    global libs
                    libs = list(map(int,value.strip().split(',')))
                    print('User Input Libs:',libs)
                
                elif param.strip() == '@QCheckStep':
                    global QCheckStep
                    QCheckStep = int(value.strip())
                    #print('User Input QCheckStep:',QCheckStep)
                
                elif param.strip() == '@preProGraphsStep':
                    global preProGraphsStep
                    preProGraphsStep = int(value.strip())
                    #print('User Input preProGraphsStep:',preProGraphsStep)
                
                elif param.strip() == '@trimLibsStep':
                    global trimLibsStep
                    trimLibsStep = int(value.strip())
                    #print('User Input preProGraphsStep:',trimLibsStep)
                
                elif param.strip() == '@chopLibsStep':
                    global chopLibsStep
                    chopLibsStep = int(value.strip())
                    #print('User Input preProGraphsStep:',chopLibsStep)
                
                elif param.strip() == '@fastQ2CountStep':
                    global fastQ2CountStep
                    fastQ2CountStep = int(value.strip())
                    #print('User Input preProGraphsStep:',fastQ2CountStep)
                
                elif param.strip() == '@mapperStep':
                    global mapperStep
                    mapperStep = int(value.strip())
                    #print('User Input preProGraphsStep:',mapperStep)

                elif param.strip() == '@summaryFileStep':
                    global summaryFileStep
                    summaryFileStep = int(value.strip())
                    #print('User Input preProGraphsStep:',mapperStep)
                
                elif param.strip() == '@cleanupStep':
                    global cleanupStep
                    cleanupStep = int(value.strip())
                    #print('User Input preProGraphsStep:',cleanupStep)
                
                
            else:
                #print("Missed line:",line)
                pass
    print('####################################')
    
    return genomeDB,libs

def QCheck(aninput):
    '''
    Output: "libName_fastqc"  
    '''
    print(aninput)
    lib,ext,nthread,infile = aninput
    print('****Checking quality of %s library****' % (lib))    
    # toolPath = "%s/tools/FastQC/fastqc" % (os.getenv('HOME'))

    outDir = "%s" % (os.getcwd())
    print(outDir)
    x = subprocess.Popen("%s --outdir=%s %s" % (fastqc,outDir,infile),shell=True)
    x.communicate() # now wait
    ## retcode2 = subprocess.call([toolPath,"--outdir"infile])

    ## Change result folder and file name to library name
    # existFolder = "%s_fastqc" % infile.split("/")[-1].split(".")[0]
    # existFile =  "%s_fastqc.zip" % infile.split("/")[-1].split(".")[0]

    # newFolder = "%s"  % lib
    # newFile =  "%s.zip" % lib
    # # print(newFile,newFolder)

    # ## Remove folder from earlier run and then rename
    # if os.path.isdir(newFolder):
    #     print("Existing FASTQC report will be deleted for lib:%s" % (lib))
    #     shutil.rmtree(newFolder, ignore_errors=True)
    #     os.rename(existFolder,newFolder)
    # else:
    #     os.rename(existFolder,newFolder)

    # ## Remove files from earlier run and then rename
    # if os.path.isfile(newFile):
    #     print("Existing FASTQC zipped report will be deleted for lib:%s" % (lib))
    #     os.remove(newFile)
    #     os.rename(existFile,newFile)
    # else:
    #     os.rename(existFile,newFile)

    return None

def trimLibs(aninput):
    '''
    Output: "libName.trimmed.fastq"
    '''
    print(aninput)
    lib,ext,nthread,infile,adp_5p,adp_3p,minTagLen = aninput
    # print (aninput)
    print('\n****Trimming %s library with min length %s****' % (lib,minTagLen))
    # toolPath = "%s/tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    
    #### LOCAL ##########################
    if Local == 1: 

        ## Single End ###################
        if seqType == 0:
            adapterFile = adp_5p ## For local analysis adp_3p is actually the adapterfile - see main
            trimmedFile = '%s.trimmed.%s' % (lib,ext) ## Output
            
            trimLog     = '%s.trim.log' % (lib) ## Output
            trimSumm    = '%s.trim.summ.txt' % (lib)
            
            retcode = subprocess.call(["java", "-jar", trimmomatic, "SE", "-phred33","-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:%s:2:30:10" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (minTagLen)])
            if retcode == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n****Trimming for %s complete****' % (infile) )
            else:
                print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
                sys.exit()

        ## Paired End ##################
        elif seqType == 1:
            adapterFile = adp_3p ## For local analysis adp_3p is actually the adapterfile - see main
            trimmedFileP1 = '%s.pair_1.trimmed.%s' % (lib,ext) ## Output
            trimmedFileP2 = '%s.pair_2.trimmed.%s' % (lib,ext) ## Output
            trimmedFileU1 = '%s.unpair_1.trimmed.%s' % (lib,ext) ## Output
            trimmedFileU2 = '%s.unpair_2.trimmed.%s' % (lib,ext) ## Output

            infile1 = "%s_1.%s" % (lib,ext) ## In local mode file are pre-named to have these suffix
            infile2 = "%s_2.%s" % (lib,ext) ## In local mode file are pre-named to have these suffix

            trimLog = '%s.trim.log' % (lib) ## Not used because slows down massively
            trimSumm = '%s.trim.summ.txt' % (lib) ## Used
            
            retcode = subprocess.call(["java", "-jar", trimmomatic, "PE", "-phred33", "-threads", nthread, infile1, infile2, trimmedFileP1,trimmedFileU1,trimmedFileP2,trimmedFileU2, "ILLUMINACLIP:%s:2:30:10:8:TRUE" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:%s" % (minTagLen)])
            if retcode == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n****Trimming for %s complete****' % (infile) )
            
            else:
                print('Something wrong happened while trimming library: %s - - Debug for reason' % (lib))
                sys.exit()

        ## SE or Paired End?
        else:
            print("Please input seqType - SE or PE")
            sys.exit()


    ## SERVER BASED ##########################
    elif Local == 0:
        adapter = open('%s_adapter.fa' % (lib), 'w')
        adapter.write('>adapter_5p\n%s\n>adapter_3p\n%s' % (adp_5p,adp_3p))
        adapter.close()

        ## Single End ##############
        if seqType == 0:
            trimmedFile = '%s.trimmed.%s' % (lib,ext) ## Output
            path,splt,filename = infile.rpartition("/")
            trimLog = '%s.trim.log' % (filename) ## Output
            retcode = subprocess.call(["java", "-jar", trimmomatic, "SE", "-phred33", "-trimlog",trimLog,"-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:./%s_adapter.fa:2:30:10" % (lib), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (minTagLen)])
            if retcode == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n****Trimming for %s complete****' % (infile) )
            
            else:
                print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
                sys.exit()

        ## Paired End ##################
        elif seqType == 1: 
            trimmedFileP1 = '%s.pair_1.trimmed.%s' % (lib,ext) ## Output
            trimmedFileP2 = '%s.pair_2.trimmed.%s' % (lib,ext) ## Output
            trimmedFileU1 = '%s.unpair_1.trimmed.%s' % (lib,ext) ## Output
            trimmedFileU2 = '%s.unpair_2.trimmed.%s' % (lib,ext) ## Output

            ## Regenerate names of input paired end files
            path,splt,filename = infile.rpartition("/")
            infile1 = "%s/%s_1.%s" % (path,filename.rpartition("_")[0],ext)
            infile2 = "%s/%s_2.%s" % (path,filename.rpartition("_")[0],ext)

            trimLog = '%s.trim.log' % (filename) ## Output
            
            retcode = subprocess.call(["java", "-jar", trimmomatic, "PE", "-phred33", "-trimlog",trimLog,"-threads", nthread, infile1, infile2, trimmedFileP1,trimmedFileU1,trimmedFileP2,trimmedFileU2, "ILLUMINACLIP:./%s_adapter.fa:2:30:10:8:TRUE" % (lib), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (minTagLen)])
            if retcode == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n****Trimming for %s complete****' % (infile) )
            else:
                print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
                sys.exit()
        
        ## SE or Paired End?
        else:
            print("Please input seqType - SE or PE")
            sys.exit()
    
    ## Local or Server mode ?       
    else:
        print("Please choose correct value for Local variable: Local = [0/1]")

    ## CLEANUP ##
    if unpairDel == 1:
        garbage = [afile for afile in os.listdir('./') if afile.endswith (('unpair_1.trimmed.fastq','unpair_2.trimmed.fastq'))] ## Excluded-'chopped.trimmed.fastq' as used by RNA Runner
        for afile in garbage:
            if os.path.isfile(afile): ## Check to see its a file from bowtie and not tophat mapped folder - Untested
                print("Deleting %s" % (afile))
                os.remove(afile)
            else:
                print("Skiping cleanup, as its a directory %s" % (afile))

    ### Make plot
    #charts(mappedList,mappedAbunList,allTagsList,allAbunList,mode)
    
    return None

def chopLibs(aninput):
    ''' 
    Reverse read set of paired end lib is chopped from right end (5' in actuality) as done by using -
    if we do chopping with trimming using PE mode it is still chopped the same way 
    Output: "libName.chopped.fastq" 
    '''

    print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('****Chopping %s library to length %s****' % (lib,maxTagLen))

    trimmedInFile  = '%s.%s' % (lib,ext)
    choppedOutFile = '%s.chopped.%s' % (lib,ext)
    print("\n")
    
    if seqType == 0:
        retcode = subprocess.call(["java", "-jar", trimmomatic, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (maxTagLen)])
    elif seqType == 1:
        retcode = subprocess.call(["java", "-jar", trimmomatic, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (maxTagLen)])
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Chopping for %s complete****' % (trimmedInFile) )
    else:
        print('Something wrong happened while chopping library: %s - - Debug for reason' % (lib))
        sys.exit()
    
    return None

def fastQ2Count(aninput):
    '''
    Output: "libName.processed.txt"
    '''
    print(aninput)
    lib,ext,nthread = aninput
    print('\n****Converting %s.%s file to tag count****\n' % (lib,ext))
    infile = '%s.%s' % (lib,ext)
    outfile = '%s.%s.processed.txt' % (lib,ext.replace(".fastq",""))
    print("This is outfile:%s" % (outfile))
    # sys.exit()
    retcode = subprocess.call([tally, "-i", infile, "-o", outfile, "--nozip", "-format","%R%t%X%n"])
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Conversion to tag count format for %s complete****' % (infile) )
    else:
        print("Something wrong happened while converting to %s library tag count - Debug for reason" % (lib))
        sys.exit()

def mapper(rawInputs,mode):
    '''
    Output: "libName.map"
    '''
    
    ## For all libs - One by one
    for aninput in rawInputs:
        print('\nInput:',(aninput))
        lib,ext,nthread,maxTagLen,genoIndexPrePro = aninput
        
        # Resolve index path #################################
        if Local == 0:
            cur = con.cursor()
            cur.execute("SELECT bowtie_index_path FROM master.genome_db WHERE genome_db like '%s'" % (genomeDB)) ##bowtie_index_path
            genoIndexPath   = cur.fetchall()
            genoIndexPrePro = genoIndexPath[0][0].replace('$ALLDATA', '/alldata') ### Index file

        else:
            pass

        print ('Genomic index being used for mapping: %s\n'% (genoIndexPrePro))
        #genoIndex = 'ASPARAGUS_UGA1_genome' ## Test
        
        ### Prepare ###########################################
        inFile = '%s.%s' % (lib,ext)
        print ('Processing %s for mapping to genome' % (inFile))
        fastaFile = tagCount2FASTA(inFile,'N') ## Unique reads to FASTA format 

        mapFile   = ('./%s.%s.map' % (lib,ext.rpartition('.')[0]))
        print(genoIndexPrePro,inFile,fastaFile,mapFile)
        
        ## Map to index ##########################################
        print ('Mapping %s processed file to genome' % (lib))
        nproc2 = str(nproc)

        if int(maxTagLen)   > 60:
            mismat = str(2)
        elif int(maxTagLen) <= 60 and maxTagLen > 40:
            mismat = str(1)
        elif int(maxTagLen) <= 40:
            mismat = str(0)
        else:
            print(f"Mismatches not assigned for this read length {maxTagLen}")
            sys.exit()
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])
        
        ## Bowtie 1 - So as to be compatible with current indexes
        retcode = subprocess.call([bowtie, "-f", "-n", mismat, "-p", nproc2, "-t" , genoIndexPrePro, fastaFile, mapFile])
        
        if retcode == 0:    ## Mapping exited with status 0, all is well
            print('\nBowtie mapping for %s complete' % (inFile) )
        else:
            print ("There is some problem with mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()
        
        ### Prepare lists for plotting
        mappedList,mappedAbunList   = mappedStats([lib,ext,nthread,maxTagLen], mode)
        print('Mapped Reads:',mappedList)
        print('Abundance of mapped:',mappedAbunList)
        allTagsList,allAbunList     = tagCountStats([lib,ext,nthread,maxTagLen],mode)
        print('\nAll Reads:',allTagsList)
        print('Abundance of all sizes:',allAbunList)

        ## write stats
        lib_counts_writer(lib, mappedList, mappedAbunList, allTagsList, allAbunList, mode)


        ## Plot
        charts(lib,ext,mappedList,mappedAbunList,allTagsList,allAbunList,mode) ## Mode 1 - preprocess graphs 2: postprocessed files graphs
    
    return None

def lib_counts_writer(lib, mappedList, mappedAbunList, allTagsList, allAbunList, mode):
    '''
    write temp stats files
    '''
    if mode == 1:
        afile = "%s_allBefore.temp"     % (lib)
        mfile = "%s_mappedBefore.temp"  % (lib)
    elif mode == 2:
        afile = "%s_allAfter.temp"      % (lib)
        mfile = "%s_mappedAfter.temp"   % (lib)
    else:
        print(f"The mode value:{mode} for 'mapper' is not recongnized - Debug.")
        sys.exit()

    afh = open(afile, 'w')
    mfh = open(mfile, 'w')

    ## write stats
    afh.write("%s\n" % (','.join(str(i) for i in allTagsList)))
    afh.write("%s" % (','.join(str(i) for i in allAbunList)))
    mfh.write("%s\n" % (','.join(str(i) for i in mappedList)))
    mfh.write("%s" % (','.join(str(i) for i in mappedAbunList)))

    afh.close()
    mfh.close()

    return None


## ASSEMBLY #######################
def indexBuilder(genoFile):
    print ("\n**Deleting old index 'folder' !!!!!!!!!!!***\n")
    print("If its a mistake cancel now by pressing ctrl+D and continue from index step by turning off earlier steps- You have 30 seconds")
    time.sleep(30)
    shutil.rmtree('./index', ignore_errors=True)
    os.mkdir('./index')
    
    nproc2      = str(nproc) 
    genoIndex   = '%s/index/%s' % (os.getcwd(),genoFile.rpartition('/')[-1].rpartition('.')[0]) ## Can be merged with genoIndex from earlier part if we use bowtie2 earlier
    print('#### Creating index of cDNA/genomic sequences:%s**\n' % (genoIndex))
    retcode     = subprocess.call(["hisat2-build","-p",nproc2,"-f", genoFile, genoIndex])
    
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        # print("Reference index prepared sucessfully")
        pass
    else:
        print("There is some problem preparing index of reference '%s'" %  (genoFile))
        print("Is 'Bowtie' installed? And added to environment variable?")
        print("Script will exit now")
        sys.exit()


    return genoIndex

def splicedMapper(rawInputs,genoIndex, gtfFile):
    '''
    Map processed files to genome index
    '''
    
    for aninput in rawInputs:
        print('\nInput:',(aninput))
        lib,ext,nthread = aninput
        
        ### Prepare ###########################################
        if seqType == 0:
            inFile = '%s.%s' % (lib,ext)
            print ('Mapping %s to genome' % (inFile))
            samFile = ('./%s.hisat.sam' % (lib))
            sumFile = ('./%s.hisat2.sum' % (lib))
            print(genoIndex,inFile,samFile)
        
        else:
            inFile1 = '%s.chopped.pair_1.trimmed.fastq' % (lib)
            inFile2 = '%s.chopped.pair_2.trimmed.fastq' % (lib)
            inFiles = [inFile1,inFile2] ## Just for print statement
            print ('#### Mapping %s processed file to genome' % (inFiles))
            samFile = ('./%s.hisat2.sam' % (lib))
            sumFile = ('./%s.hisat2.sum' % (lib))
            print(genoIndex,inFile1,inFile2,samFile)
        
        ## Map to index ##########################################
        nproc2 = str(nproc)
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])

        #### Select Library type mode
        if libType      == 0:
            atype       = 'F'
        elif libType    == 1:
            atype       = 'R'
        elif libType    == 2:
            atype       = 'RF'
        elif libType    == 3:
            atype       = 'FR'
        else:
            print("Value for 'libType' parameter is not recognized:%s" % (libType))
            print("Please check allowed modes - script will exit for now..")
            sys.exit()

        #### Run Spliced Mapper ####
        ############################
        if referenceGTF == 'T':
            ## Use user supplied GTF file
            if seqType == 0:
                ## Single End
                retcode = subprocess.call([hisat2,"-p",nproc2,"--dta","-x",genoIndex,"--known-splicesite-infile",ssFile,"--rna-strandness",atype,"-S",samFile,"-q","-U",inFile,"--summary-file",sumFile])
            else:
                ## Paired End
                retcode = subprocess.call([hisat2,"-p",nproc2,"--dta","-x",genoIndex,"--known-splicesite-infile",ssFile,"--rna-strandness",atype,"-S",samFile,"-q","-1",inFile1,"-2",inFile2,"--summary-file",sumFile])
        
        elif referenceGTF == 'F':
            ## No GTF file supplied by user
            if seqType == 0:
                ## Single End
                retcode = subprocess.call([hisat2,"-p",nproc2,"--dta","-x",genoIndex,"--rna-strandness",atype,"-S",samFile,"-q","-U",inFile,"--summary-file",sumFile])
                ## Paired End
                retcode = subprocess.call([hisat2,"-p",nproc2,"--dta","-x",genoIndex,"--rna-strandness",atype,"-S",samFile,"-q","-1",inFile1,"-2",inFile2,"--summary-file",sumFile])

        else:
            print("Please choose correct option for 'referenceGTF' and provide 'GTF' file if selected - Script will exit for now")
            sys.exit()

        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\nHISAT2 mapping for %s complete' % (inFiles) )
        else:
            print ("There is some problem with HISAT2 mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFiles))
            print ("Script exiting.......")
            sys.exit()

        ## Convert output to sorted BAM
        bamSort = convertToBAM(samFile,nproc2)

    return None

def convertToBAM(samFile,nproc):
    '''
    converts a SAM to sorted BAM
    '''

    #### Prepare BAM Files ######
    #############################
    print("#### Converting and sorting SAM to BAM")
    # bamFile     = "%s.bam" % samFile.rpartition(".")[0]
    bamSort     = "%s.sorted.bam" % samFile.rpartition(".")[0] ### Extention added automatically to BAM

    # retcode     = subprocess.call([samtools, "view", "-Su", samFile, "-o", bamFile])
    # if retcode == 0:
    #     print("%s converted to BAM format" % (samFile))
    # else:
    #     print("conversion of SAM to BAM failed:%s" % (samFile))
    #     sys.exit()


    retcode = subprocess.call([samtools, "sort", "-@", nproc, "-o", bamSort, samFile])
    if retcode == 0:
        print("Sorted %s " % (bamSort))
        pass
    else:
        print("Converion or sorting of SAM file failed:%s" % (samFile))
        sys.exit()


    return bamSort

def stringTie(aninput):
    '''
    generates output for edgeR and ballgown analyses
    '''

    print("\n#### Fn: StringTie Analysis")
    print('Inputs:',(aninput))
    lib,ext,nthread,amode   = aninput
    nthread2                = str(nthread)
    
    if amode == 1:
        ## used to generate gene/transcript 
        ## assembly
        inFile  = '%s.%s' % (lib,ext) 
        outFile = ('%s.stringtie.gtf' % (lib))
        print("\n Generating library-specific assemblies")
        print ('Input:%s | Output:%s' % (inFile,outFile))
        retcode = subprocess.call([stringtie, "-p", nthread2, "-G", gtfFile, "-o", outFile, "-l", lib, inFile])

    elif amode == 2:
        ## used to map reads from different 
        ## libraries to merged assembly from
        ## first round of stringtie on all
        ## different libraries
        inFile      = '%s.%s' % (lib,ext) 
        outFile     = '%s.stringtie.quant.gtf' % (lib)
        mergedgtf   = 'stringtie_merged.gtf' 
        print("\n Generating summaries for the merged assembly (in accordance with `referenceGTF` value)")
        print ('Input:%s | Output:%s' % (inFile,outFile))
        retcode = subprocess.call([stringtie, "-p", nthread2, "-G", mergedgtf, '-e', "-o", outFile, "-l", lib, inFile])

    else:
        print("Wrong mode '%s'provided for stringTieStep" % (amode))
        print("Script exiting")
        sys.exit()

    #### Sanity Check ######
    if retcode == 0:
        print("StringTie finished for:%s " % (inFile))
        pass
    else:
        print("StringTie analysis failed:%s" % (inFile))
        sys.exit()


    return None

def stringMerge(rawInputs,gtfFile,genoFile):
    
    '''
    merge all output for stringtie from pheno file; first generate a merge list and then 
    supply that for merging
    https://github.com/griffithlab/rnaseq_tutorial/wiki/Transcript-Assembly-Merge
    '''
    
    ## Make a list of lib-wise assemblies to merge
    print ("\n#### Fn: StringMerge")
    assemblyFile    = 'mergelist.txt'
    assemblyOut     = open(assemblyFile,'w')
    for aninput in rawInputs:  
        lib,ext,nthreads    = aninput
        assemblyPath        = './%s.%s\n' % (lib,ext)
        print("Adding lib to list for merging:%s" % (assemblyPath.strip('\n')))
        assemblyOut.write(assemblyPath)

    assemblyOut.close()
    
    #### Use mergelist to merge assemblies
    print ("\n#### Merging lib-specific assemblies from StringTie")
    mergedAssembly  = "stringtie_merged.gtf"
    nproc2          = str(nproc)
    
    
    if referenceGTF == 'T':
        print("CMD:", stringtie, "--merge", "-p", nproc2, "-G", gtfFile, "-o", mergedAssembly, assemblyFile)
        retcode = subprocess.call([stringtie, "--merge", "-p", nproc2, "-G", gtfFile, "-o", mergedAssembly, assemblyFile])
    else:
        print("CMD:", stringtie, "--merge", "-p", nproc2, "-o", mergedAssembly, assemblyFile)
        retcode = subprocess.call([stringtie, "--merge", "-p", nproc2, "-o", mergedAssembly, assemblyFile])

    if retcode == 0:
        print("Lib-speciifc assemblies merged successfully\n\n")
    else:
        print("Problem merging lib-specific assemblies - System wil exit now\n")
        sys.exit()
    
    return None
 
def stringQuant(rawInputs):
    '''
    Run stringtie with merged assembly and then compute gene/transcript estimates
    
    prepDE bundled with StringTie requires python3;
    download python3 version from here: https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
    direct link: https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3
    '''
    print('\n#### Fn: StringQuants\n')

    #### prepare samplelist
    print ("\n#### Generating list of lib-specific assemblies to generate counts")
    quantsFile    = 'sample_lst.txt'
    fh_out        = open(quantsFile,'w')
    acount        = 0
    for aninput in rawInputs:
        acount          += 1
        lib,ext         = aninput
        assemblyPath    = '%s\t%s.%s' % (lib,lib,ext)
        print("Adding lib to list for merging:%s" % (assemblyPath))
        fh_out.write("%s\n" % (assemblyPath))
    fh_out.close()

    #### generate gene/transcript counts
    print ("\nNOTE: Use Python3 version of prepDE.py")
    print ("NOTE: it's bundled with Rocket.v4, also can be downloaded from StringTie webpage")
    print ("NOTE: see notes for `stringQuant` finction for links")  
    retcode     = subprocess.call(["python",prepDE,"-i", quantsFile])
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        print('\nCounts file generated from %s libraries' % (acount))
        print("See 'sample_lst.txt' for list of libraries\n")
        print("Output: 'gene_count' and 'transcript_coun' CSV files\n")
    else:
        print ("Problem generating counts table")
        sys.exit()
        
    return None


## HELPERS ########################
def sampleInfoRead(sampleInfo):
    ''' This module reads a sample info file to make a list of 
    libraries/files, replicates and groups'''
    
    print("\nFunction - sampleInfoRead")

    fh_in = open(sampleInfo,'r')
    fh_in.readline() ## Remove header
    sampleRead = fh_in.readlines()

    libs = []           ## List to hold libraries (server mode) or files (local mode)
    reps = []           ## List to hold replicates
    for i in sampleRead:
        ent     = i.strip("\n").split("\t")
        arep    = ent[2]
        agroup  = ent[3]

        ## Filename or library id as an identifier
        if Local == 1:
            anid = ent[1]
        elif Local == 0:
            anid = ent[0]
        else:
            print("Please correct input mode - Server or Local")
            print("System will exit now")
            sys.exit()

        libs.append(anid)

        ## Group on replicates or some other group
        if groupBy == 'R':
            reps.append(arep)
        elif groupBy == 'G':
            reps.append(agroup)
        else:
            print("Please choose correct sample grouping method")
            print("System will exit now")
            sys.exit()

    print("This is 'libs':",libs)
    print("These are 'reps':",reps)

    print("Total files from sampleInfo:%s | Total number of replicates:%s" % (str(len(libs)),str(len(reps))))
    print("Exiting function - sampleInfoRead\n")
    
    return libs,reps

def PP(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(module, alist)

def PPBalance(module,alist):
    #print('***********Parallel instance of %s is being executed*********' % (module))
    start   = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = int(nproc/int(nthread)) ## 1 added so as to avoid 0 processor being allocated in serial mode

    ## Sanity check
    if nprocPP == 0:
        nprocPP = 1
    else:
        pass

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

def mappedStats(aninput,mode):
    '''Parse map file and collect statics for graph generation'''
    print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('\nCollecting statistics of matched reads for Lib:%s' % (lib))
    
    inFile = '%s.%s.map' % (lib,ext.rpartition('.')[0])
    fh_in = open(inFile,'r')
    mapFile = fh_in.read().split('\n')
    
    if mode == 1:
        mappedList      = [0]*(maxReadLen) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        mappedAbunList  = [0]*(maxReadLen) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
    elif mode == 2:
        mappedList      = [0]*(maxTagLen+1) ## List lenght equal to max size of fragment (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
        mappedAbunList  = [0]*(maxTagLen+1) ## List length equal to max size of fragment +1 (to compensate the python indexing) for ex. max chopped length = 34 than list should have 35 slots - Can be put in settings
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

def tagCountStats(aninput,mode):

    '''Get stats for all the reads from tagCount file'''
    
    #print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('\nCollecting statistics of total reads for Lib:%s' % (lib))
    
    inFile = '%s.%s' % (lib,ext)
    print(inFile)
    fh_in           = open(inFile,'r')
    tagCountFile    = fh_in.read().split('\n')
    
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

def charts(lib,ext,mappedList,mappedAbunList,allTagsList,allAbunList,mode):
    
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
        plotFile = ('%s.%s_distinct_before_chop.png' % (lib,ext.rsplit('.',2)[0]) )##Plot results file
        
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
        plt.xticks(np.arange(N),    np.arange(int(minLen),int(maxLen)+1),   rotation = 45,  fontproperties=font_manager.FontProperties(size=6))
        plt.yticks(np.arange(0,maxAbun+ybreak,ybreak),  fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))
        
        plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file
        
        
        
        ### CHART 2:  Mapped abundance vs total abundance ######
        plotFile2 = ('%s.%s_abund_before_chop.png' % (lib,ext.rsplit('.',2)[0]))##Plot results file
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
        plt.xticks(np.arange(N), np.arange(int(minLen),int(maxLen)+1),  rotation = 45,  fontproperties=font_manager.FontProperties(size=6) )
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
        plotFile = ('%s.%s_distinct_after_chop.png' % (lib,ext.rsplit('.',2)[0]))##Plot results file
        
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
        plt.xticks(np.arange(N),    np.arange(int(minLen),int(maxLen)+1),   rotation = 45,  fontproperties=font_manager.FontProperties(size=6))
        plt.yticks(np.arange(0,maxAbun+ybreak,ybreak),  fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))
        
        plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file
        
        
        ### Chart 2:  All genomic reads
        plotFile2 = ('%s.%s_abund_after_chop.png' % (lib,ext.rsplit('.',2)[0]))##Plot results file
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
        plt.xticks(np.arange(N), np.arange(int(minLen),int(maxLen)+1),  rotation = 45,  fontproperties=font_manager.FontProperties(size=6) )
        plt.yticks(np.arange(0,maxAbun2+ybreak2,ybreak2),   fontproperties=font_manager.FontProperties(size=6))
        ax = plt.gca()
        ax.yaxis.grid(True)
        plt.legend((p1[0], p2[0]), ('Mapped reads','Total Reads'), loc=1, prop=font_manager.FontProperties(size=6))

        plt.savefig(plotFile2, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
        plt.clf() ## Clear figure so that next plot is different file

    else:
        print('\nThe mode selected for generating graph is not correct - Debug for reason')
        
    
    return None

def writeStats(aninput):
    ''' 
    Write a lib-specific summary file
    '''
    lib,minTagLen,maxTagLen = aninput
    print (aninput)

    ### Read the library temp files with counts and abundances before processing
    print("Reading stats from: %s_allBefore.temp" % lib)
    fh_before           = open("%s_allBefore.temp" % lib,'r')
    aread               = fh_before.read().split('\n')
    allTags,allAbun     = aread
    allTagsList,allAbunList = list(map(int,allTags.split(','))),list(map(int,allAbun.split(',')))
    # print("allTagsList:",allTagsList,"\nallAbunList:",allAbunList)
    fh_before.close()

    ### Read the library temp files with counts and abundances after processing
    print("Reading stats from: %s_allAfter.temp" % lib)
    fh_after            = open("%s_allAfter.temp" % (lib),'r')
    aread2              = fh_after.read().split('\n')
    allTags2,allAbun2   = aread2
    allTagsList2,allAbunList2 = list(map(int,allTags2.split(','))),list(map(int,allAbun2.split(',')))
    fh_after.close()

    fh_map_before = open("%s_mappedBefore.temp" % (lib), 'r')
    aread3 = fh_map_before.read().split('\n')
    mapped,mappedAbun = aread3
    mappedList,mappedAbunList = list(map(int,mapped.split(','))),list(map(int,mappedAbun.split(',')))
    fh_map_before.close()
    
    fh_map_after = open("%s_mappedAfter.temp" % (lib),'r')
    aread4 = fh_map_after.read().split('\n')
    mapped2,mappedAbun2 = aread4
    mappedList2,mappedAbunList2 = list(map(int,mapped2.split(','))),list(map(int,mappedAbun2.split(',')))
    fh_map_after.close()

    ### Prepare to write
    summFile = "%s_chopinfo.txt" % lib
    fh_out = open(summFile,'w')
    fh_out.write("Date - %s | Genome - %s\n" % (time.strftime("%d/%m/%Y"),genomeDB))
    fh_out.write("Lib-%s\tBeforeProcessing\tAfterProcessing\tMappedBeforeProcessing\tmappedAfterProcessing\n" % (lib))

    # print("allTagsList length:",len(allTagsList),"\nallAbunList length:",len(allAbunList))
    indexList   = [i for i,x in enumerate(allTagsList) if x != 0]
    # print ('indexList:',indexList)
    minLen      = min(indexList)
    maxLen      = max(indexList)

    indexList2  = [i for i,x in enumerate(allAbunList) if x != 0]
    # print ('indexList2:',indexList)
    minLenAbun  = min(indexList2)
    maxLenAbun  = max(indexList2)
    # print("Allowed min len:%s | Allowed max len:%s" % (minTagLen,maxTagLen))
    # print("allTags min len:%s | allTags max len = %s | allAbun max len:%s | allAbun max len:%s\n" % (minLen,maxLen,minLenAbun,maxLenAbun))

    countsAllBefore     = list(allTagsList[minLen:maxLen+1])
    countsAllAfter      = list(allTagsList2[minLen:maxLen+1])
    abunAllBefore       = list(allAbunList[minLen:maxLen+1])
    abunAllAfter        = list(allAbunList2[minLen:maxLen+1])
    print(countsAllBefore,countsAllAfter,abunAllBefore,abunAllAfter)

    countsMappedBefore  = list(mappedList[minLen:maxLen+1])
    countsMappedAfter   = list(mappedList2[minLen:maxLen+1])
    abunMappedBefore    = list(mappedAbunList[minLen:maxLen+1])
    abunMappedAfter     = list(mappedAbunList2[minLen:maxLen+1])

    ## Write file for every tag size
    indBefore = 0 ## Index to keep track of psoition in all tags (before processing list)
    indAfter = 0 ## Index to keep track of position in processed lists
    for i in indexList:
        # print("Size of tag:%s" % (i))
        ## use size info from wishlist, to nter zero if size has been filtered out in processing
        if i >= minTagLen and i <= maxTagLen:
            fh_out.write("%snt-Count\t%s\t%s\t%s\t%s\n" % (i,countsAllBefore[indBefore],countsAllAfter[indAfter],countsMappedBefore[indBefore],countsMappedAfter[indAfter]))
            fh_out.write("%snt-Abundance\t%s\t%s\t%s\t%s\n" % (i,abunAllBefore[indBefore],abunAllAfter[indAfter],abunMappedBefore[indBefore],abunMappedAfter[indAfter]))
            indAfter += 1
            indBefore+=1

        else:
            fh_out.write("%snt-Count\t%s\t0\t%s\t0\n" % (i,countsAllBefore[indBefore],countsMappedBefore[indBefore]))
            fh_out.write("%snt-Abundance\t%s\t0\t%s\t0\n" % (i,abunAllBefore[indBefore],abunMappedBefore[indBefore]))
            indBefore+=1

    fh_out.close()

    return None

def parseCoordFile(coordFile):

    '''Parse Rocket merged GTF file and prepare
    a list of coords'''

    print("\nModule:parseCoordFile")
    print("Parsing GTF file and preparing coords dictionary\n")

    coordsList  = [] ## List to store coords for sequence extraction
    coordsDict  = {} ## Dictionary with unique gene_id ans key and other info as values
    transList   = [] ##List of transcript coords - Generated only in mode 2

    geneSet     = set()
    transSet    = set()

    fh_in = open(coordFile,'r')
    
    if (mode == 1) or (mode == 2):
        fileRead = fh_in.readlines()
        faultCount = 0

        for ent in fileRead:
            # print(ent)
            chr_id,trash1,trash2,start,stop,dot1,strand,dot2,info = ent.strip('\n').split('\t')
            # print(chr_id,start,stop,strand,info)
            info_splt       = info.split(';')[:4]
            gene_id         = info_splt[0].split('"')[1].replace('"','')
            transcript_id   = info_splt[1].split('"')[1].replace('"','')
            exon            = info_splt[2].split('"')[1].replace('"','')
            gene_name       = info_splt[3].split('"')[1].replace('"','')
            # print(chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name)
            geneSet.add(gene_id)
            transSet.add(transcript_id)
            
            ## Change strand format
            if strand == '+':
                strand = 'w'
            elif strand == '-':
                strand = 'c'
            else:
                # print('Wrong strand encountered while splitting: %s' % strand)
                # print(ent)
                faultCount += 1

            ## Record
            coordsList.append((chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name))
            coordsDict[gene_id] = ((chr_id,strand,gene_name))
    ## Extract entries and add gene strand from coordsList
    
    if mode == 2:
        print("Preparing coords for transcript extraction")
        fh_in = open(transFile,'r')
        fh_in.readline() ## Remove header
        readFile = fh_in.readlines()

        for line in readFile:
            ent         = line.strip('\n').split('\t')
            # print(ent)
            gene_id     = ent[0]
            gene_name   = ent[4]
            coords      = ent[6]
            # print(coords)
            chr_id      = coords.split(':')[0]
            start2,stop2 = coords.split(':')[1].split('-')

            ## Get strand from the coordsDict
            chr_id,strand,gene_name = coordsDict[gene_id]
            transList.append((chr_id,start2,stop2,strand,gene_id,gene_name))

    ##Sort the list on chromosome,strand and stop
    sorted(coordsList, key=operator.itemgetter(0,3,1))
    sorted(transList, key=operator.itemgetter(0,3,1))

    # print(coordsList[:3])
    print('Total entries recorded:%s | Entries with problem:%s\n'% (len(coordsList),faultCount))
    print('Total genes recorded:%s | transcripts recorded:%s' % (len(geneSet),len(transSet)))
    time.sleep(2)
    fh_in.close()

    return coordsList,transList

def makeGenomeDict(genoFile):

    ''' Prepare a dictionary of genome with chr/scaffold as key
    and seq as value'''

    print("\nModule:makeGenomeDict")
    print("Preparing genome dictionary\n")

    fh_in       = open(genoFile,'r')
    readFile    = fh_in.read().split('>')
    genomeDict  = {} 

    for i in readFile[1:]:
        # print(i.split('\n'))
        chromo  = i.split('\n')
        head    = chromo[0].split(' ')[0] ## Shorten the header to be used as key
        chrSeq  = '' ## Empty string for sequence
        for i in chromo[1:]:
            # print("seq:%s" % (i))
            chrSeq+=i
        # print(">%s\n%s" % (head,chrSeq))
        
        ## Add to dict
        genomeDict[head] = chrSeq

    print("Entries in the genomeDict:%s\n" % (len(genomeDict)))
    return genomeDict

def extractCDS(coordsList,genomeDict):

    ''' Module extracts CDS sequences from the merged gtf fileby stiching togther the exons for every transcript
    This includes isoform for same genes - I am lazy, could have used mySQLite'''

    print("\nModule:extractCDS")
    print("Constructing CDS from Exons")

    featureSeq  = '' ## Intialize empty sequence and then empty after every feature (below)
    seqList     = [] ## List to store results
    chrLoaded   = '' ## Will keep track of running chromosome, since coordsList is sorted this will reduce runtime
    exonCount   = 0 ## Keep counts of exons for a gene
    
    # for ents in zip(coordsList[0:],coordsList[1:]):
    for i in range(0,int(len(coordsList))-1):
        ent1    = coordsList[i]
        ent2    = coordsList[i+1]
        chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name = ent1
        print("Ent:",chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name)
        ## To accomodate missing one nucleotide at start
        if start == 0:
            pass ## Bug-1 fixed
        else:
            start = int(start)-1

        ## Load Chromosome 
        if chrLoaded:
            ## This is not the first entry and chromosome has been loaded
            if chr_id == chrLoaded:
                ##No need to fetch chromosme
                pass
            else:
                ## The chromosome is different from last entry - Load and update chrLoaded
                chromo = genomeDict[ent1[0]] ## Use chr_id to fetch sequence
                chrLoaded = chr_id ##update chrLoaded
        else:
            ## This is the first entry - Load chromosme and update chrLoaded
            chromo = genomeDict[ent1[0]] ## Use chr_id to fetch sequence
            chrLoaded = chr_id ## update chrLoaded

        ## Fetch Sequences
        if ent2:
            ## Not EOF
            if transcript_id == ent2[5]:
                ## Exons belong to same isoform
                # chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
                exonSeq = chromo[int(start):int(stop)]
                # print(exonSeq)
                featureSeq+=exonSeq
                exonCount +=1
            
            else:
                if exonCount >= 1:
                ## This is last exon of multi exon gene
                    # chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
                    exonSeq = chromo[int(start):int(stop)]
                    print("This is final exon of multi-exon transcript:%s" % (exonSeq))
                    featureSeq+=exonSeq
                    # print("This is final transcript:%s" % (featureSeq))
                    if strand == 'w':
                        seqList.append((gene_id,transcript_id,featureSeq))
                    elif strand == 'c':
                        seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
                    else:
                        ## Strand is uknown and must be '.' in file - So record both
                        print("GeneID:%s - Unknown strand" % (gene_id))
                        seqList.append((gene_id,transcript_id,featureSeq))
                        seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
                    ## Current transcript is complete - reintialize transcript-specific variables
                    featureSeq = '' ## Empty for next gene as current transcript is complete
                    exonCount = 0 ## Set end of gene

                else:
                    ## This is the first Exon of single exon transcript
                    featureSeq = '' ## Just to make sure
                    # chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
                    exonSeq = chromo[int(start):int(stop)]
                    print("This is single exon transcript:%s\%s" % (ent1,exonSeq))
                    featureSeq += exonSeq
                    if strand == 'w':
                        seqList.append((gene_id,transcript_id,featureSeq))
                    elif strand == 'c':
                        seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
                    else:
                        ## Strand is uknown and must be '.' in file - So record both
                        print("GeneID:%s - Unknown strand" % (gene_id))
                        seqList.append((gene_id,transcript_id,featureSeq))
                        seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))

                    ## Current transcript is complete - reintialize transcript-specific variables
                    featureSeq = '' ## Empty for next gene as current transcript is complete
                    exonCount = 0 ## Set end of gene
        else:
            ## This is EOF - Final entry - Wrap up the transcript
            print ("This is the final entry")
            if exonCount >= 1:
            ## This is last exon of multi exon gene
                # chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
                exonSeq = chromo[int(start):int(stop)]
                print("This is final exon of multi-exon transcript:%s" % (exonSeq))
                featureSeq+=exonSeq
                # print("This is final transcript:%s" % (featureSeq))
                if strand == 'w':
                    seqList.append((gene_id,transcript_id,featureSeq))
                elif strand == 'c':
                    seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
                else:
                    ## Strand is uknown and must be '.' in file - So record both
                    print("GeneID:%s - Unknown strand" % (gene_id))
                    seqList.append((gene_id,transcript_id,featureSeq))
                    seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))


                ## Current transcript is complete - reintialize transcript-specific variables
                featureSeq = '' ## Empty for next gene as current transcript is complete
                exonCount = 0 ## Set end of gene
                break

            else:
                ## This is the first Exon of single exon transcript
                featureSeq = '' ## Just to make sure
                # chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
                exonSeq = chromo[int(start):int(stop)]
                print("This is single exon transcript:%s\%s" % (ent1,exonSeq))
                featureSeq += exonSeq
                if strand == 'w':
                    seqList.append((gene_id,transcript_id,featureSeq))
                elif strand == 'c':
                    seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
                else:
                    ## Strand is uknown and must be '.' in file - So record both
                    print("GeneID:%s - Unknown strand" % (gene_id))
                    seqList.append((gene_id,transcript_id,featureSeq))
                    seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))

                ## Current transcript is complete - reintialize transcript-specific variables
                featureSeq = '' ## Empty for next gene as current transcript is complete
                exonCount = 0 ## Set end of gene
                break
    
    print("Total CDS extracted: %s" % (len(seqList)))
    return seqList

def extractTrans(genomeDict,transList):
    ''' Module takes gene attributes file from Rocket
    to extract transcripts '''

    print("\nModule:extractTrans")
    print("Constructing CDS from Exons")
    
    chrLoaded = '' ## Will keep track of running chromosome, since coordsList is sorted this will reduce runtime
    transSeqList = [] ## List to record results

    for i in transList:
        chr_id,start,stop,strand,gene_id,gene_name = i
        print("This is the ent:",i)
        ## To accomodate missing one nucleotide at start
        if start == '0':
            # print("Start is Zero")
            start = int(start) ## Bug-1 fixed
        else:
            start = int(start)-1

        ## Load Chromosome 
        if chrLoaded:
            ## This is not the first entry and chromosome has been loaded
            if chr_id == chrLoaded:
                ##No need to fetch chromosme
                pass
            else:
                ## The chromosome is different from last entry - Load and update chrLoaded
                chromo = genomeDict[chr_id] ## Use chr_id to fetch sequence
                chrLoaded = chr_id ##update chrLoaded
        else:
            ## This is the first entry - Load chromosme and update chrLoaded
            chromo = genomeDict[chr_id] ## Use chr_id to fetch sequence
            chrLoaded = chr_id ## update chrLoaded

        # print("This is the chromosome:",chr_id,chromo)

        ## Extract Transcript - Both strand in case strand is unknown
        if strand == 'w':
            featureSeq = chromo[int(start):int(stop)]
            # print("FeatureSeq:",featureSeq)
            transSeqList.append((gene_id,gene_name,featureSeq))
            
        elif strand == 'c':
            featureSeq_rc = chromo[int(start):int(stop)][::-1].translate(str.maketrans("TAGC","ATGC"))
            # print("FeatureSeq:",featureSeq_rc)
            transSeqList.append((gene_id,gene_name,featureSeq_rc))
            
        elif strand == '.':
            ## Strand is uknown and must be '.' in file - So record both
            print("GeneID:%s - Unknown strand" % (gene_id))
            featureSeq = chromo[int(start):int(stop)]
            featureSeq_rc = chromo[int(start):int(stop)][::-1].translate(str.maketrans("TAGC","ATGC"))
            # print("FeatureSeq:",featureSeq)
            # print("FeatureSeq:",featureSeq_rc)
            
            transSeqList.append((gene_id,gene_name,featureSeq))
            transSeqList.append((gene_id+'_rc',gene_name+'_rc',featureSeq_rc))

    return transSeqList

def extractProm(coords,genoFile):
    pass

def writer(resList):
    '''Write the CDS, transcripts or promoter identified above
    resList has three entries - ID1,ID2 and Seq'''

    print("\nModule:writer")
    print("Writing CDS, transcripts or promoter")

    if mode == 1:
        outFile = '%s.fa' % (coordFile.rpartition('.')[0])
        fh_out = open(outFile,'w')

        for i in resList:
            # print(i)
            gene_id,transcript_id,featureSeq = i
            fh_out.write('>%s_%s\n%s\n' % (gene_id,transcript_id,featureSeq))
        pass

    elif mode == 2:
        outFile = '%s.fa' % (transFile)
        fh_out = open(outFile,'w')

        for i in resList:
            # print(i)
            gene_id,gene_name,featureSeq = i
            fh_out.write('>%s_%s\n%s\n' % (gene_id,gene_name,featureSeq))
        pass

def coreReserve(ncores):
    '''
    Decides the core pool for machine - written to make PHASIS comaptible with machines that 
    have less than 10 cores - Will be improved in future
    '''
    print(f"\n#### Fn: coreReserve")
    totalcores = int(multiprocessing.cpu_count())

    if ncores == 0:
        ## Automatic assignment of cores selected       
        if totalcores   == 4: ## For quad core system
            ncores = 3
        elif totalcores == 6: ## For hexa core system
            ncores = 5
        elif totalcores > 6 and totalcores <= 10: ## For octa core system and those with less than 10 cores
            ncores = 7
        else:
            ncores = int(totalcores*0.85)
    else:
        ## Reserve user specifed cores
        ncores = int(ncores)

    print(f"Reserved {ncores}/{totalcores} cores")
    return ncores

def splicesitemap(gtfFile):
    '''
    converts GTF file to input required by the hisat2
    '''

    subprocess.Popen("script2.py 1", shell=True)
    inFile = "xxx"
    return None


## MAIN ###########################
def main(sampleInfo):
    
    start = time.time() 
    runLog = open('%s_run.log' % (datetime.datetime.now().strftime("%m_%d_%H_%M")), 'w')
    if userMinTagLen:
        pass
    
    #### 0. Initialize Input Register ###################################################
    register =[] ## Holds values for all different steps - master feeder
    libs,rep = sampleInfoRead(sampleInfo)
    print('\n\n**Total %s libraries provided for pre-processing**\n' % (len(libs)))
    
    if Local == 1: ## Local analysis
            for i,k in zip(libs,rep):
                ext = 'fastq'
                minTagLen = userMinTagLen
                maxTagLen = userMaxTagLen
                filePath = './%s.%s' % (i.strip(),ext) ## Mayumi add a space before lib_id/lib_code, strip removes empty space 
                register.append((str(i),ext,str(nthread),str(filePath),adapterFileSE,adapterFilePE,minTagLen,maxTagLen,str(k))) ## File ID, ext, nthread, raw file path,None (added to maintain same structure as remote analysis register),
                                                                                                            ## adapter file, min tag len, max tag len, replicate

    elif Local == 0: ## Remote analysis
        global con
        con = ConnectToDB(dataServer,0)
        for i,k in zip(libs,rep):
            cur = con.cursor()
            cur.execute("SELECT raw_path,raw_file_format,adapter_5p,adapter_3p,minimum_tag_length,maximum_tag_length FROM master.library_info WHERE lib_id like '%s'" % (i)) ##bowtie_index_path
            fileInfo = cur.fetchall()
            print ('Library:', (fileInfo[0]))
            
            filePath = fileInfo[0][0].replace('$ALLDATA', '/alldata') 
            ext = fileInfo[0][1]
            adp_5p = fileInfo[0][2]
            adp_3p = fileInfo[0][3]
            if hardMinTagLen == 'Y':
                minTagLen = userMinTagLen
                maxTagLen = userMaxTagLen
            else:
                minTagLen = fileInfo[0][4]
                maxTagLen = fileInfo[0][5]
            register.append((str(i),ext,str(nthread),str(filePath),adp_5p,adp_3p,minTagLen,maxTagLen,str(k))) ## Lib, ext, nthread, raw file path, adapter 5p, adapter 3p, min tag len, max tag len, replicate

    else:
        print("Please choose correct value for Local variable: Local = [0/1]")
    
    #####################################################################################
    #### 1. QC Raw files ################################################################

    ## SE and PE on Local and Server Ready
    if QCheckStep == 1:
        print('\n**Quality check of the raw files will be performed now**\n')

        if Local == 1: # Provide local file name
            if seqType == 0: ## Single-end
                rawInputs = [(i[0],i[1],i[2],'%s.fastq' % (i[0])) for i in register] ## ## Library/Filename, extension, nthread, local filepath
            elif seqType == 1: ## Paired end
                inputsR = [(i[0],i[1],i[2],'%s_1.fastq' % (i[0])) for i in register] ## ## Library/Filename, extension, nthread, local filepath
                inputsL = [(i[0],i[1],i[2],'%s_2.fastq' % (i[0])) for i in register] ## ## Library/Filename, extension, nthread, local filepath
                rawInputs = inputsL+inputsR
            else:
                print("Please input seqType - SE or PE")

        elif Local == 0: ## Provide input file path
            if seqType == 0: ## Single End
                rawInputs = [(i[0],i[1],i[2],i[3]) for i in register] ## ## Library/Filename, nthread, extension, and alldata filepath
            elif seqType == 1: ## PairedEnd
                inputsR = [(i[0],i[1],i[2],'%s_1.fastq' % (i[3].rpartition("_")[0])) for i in register] ## ## Library/Filename, nthread, extension, alldata filepath
                inputsL = [(i[0],i[1],i[2],'%s_2.fastq' % (i[3].rpartition("_")[0])) for i in register] ## ## Library/Filename, nthread, extension, alldata filepath
                rawInputs = inputsL+inputsR
            else:
                print("Please input seqType - SE or PE")

        else:
            print("Please input correct run mode - Local or Server")
            rawInputs = None
            sys.exit()

        PPBalance(QCheck,rawInputs)
    
    else:
        print ('\n**Quality check of the raw files will be skipped as selected**\n')
    
    ####################################################################################
    #### 2. TRIM RAW files #############################################################
    
    ## SE and PE on Local and Server Ready
    if trimLibsStep == 1:

        rawInputs = [(i[0],'fastq',i[2],i[3],i[4],i[5],i[6]) for i in register] ## Library/Filename, extension, nthread, raw file path, adapter 5p/None(local), adapter 3p/adapter file (local), min tag len
        print('\n**Trimming of the libraries will be performed now**\n')
        # for i in rawInputs:
        #    trimLibs(i)
        PPBalance(trimLibs,rawInputs)
    else:
        print('\n**Trimming of the libraries will be skipped as selected**\n')
        pass

    #####################################################################################
    #### 3. Prepare graphs of trimmed files #############################################
    
    if preProGraphsStep == 1:

        if seqType == 0: ## SingleEnd
            rawInputs = [(i[0],'trimmed.fastq',i[2]) for i in register] ## ## Library/Filename, extension
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.trimmed.fastq',i[2]) for i in register] ## ## Library/Filename, extension
            inputsL = [(i[0],'pair_2.trimmed.fastq',i[2]) for i in register] ## ## Library/Filename, extension
            rawInputs = inputsL+inputsR
        print('\n**Converting trimmed files to tag count format for quality graphs**\n')
        PP(fastQ2Count,rawInputs)
        

        if seqType  == 0: ## SingleEnd
            rawInputs = [(i[0],'trimmed.processed.txt',i[2],i[7], genoIndexPrePro) for i in register] ## ## Library/Filename, extension, max tag len
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.trimmed.processed.txt',i[2],i[7],genoIndexPrePro) for i in register] ## ## Library/Filename, extension, max tag len
            inputsL = [(i[0],'pair_2.trimmed.processed.txt',i[2],i[7],genoIndexPrePro) for i in register] ## ## Library/Filename, extension, max tag len
            rawInputs = inputsL+inputsR
        print('\n**Mapping to generate pre-chopping quality graphs graphs**\n')
            # maps = mapper(rawInputs,1)
        
        if Local == 0:
            # con = ConnectToDB(dataServer,0)
            maps = mapper(rawInputs,1)
        else:
            # con = ConnectToDB(dataServer,0)
            maps = mapper(rawInputs,1)

        
        ### Delete tag count, fasta files and map files to make way for real processed files
        print ("** Cleaning temp files **")
        garbage = [file for file in os.listdir('./') if file.endswith (('.map','trimmed.processed.txt','.processed.fa'))]
        for file in garbage:
            print("Deleting %s" % (file))
            os.remove(file)
    else:
        pass

    ######################################################################################    
    #### 3. Chop trimmed files ###########################################################
    
    ## SE and PE on Local and Server Ready
    if seqType == 0: ## Single end
        rawInputs = [(i[0],'trimmed.fastq',i[2],i[7]) for i in register] ## ## Library/Filename, 'trimmed.fastq', max tag len
    else: ## Paired end
        inputsR = [(i[0],'pair_1.trimmed.fastq',i[2],i[7]) for i in register] ## ## Library/Filename, 'trimmed.fastq', max tag len
        inputsL = [(i[0],'pair_2.trimmed.fastq',i[2],i[7]) for i in register] ## ## Library/Filename, 'trimmed.fastq', max tag len
        rawInputs = inputsL+inputsR

    if chopLibsStep == 1:
        print('\n**Chopping of the libraries will be performed now**\n')
        PPBalance(chopLibs,rawInputs)
    else:
        print('\n**Chopping of the libraries will be skipped as selected**\n')
        pass

    ####################################################################################    
    #### 3B. QC chopped-trimmed files ##################################################
    
    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],i[1],i[2],'%s.chopped.trimmed.fastq' % (i[0],)) for i in register] ## ## Library/Filename,nthread, input file
    else: ## PairedEnd
        inputsR = [(i[0],i[1],i[2],'%s.chopped.pair_1.trimmed.fastq' % (i[0],)) for i in register] ## ## Library/Filename, nthread, input file
        inputsL = [(i[0],i[1],i[2],'%s.chopped.pair_2.trimmed.fastq' % (i[0],)) for i in register] ## ## Library/Filename, nthread,  input file
        rawInputs = inputsL+inputsR
    
    if QCheckStep == 1:
        print('\n**Quality check of the trimmed-chopped files will be performed now**\n')
        PPBalance(QCheck,rawInputs)
    else:
        print ('\n**Quality check of the trimmed-chopped files will be skipped as selected**\n')

    
    #### 4. Convert  chopped files to tagcount format ###################################

    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],'chopped.trimmed.fastq',i[2]) for i in register] ## ## Library/Filename,extension, nthread
    else: ## PairedEnd
        inputsR = [(i[0],'chopped.pair_1.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension, nthread
        inputsL = [(i[0],'chopped.pair_2.trimmed.fastq', i[2]) for i in register] ## ## Library/Filename, extension, nthread
        rawInputs = inputsL+inputsR

    if fastQ2CountStep == 1:
        print('\n**Converting processed files to tag count format**\n')
        PP(fastQ2Count,rawInputs)
    else:
        pass

    ####################################################################################
    #### 5. Map tagcount to genome index ################################################

    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],'chopped.trimmed.processed.txt',i[2],i[7],genoIndexPrePro) for i in register] ## ## Library/Filename,extension, nthread
    else: ## PairedEnd
        inputsR = [(i[0],'chopped.pair_1.trimmed.processed.txt', i[2],i[7],genoIndexPrePro) for i in register] ## ## Library/Filename, extension, nthread, max tag len
        inputsL = [(i[0],'chopped.pair_2.trimmed.processed.txt', i[2],i[7],genoIndexPrePro) for i in register] ## ## Library/Filename, extension, nthread, max tag len
        rawInputs = inputsL+inputsR
    
    if mapperStep == 1:
        print('\n**Mapping chopped files to the genome index**\n')
        # con = ConnectToDB(dataServer,0)
        mapper(rawInputs,2)
    else:
        print("\n**Quality-Mapping Step skipped as selected - What are your intentions?\n")
        pass


    #### 6. Write summary file #########################################################
    print("Writing summary files")

    if summaryFileStep == 1:
        rawInputs = [(i[0],i[6],i[7]) for i in register] ## lib, minTagLen,maxTagLen
        ## Test - Serial
        for i in rawInputs:
            writeStats(i)

    ####################################################################################
    #### 6. Clean up ###################################################################
    if cleanupStep == 1:
        print ("**Cleaning temp files last time**")
        garbage = [afile for afile in os.listdir('./') if afile.endswith (('trim.log','processed.fa','.zip'))] ## Excluded 'chopped.trimmed.fastq' as used by RNA Runner
        for afile in garbage:
            if os.path.isfile(afile): ## Check to see its a file from bowtie and not tophat mapped folder - Untested
                print("Deleting %s" % (afile))
                os.remove(afile)
            else:
                print("Skiiping cleanup, as its a directory %s" % (afile))

    else:
        pass
    
    end = time.time()
    runLog.write('Complete run time is %s' % (round(end-start,2)))
    runLog.close()
    
    print ('Complete run time is %s seconds' % (round(end-start,2)))
    
    ####################################################################################
    ####7. RNA-SEQ runner ######################################
    if indexBuilderStep     == 1:
        genoIndex = indexBuilder(genoFile)
    elif indexBuilderStep   == 0:
        if genoIndexHS:
            ## User provided index value 
            ## is not None or False; use
            ## provided index
            genoIndex = genoIndexHS
            pass
        else:
            ## Use index from the
            ## previous runs
            genoIndex = '%s/index/%s' % (os.getcwd(),genoFile.rpartition('/')[-1].rpartition('.')[0]) ## Can be merged with genoIndex from earlier part if we use bowtie2 earlier  
            print("This index will be used:%s" % (genoIndex))
        # sys.exit()
    else:
        print("Please check the 'indexBuilderStep' settings it should be 0/1")
        print("System will exit now")
        sys.exit()
       
    if spliceMapperStep == 1:
        rawInputs = [(i[0],'chopped.trimmed.fastq',i[2]) for i in register]     ## lib/filename, ext, nthreads
        print('Index being used: %s' % (genoIndex))
        splicedMapper(rawInputs,genoIndex,gtfFile)
    
    if stringTieStep == 1:
        rawInputs       = [(i[0],'hisat2.sorted.bam',i[2],1) for i in register] ## 'tophat.map' is a folder
        PPBalance(stringTie,rawInputs)
    
    if stringMergeStep == 1:
        rawInputs       = [(i[0],'stringtie.gtf',i[2]) for i in register]       ## 'cufflink.out' is a folder
        finalAssembly   = stringMerge(rawInputs,gtfFile,genoFile)

    if stringCountStep == 1:
        rawInputs       = [(i[0],'hisat2.sorted.bam',i[2],2) for i in register] ## 'cufflink.out' is a folder
        PPBalance(stringTie,rawInputs)

    if stringQuantStep == 1:
        rawInputs       = [(i[0],'stringtie.quant.gtf') for i in register]      ## 'cufflink.out' is a folder
        finalAssembly   = stringQuant(rawInputs)

if __name__ == '__main__':
    
    #### Assign Cores
    if ncores == 0:
        nproc = int(coreReserve(ncores))
    else:
        nproc = int(ncores)

    main(sampleInfo)
    print("\nScript finished sucessfully - you owe a beer !_! to Atul\n")
    
    sys.exit()

### CHANGELOG ---------------------------------------------------
## Version-1 preprocessing
##Screw losers i.e. Sa ta- I will write myself

##v02 -> v05
## First working copy of the script

## v05 -> v06
## Automatic resolve path for raw fasta file - so no need to copy in working directory
## Get tha adapter sequences for each library from master DB
## Get the max and min tag len from master DB
## Fixed maxTagLen issue - it was being swapped with minTagLen

## v06 -> v07
## Added RNA seq analysis till cuffmerge

## v07 -> v01 of seqAnalyze
## Completed the pipeline with most default parameters

## v07 -> Rocketv2
## Fixed tool Path for chopping module
## Organized the settings
## Added FASTQC analysis after trimming and chopping
## Fixed missing argumnet required - error for preprocessing graphs

##V2.0 -> v2.2
## Added option to generate assemblies based on just the transcriptome i.e. by excluding GTF file
## Added file check at cleanup step so as to avoid error while deleting bowtie map files which overlap with tophat map directories

## v2.2 - v2.3 [Major changes, Unstable] Feb 27
## Added a local mode to pre-process files witout using information from server. This is useful for cases when pre-processings differes from standard-lab procedures
## Paired end capabilities added
## Mismatches assigned based on MaxLength parameter for charts
## Need to test server mode??? i.e. when parameters are fetched from server

## v2.3 -> v2.4 [Stable]
## Fixed a minor filename error in quality check module
## Fixed input file  extension in all modules for single end- needs testing for paired end

## v2.4 -> v2.5
## Added function to read from summaryInfo file with three mandatory coulmns - id,rep and group
## Added function to club files either based on replicates or user defined group

##v2.5 -> v2.6
## Fixed a small error that stopped cuffquants module working, due to wrong structure of inputs - Trivial fix, didn't effected any results

## v2.6 -> v2.7
## Modified sample info file reading, first column in lib_id/lib_num , filename/lib_code, rep, group - In case of server mode lib_id is used and in case of local mode filename is used
## Added support for server-based paired end analysis

## v2.7b -> v2.8a
## Remove white space from lib_code while reading from a sampleinfo file
## Added full path to index
## Fixed bug in splice mapper, it was missing functionality if referenceGTF is "F"

## v2.8a -> v3.0 [Major][stable]
## Fixed trimming for paired end reads where most reverse reads were being thrown out. The reason being these shorter reads have same seqeunce in both left and right, so trimmomatic retains only one.
## .... see KeepBothReads option of trimmomatic

## v3.0 -> v3.1[major][stable]
## Fixed tophat mapping bug for stranded reads. Tophat was not using parameter libtype
## Old index can now be used

## v3.3 -> v3.4[major][stable]
## Cufflinks parallelized
## Cuffquant parallelization turned on
## Fixed PPBalance bug allocating process that require more then specified cores

## v3.4 -> v3.5
## Unkwnon changes while in India
## Update the TOPHAT version

## v3.5 -> v3.6
## Added "--max-bundle-frags" to cufflinks

## v3.6 -> v4.0b
## updated index to be made by hisat2
## Splice sites file generated before running script using: hisat2_extract_splice_sites.py genes.gtf > splicesites.txt
## added a function to sort bamfile to sam; this can be arallelized by providing list of libnames with bam extentsion

## v4.0b --> v4.0b1 [01/04/2019]
## stringTie function added
## added stringMerge, stringCount and stringQuant function
## removed cufflinks, cuffmerge, cuffquants and cuffnorm functions
## Note use AnnotationDBI package to add annotations: https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

## v4.0b1 --> v4.0b2 [01/04/2021]
## Moved all settings to a seprate YAML file
## Generalized paths for tools, including fastqc and trimmomatic
## removed redundant parameters "minLen" and "maxLen"
#### these were included as default values if the Meyers DB doesn't have these in lib info
## added library-specific reads stats files for 'stats_writer'

## v4.0b2 -> v4.0b3 [NEXT RELEASE]
## [Replace Bowtie1 for quality charts with HiSat2; this will remove requirement of older, redundant indexes]
#### [See Bowtie1 and Hisat outputs to make sure HiSat is providing required information]
## [Add index integrity checker from PHASIS]
## [UPDATE Stringtie and HiSat versions]


### Future fixes -----------------------------------------------
## Which quality scores used for sequencing: phred33 or phred66
## Adpaters need to selected or specified at begning of preprocessing
## Xlabels for abundance graph needs to be fixed
## Add splitting of paired end read capability - by adding script given to Ayush
## Improve naming - Extension added to end and before fastq after every step. This will improve deletion of specific files
## Add functionality of inputting lib names and replicates info from sample information file
## EC2 for running this script: trimmomatic is compute-intensive, 
#### Bowtie and HiSat2 mapping is compute and memory intensive 
#### (uses memory for page caching on 240 MB of 30.0 GB was free in test when mapping one paired end library)
#### StringTie is memory intensive (uses page cache, 9.4/32.0 GB free when running on two paired end libs)

## Sample info file format [Local]
# lib_id  filename  rep group tissue  stage
# 1   SRR531901   1   1   seed    4
# 2   SRR531902   1   1   seed    4
# 3   SRR531903   1   1   seed    4
# 4   SRR940300   2   1   seed    18
# 5   SRR531910   2   1   seed    18
# 6   SRR531912   2   1   seed    18

## Sample info file format [server]
# lib_id   code               rep     group
# 5170     Asp_le  Asp_le       1     1
# 5110     Asp_0_5_ant_budr     2     2
# 5158     Asp_0_5_ant_bu       3     2
# 5168     Asp_0_5_ant_b        3     2
# 5111     Asp_1_ant_budr       4     3
# 5159     Asp_1_ant_bu         5     3
# 5169     Asp_1_ant_b          5     3