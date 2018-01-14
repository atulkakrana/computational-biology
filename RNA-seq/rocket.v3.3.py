#!/usr/bin/python3

## Requires fastQC, trimmomatic and tally in $PATH variable
## Script written by ATUL Kakrana: kakrana@udel.edu

## Run : python3 ScriptName.py

## IMPORTS #####################################
import sys,os,re,time,timeit,csv,glob,string
import shutil,datetime,operator,subprocess,multiprocessing,matplotlib
import itertools as it
from multiprocessing import Process, Queue, Pool
import mysql.connector as sql
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

## PRE-PROCESSING SETTINGS ####################
Local           = 1                               ## 1: Local 0: Remote | NOTE: Local Analysis: Requires - maxLen,minLen,maxReadLen and adpaters.fa and libraries
                                                  ## Requires SampleInfo file with sampleNum, Filename/Lib code, reps, group
genomeDB        = 'None'                            ## [Server mode]For Bowtie index path used for mapping for graphs
## ROCKET ####################################
genoFile        = "/alldata/Genomic/Asparagus/UGA/v3/GenomicSEQ/AsparagusCHR_V1.1_modified2.fa"
gtfFile         = "AsparagusCHR_V1.1.evmRun2.withintrons.mod.gtf"
genoIndex       = "Aspa.v3.index"               ## If index is not in $ALLDATA i.e. local analysis, then specify bowtie1 index here for pre-processing. For Seq-analysis a Bowtie2 index  will be made using 'indexBuilderStep'
sampleInfo      = "sampleInfoALL.txt"           ## [mandatory] Tab delimted file with three mandatory columns - num (sample numbers), id (filename,library id), rep (same number if replicates
                                                ## And one optional columns group (sample grouped for edgeR analysis). See end of code for format.


referenceGTF    = 'T'                   ## [optional] T: True - use reference gtf file for merging assembling and annotation | F: Do not use GTF file and report transcripts based on transcriptome
                                        ## http://plants.ensembl.org/info/website/ftp/index.html OR USE gffread my.gff3 -T -o my.gtf (GFFREAD IS PART OF CUFFLINKS/TUXEDO)
libType         = 1                     ## [mandatory] From cuffNorm manual 1) 0 = fr-unstranded 2) 1 = fr-firststrand 3) 2 = fr-secondstrand
seqType         = 1                     ## [mandatory] 0: Single End; 1:Paired end (requires splitted reads - see fastq dump --split-reads for lib/or custom script)

groupBy         = 'R'                   ## [mandatory]   R: Group Samples by replicates, G: By user specified groups in sampleInfo 'group' column

hardMinTagLen   = 'Y'                   ## [server] Override Chopping values from server 
userMinTagLen   = 35                    ## [server] Used if 'hardMinTagLen' is ON
userMaxTagLen   = 140                   ## [server] Used if 'hardMinTagLen' is ON

## PRE_PROCESSING - OPTIONAL STEPS [Value: 0/1] ###############
QCheckStep       = 0                    ## Optional -Performs preliminary QC

## PRE_PROCESSING - REQUIRED STEPS [Value: 0/1] ##############
trimLibsStep     = 1                    ## Trim fastq files
preProGraphsStep = 1                    ## Generates before chopping graphs
chopLibsStep     = 1                    ## Chops adapter trimmed files
fastQ2CountStep  = 1                    ## Converts chopped to tag count
mapperStep       = 1                    ## Maps final chopped files and generates graphs
summaryFileStep  = 1                    ## Generates mapped summary - Never tested in Rocket, actually imported from prepro
cleanupStep      = 1                    ## Final cleanup

## SEQ-ANALYSIS - REQUIRED STEPS [Value: 0/1] ##############
indexBuilderStep = 0                    ## Build index for all the mappings
spliceMapperStep = 0                    ## Real Tophat mapping step
cuffLinksStep    = 0                    ## Cufflinks
cuffMergeStep    = 0                    ## Merge GTF and new assembly
cuffQuantStep    = 0
cuffNormStep     = 0

## ADVANCED SETTINGS #######################
minLen          = 35                    ## [mandatory] Min length of tag allowed
maxLen          = 140                    ## [mandatory] Max length of the tag allowed. Based on maxLen mismatches are allowed for mapping
unpairDel       = 1                     ## [Only for paired end analysis] 0: Retain unpaired read files after trimming 1: Delete these files

numProc         = 56                    ## [developer]  Coarse grain PP [0: Maximize parallel processing | [1-64]: Number of Cores]
nthread         = 16                    ## [developer]  Fine grain PP
maxReadLen      = 1000                  ## [developer]  Max allowed unchopped read length for graph generation

masterDB        = 'master'              ## [server]
dataServer      = 'raichu.ddpsc.org' ## [server]

## TOOL/FILE PATH ###################################
adapterFileSE   = '/data1/homes/kakrana/tools/Trimmomatic-0.32/adapters/TruSeq-SE.fa'      ## [mandatory] Sequence adapter file in FASTA format - Trimmomatic has files for different kits - Copy from there
adapterFilePE   = '/data1/homes/kakrana/tools/Trimmomatic-0.32/adapters/TruSeq-PE.fa'      ## You can try merged file with SE and PE adapter too, it gave best results to Atul, because adapters are sometimes not in pairs
tophat          = '/home/kakrana/tools/tophat-2.0.12/tophat2'                               ## [mandatory]
cufflinks       = '/home/kakrana/tools/cufflinks-2.2.1/cufflinks'                           ## [mandatory]
cuffmerge       = '/home/kakrana/tools/cufflinks-2.2.1/cuffmerge'                           ## [mandatory]
cuffquant       = '/home/kakrana/tools/cufflinks-2.2.1/cuffquant'                           ## [mandatory]
cuffnorm        = '/home/kakrana/tools/cufflinks-2.2.1/cuffnorm'                            ## [mandatory]
#################################################


## PREPARE ##################################

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

## PRE-PRO ##################################

## Output: "libName_fastqc"  
def QCheck(aninput):
    print(aninput)
    lib,ext,nthread,infile = aninput
    print('****Checking quality of %s library****' % (lib))    
    toolPath = "%s/svn/Tools/FastQC/fastqc" % (os.getenv('HOME'))

    outDir = "%s" % (os.getcwd())
    print(outDir)
    x = subprocess.Popen("%s --outdir=%s %s" % (toolPath,outDir,infile),shell=True)
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

## Output: "libName.trimmed.fastq"
def trimLibs(aninput):
    print(aninput)
    lib,ext,nthread,infile,adp_5p,adp_3p,minTagLen = aninput
    # print (aninput)
    print('\n****Trimming %s library with min length %s****' % (lib,minTagLen))
    toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    
    #### LOCAL ##########################
    if Local == 1: 

        ## Single End ###################
        if seqType == 0:
            adapterFile = adp_5p ## For local analysis adp_3p is actually the adapterfile - see main
            trimmedFile = '%s.trimmed.%s' % (lib,ext) ## Output
            
            trimLog     = '%s.trim.log' % (lib) ## Output
            trimSumm    = '%s.trim.summ.txt' % (lib)
            
            retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33","-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:%s:2:30:10" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (minTagLen)])
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
            
            retcode = subprocess.call(["java", "-jar", toolPath, "PE", "-phred33", "-threads", nthread, infile1, infile2, trimmedFileP1,trimmedFileU1,trimmedFileP2,trimmedFileU2, "ILLUMINACLIP:%s:2:30:10:8:TRUE" % (adapterFile), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:%s" % (minTagLen)])
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
            retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-trimlog",trimLog,"-threads", nthread, infile, trimmedFile, "ILLUMINACLIP:./%s_adapter.fa:2:30:10" % (lib), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (minTagLen)])
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
            
            retcode = subprocess.call(["java", "-jar", toolPath, "PE", "-phred33", "-trimlog",trimLog,"-threads", nthread, infile1, infile2, trimmedFileP1,trimmedFileU1,trimmedFileP2,trimmedFileU2, "ILLUMINACLIP:./%s_adapter.fa:2:30:10:8:TRUE" % (lib), "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:10", "MINLEN:%s" % (minTagLen)])
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

## Output: "libName.chopped.fastq" 
def chopLibs(aninput):
    ''' Reverse read set of paired end lib is chopped from right end (5' in actuality) as done by using -
    if we do chopping with trimming using PE mode it is still chopped the same way '''

    print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('****Chopping %s library to length %s****' % (lib,maxTagLen))

    trimmedInFile = '%s.%s' % (lib,ext)
    choppedOutFile = '%s.chopped.%s' % (lib,ext)
    print("\n")
    toolPath = "%s/svn/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar" % (os.getenv('HOME'))
    
    if seqType == 0:
        retcode = subprocess.call(["java", "-jar", toolPath, "SE", "-phred33", "-threads", nthread, trimmedInFile, choppedOutFile, "CROP:%s" % (maxTagLen)])
    elif seqType == 1:
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
    print('\n****Converting %s.%s file to tag count****\n' % (lib,ext))
    infile = '%s.%s' % (lib,ext)
    outfile = '%s.%s.processed.txt' % (lib,ext.replace(".fastq",""))
    print("This is outfile:%s" % (outfile))
    # sys.exit()
    retcode = subprocess.call(["tally", "-i", infile, "-o", outfile, "--nozip", "-format","%R%t%X%n"])
    if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n**** Conversion to tag count format for %s complete****' % (infile) )
    else:
        print("Something wrong happened while converting to %s library tag count - Debug for reason" % (lib))
        sys.exit()

## ASSEMBLY ###################################

## Output: "libName.map"
def mapper(rawInputs,mode):
    
    ## For all libs - One by one
    for aninput in rawInputs:
        print('\nInput:',(aninput))
        lib,ext,nthread,maxTagLen = aninput
        
        # Resolve index path #################################
        if Local == 0:
            cur = con.cursor()
            cur.execute("SELECT bowtie_index_path FROM master.genome_db WHERE genome_db like '%s'" % (genomeDB)) ##bowtie_index_path
            genoIndexPath = cur.fetchall()
            
            genoIndexPrePro = genoIndexPath[0][0].replace('$ALLDATA', '/alldata') ### Index file

        else:
            genoIndexPrePro = genoIndex

        print ('Genomic index being used for mapping: %s\n'% (genoIndexPrePro))
        #genoIndex = 'ASPARAGUS_UGA1_genome' ## Test
        
        ### Prepare ###########################################
        inFile = '%s.%s' % (lib,ext)
        print ('Processing %s for mapping to genome' % (inFile))
        fastaFile = tagCount2FASTA(inFile,'N') ## Unique reads to FASTA format 

        mapFile = ('./%s.%s.map' % (lib,ext.rpartition('.')[0]))
        print(genoIndexPrePro,inFile,fastaFile,mapFile)
        
        ## Map to index ##########################################
        print ('Mapping %s processed file to genome' % (lib))
        nproc2 = str(nproc)

        if int(maxTagLen) > 60:
            mismat = str(2)
        elif int(maxTagLen) <= 60 and maxTagLen > 40:
            mismat = str(1)
        elif int(maxTagLen) <= 40:
            mismat = str(0)
        else:
            pass
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])
        
        ## Bowtie 1 - So as to be compatible with current indexes
        retcode = subprocess.call(["bowtie","-f","-n",mismat,"-p", nproc2,"-t" ,genoIndexPrePro, fastaFile, mapFile])
        
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
        ###mappedList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95, 166, 278, 869, 735, 1515, 7389, 694, 122, 600, 460, 39, 33, 26, 86]
        ###mappedAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 594, 805, 1025, 5581, 4444, 4992, 20590, 1714, 805, 732, 595, 476, 392, 299, 180]
        ###
        ###allTagsList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 126, 220, 370, 1103, 954, 1886, 9010, 1012, 274, 140, 105, 913, 825, 835, 644]
        ###allAbunList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 660, 944, 1217, 6338, 5015, 5567, 23419, 2145, 1029, 858, 709, 672, 658, 416, 294]

        ## Plot
        charts(lib,ext,mappedList,mappedAbunList,allTagsList,allAbunList,mode) ## Mode 1 - Preprocess graphs 2: processed files graphs
        
    
    return None

def indexBuilder(genoFile):
    print ("\n**Deleting old index 'folder' !!!!!!!!!!!***\n")
    print("If its a mistake cancel now by pressing ctrl+D and continue from index step by turning off earlier steps- You have 30 seconds")
    time.sleep(30)
    shutil.rmtree('./index', ignore_errors=True)
    os.mkdir('./index')
    
    genoIndex = '%s/index/%s' % (os.getcwd(),genoFile.rpartition('/')[-1].rpartition('.')[0]) ## Can be merged with genoIndex from earlier part if we use bowtie2 earlier
    print('**Creating index of cDNA/genomic sequences:%s**\n' % (genoIndex))
    retcode = subprocess.call(["bowtie2-build", genoFile, genoIndex])
    return genoIndex

def splicedMapper(rawInputs,genoIndex,gtfFile):
    
    for aninput in rawInputs:
        print('\nInput:',(aninput))
        lib,ext,nthread = aninput
        
        ### Prepare ###########################################
        if seqType == 0:
            inFile = '%s.%s' % (lib,ext)
            print ('Mapping %s for to genome for RNA-Seq analysis' % (inFile))
            mapFile = ('./%s.tophat.map' % (lib))
            print(genoIndex,inFile,mapFile)
        
        else:
            inFile1 = '%s.chopped.pair_1.trimmed.fastq' % (lib)
            inFile2 = '%s.chopped.pair_2.trimmed.fastq' % (lib)
            inFile = [inFile1,inFile2] ## Just for print statement
            mapFile = ('./%s.tophat.map' % (lib))
            print(genoIndex,inFile1,inFile2,mapFile)
        
        ## Map to index ##########################################
        print ('Mapping %s processed file to genome' % (lib))
        nproc2 = str(nproc)
        
        ## Bowtie2 for future - Needs retest for speed before switching
        #retcode = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread,"-f", genoIndex,fastaFile,"-S",mapFile])

        #### Select
        if libType == 0:
            atype = 'fr-unstranded'
        if libType == 1:
            atype = 'fr-firststrand'
        if libType == 2:
            atype = 'fr-secondstrand'

        if referenceGTF == 'T':
            if seqType == 0:
                retcode = subprocess.call([tophat,"-p",nproc2,"-G",gtfFile,"--library-type",atype,"-o",mapFile,genoIndex,inFile])
            else:
                retcode = subprocess.call([tophat,"-p",nproc2,"-G",gtfFile,"--library-type",atype,"-o",mapFile,genoIndex,inFile1,inFile2])
        
        elif referenceGTF == 'F':
            if seqType == 0:
                retcode = subprocess.call([tophat,"-p",nproc2,"--library-type",atype,"-o",mapFile,genoIndex,inFile])
            else:
                retcode = subprocess.call([tophat,"-p",nproc2,"--library-type",atype,"-o",mapFile,genoIndex,inFile1,inFile2])

        else:
            print("Please choose correct option for 'referenceGTF' and provide 'GTF' file if selected - Script will exit for now")
            sys.exit()

        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\nTophat mapping for %s complete' % (inFile) )
        else:
            print ("There is some problem with Tophat mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()

def cuffLinks(rawInputs):

    '''
    create a cufflink assembly for each libarary using tophat generated bam files
    '''
    
    print ('\n**Cufflink anlysis in progress**\n')
    
    # assemblyOut = open('assemblies.txt','w')
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

def cuffMerge(rawInputs,gtfFile,genoFile):
    
    '''
    Create a merged transcript from cufflink assemblies
    '''
    
    ## Make a list of lib-wise assemblies
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
    print('Input files: %s,%s,%s,%s\n' % (gtfFile,genoFile,nproc2,assemblyFile))
    if referenceGTF == 'T':
        print ('System:', cuffmerge, "-g", gtfFile, "-s", genoFile,"-p", nproc2, assemblyFile)
        retcode = subprocess.call([cuffmerge, "-g", gtfFile, "-s", genoFile,"-p", nproc2, assemblyFile])
    else:
        print ('System:', cuffmerge,"-s", genoFile,"-p", nproc2, assemblyFile)
        retcode = subprocess.call([cuffmerge,"-s", genoFile,"-p", nproc2, assemblyFile])

    if retcode == 0:
        print("Cuffmerge step complete\n\n")
    else:
        print("There is some problem with cuffMerge - System wil exit now\n")
        sys.exit()

    
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

def cuffNorm(rawInputs):

    print ('\n** Executing cuffNorms **')
    
    libs = []
    reps =[]
    for i in rawInputs:
        libs.append(i[0])
        reps.append(i[1])

    print('These are libs:', libs, '| These are reps',reps)
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

    ## Write a file with order of columns in results - See samples.table file too
    fh_out = open('./Final.cuffNorm.out/sampleHeader.txt','w')
    fh_out.write('%s\n' % ('\t'.join(sampleDict.keys()) )) 
    fh_out.close()
        
    return None

def compileResults():
    ## Combine 
    pass

############ STANDARAD FUNCTIONS ###############

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

def mappedStats(aninput,mode):
    '''Parse map file and collect statics for graph generation'''
    print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('\nCollecting statistics of matched reads for Lib:%s' % (lib))
    
    inFile = '%s.%s.map' % (lib,ext.rpartition('.')[0])
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

def tagCountStats(aninput,mode):

    '''Get stats for all the reads from tagCount file'''
    
    #print(aninput)
    lib,ext,nthread,maxTagLen = aninput
    print('\nCollecting statistics of total reads for Lib:%s' % (lib))
    
    inFile = '%s.%s' % (lib,ext)
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
    ''' Write a lib-specific summary file'''
    lib,minTagLen,maxTagLen = aninput
    print (aninput)

    ### Read the library temp files with counts and abundances before and after processing
    print("Reading stats from: %s_allBefore.temp" % lib)
    fh_before = open("%s_allBefore.temp" % lib,'r')
    aread = fh_before.read().split('\n')
    allTags,allAbun = aread
    allTagsList,allAbunList = list(map(int,allTags.split(','))),list(map(int,allAbun.split(',')))

    # print("allTagsList:",allTagsList,"\nallAbunList:",allAbunList)

    fh_after = open("%s_allAfter.temp" % (lib),'r')
    aread2 = fh_after.read().split('\n')
    allTags2,allAbun2 = aread2
    allTagsList2,allAbunList2 = list(map(int,allTags2.split(','))),list(map(int,allAbun2.split(',')))

    fh_map_before = open("%s_mappedBefore.temp" % (lib), 'r')
    aread3 = fh_map_before.read().split('\n')
    mapped,mappedAbun = aread3
    mappedList,mappedAbunList = list(map(int,mapped.split(','))),list(map(int,mappedAbun.split(',')))
    
    fh_map_after = open("%s_mappedAfter.temp" % (lib),'r')
    aread4 = fh_map_after.read().split('\n')
    mapped2,mappedAbun2 = aread4
    mappedList2,mappedAbunList2 = list(map(int,mapped2.split(','))),list(map(int,mappedAbun2.split(',')))

    ### Prepare to write
    summFile = "%s_chopinfo.txt" % lib
    fh_out = open(summFile,'w')
    fh_out.write("Date - %s | Genome - %s\n" % (time.strftime("%d/%m/%Y"),genomeDB))
    fh_out.write("Lib-%s\tBeforeProcessing\tAfterProcessing\tMappedBeforeProcessing\tmappedAfterProcessing\n" % (lib))

    # print("allTagsList length:",len(allTagsList),"\nallAbunList length:",len(allAbunList))
    indexList = [i for i,x in enumerate(allTagsList) if x != 0]
    # print ('indexList:',indexList)
    minLen = min(indexList)
    maxLen = max(indexList)

    indexList2 = [i for i,x in enumerate(allAbunList) if x != 0]
    # print ('indexList2:',indexList)
    minLenAbun = min(indexList2)
    maxLenAbun = max(indexList2)
    # print("Allowed min len:%s | Allowed max len:%s" % (minTagLen,maxTagLen))
    # print("allTags min len:%s | allTags max len = %s | allAbun max len:%s | allAbun max len:%s\n" % (minLen,maxLen,minLenAbun,maxLenAbun))

    countsAllBefore = list(allTagsList[minLen:maxLen+1])
    countsAllAfter = list(allTagsList2[minLen:maxLen+1])
    abunAllBefore = list(allAbunList[minLen:maxLen+1])
    abunAllAfter = list(allAbunList2[minLen:maxLen+1])
    print(countsAllBefore,countsAllAfter,abunAllBefore,abunAllAfter)

    countsMappedBefore = list(mappedList[minLen:maxLen+1])
    countsMappedAfter = list(mappedList2[minLen:maxLen+1])
    abunMappedBefore = list(mappedAbunList[minLen:maxLen+1])
    abunMappedAfter = list(mappedAbunList2[minLen:maxLen+1])

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

##### DOWNSTREAM ###############################

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

################# OBSELETE #####################

def cuffMergeQuant2(rawInputs,genoFile):
    
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
    -u/--multi-read-correct  Tells Cuffquant to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome.
    '''
    
    command = "%s -p %s -u %s %s" % (cuffquant,nproc2,mergedAssembly,mappedFilesPath) ## Had problem supplying multiple files seprated with space using normal switch based approach approach
    print(command)
    retcode =subprocess.call(command,shell=True)
        
    return None

############## MAIN ###########################
def main(sampleInfo):
    
    start = time.time() 
    runLog = open('%s_run.log' % (datetime.datetime.now().strftime("%m_%d_%H_%M")), 'w')
    
    #### 0. Initialize Input Register ###################################################
    register =[] ## Holds values for all different steps - master feeder
    libs,rep = sampleInfoRead(sampleInfo)
    print('\n\n**Total %s libraries provided for pre-processing**\n' % (len(libs)))
    
    if Local == 1: ## Local analysis
            for i,k in zip(libs,rep):
                ext = 'fastq'
                minTagLen = minLen
                maxTagLen = maxLen
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
        

        if seqType == 0: ## SingleEnd
            rawInputs = [(i[0],'trimmed.processed.txt',i[2],i[7]) for i in register] ## ## Library/Filename, extension, max tag len
        else: ## PairedEnd
            inputsR = [(i[0],'pair_1.trimmed.processed.txt',i[2],i[7]) for i in register] ## ## Library/Filename, extension, max tag len
            inputsL = [(i[0],'pair_2.trimmed.processed.txt',i[2],i[7]) for i in register] ## ## Library/Filename, extension, max tag len
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
        print ("**Cleaning temp files**")
        garbage = [file for file in os.listdir('./') if file.endswith (('.map','trimmed.processed.txt','.processed.fa'))]
        for file in garbage:
            print("Deleting %s" % (file))
            os.remove(file)
    else:
        pass

    ######################################################################################    
    #### 3. Chop trimmed files ###########################################################
    
    ## SE and PE on Local and Server Ready
    if seqType == 0: ## SingleEnd
        rawInputs = [(i[0],'trimmed.fastq',i[2],i[7]) for i in register] ## ## Library/Filename, 'trimmed.fastq', max tag len
    else: ## PairedEnd
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
        rawInputs = [(i[0],'chopped.trimmed.processed.txt',i[2],i[7]) for i in register] ## ## Library/Filename,extension, nthread
    else: ## PairedEnd
        inputsR = [(i[0],'chopped.pair_1.trimmed.processed.txt', i[2],i[7]) for i in register] ## ## Library/Filename, extension, nthread, max tag len
        inputsL = [(i[0],'chopped.pair_2.trimmed.processed.txt', i[2],i[7]) for i in register] ## ## Library/Filename, extension, nthread, max tag len
        rawInputs = inputsL+inputsR
    
    if mapperStep == 1:
        print('\n**Mapping chopped files to the genome index**\n')
        # con = ConnectToDB(dataServer,0)
        mapper(rawInputs,2)
    else:
        print("\n**Mapping Step skipped as selected - What are you intentions?\n")
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
    if indexBuilderStep == 1:
        genoIndex = indexBuilder(genoFile)
    elif indexBuilderStep == 0:
        genoIndex = '%s/index/%s' % (os.getcwd(),genoFile.rpartition('/')[-1].rpartition('.')[0]) ## Can be merged with genoIndex from earlier part if we use bowtie2 earlier  
        print("This index will be used:%s" % (genoIndex))
        # sys.exit()
    else:
        print("Please check the 'indexBuilderStep' settings it should be 0/1")
        print("System will exit now")
        sys.exit()
       
    if spliceMapperStep == 1:
        rawInputs = [(i[0],'chopped.trimmed.fastq',i[2]) for i in register] ## lib/filename, ext, nthreads
        print('Index being used: %s' % (genoIndex))
        splicedMapper(rawInputs,genoIndex,gtfFile)
    
    if cuffLinksStep == 1:
        rawInputs = [(i[0],'tophat.map',i[2]) for i in register] ## 'tophat.map' is a folder
        assemblyOut = cuffLinks(rawInputs)
    
    if cuffMergeStep == 1:
        rawInputs = [(i[0],'cufflink.out',i[2]) for i in register] ## 'cufflink.out' is a folder
        finalAssembly = cuffMerge(rawInputs,gtfFile,genoFile)
    
    if cuffQuantStep == 1:
        rawInputs = [(i[0],'tophat.map',i[2]) for i in register] ##  
        for i in rawInputs:
            cuffQuant(i)
        # PPBalance(cuffQuant,rawInputs)
    
    if cuffNormStep == 1:
        rawInputs = [(i[0],i[8],'cuffquant.out') for i in register]
        cuffNorm(rawInputs)

################ Execute #######################
if __name__ == '__main__':
    
    #### Assign Cores
    if numProc == 0:
        nproc = int(multiprocessing.cpu_count()*0.95)
    else:
        nproc = int(numProc)
    
    #### Execute modules

    main(sampleInfo)
    print("\n\n-------Script finished sucessfully - CHEERS!!!! - You owe a beer !_! to Atul--------\n")
    
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
## Fixed trimming for paired end reads where most reverse reads were being thrown out. The reason being these shorter reads have same seqeunce in both left and right, so timmomatic retails only one.
## .... see KeepBothReads option of trimmomatic

## v3.0 -> v3.1[major][stable]
## Fixed tophat mapping bug for stranded reads. Tophat was not using parameter libtype
## Old index can now be used

### Future fixes -----------------------------------------------
## Which quality scores used for sequencing: phred33 or phred66
## Adpaters need to selected or specified at begning of preprocessing
## Xlabels for abundance graph needs to be fixed
## Add splitting of paired end read capability - by adding script given to Ayush
## Improve naming - Extension added to end and before fastq after every step. This will improve deletion of specific files
## Add functionality of inputting lib names and replicates info from sample information file
## Parallelize cufflinks


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
