#!/usr/local/bin/python3
##sPARTA: small RNA-PARE Target Analyzer Version 
##Updated: version-1.02 6/12/2014

#### PYTHON FUNCTIONS ##############################
import argparse
import sys,os,re,time,glob
import subprocess, multiprocessing
import shutil
import operator
from multiprocessing import Process, Queue, Pool
from operator import itemgetter
import pyfasta
import rpy2.robjects as robjects
from scipy import stats
import numpy
import datetime
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
RStats = importr('stats')

nthread = 6
#### USER SETTINGS ################################
parser = argparse.ArgumentParser()
parser.add_argument('-gffFile',  default='', help='GFF file for the species '\
    'being analyzed corresponding to the genome assembly being used')
parser.add_argument('-genomeFile', default='', help='Genome file in FASTA '\
    'format')
parser.add_argument('-featureFile', default='', help='Feature file in FASTA '\
    'format')
parser.add_argument('-genomeFeature', required=True, help='0 if prediction is '\
    'to be done in genic region. 1 if prediction is to be done in intergenic '\
    'region')
parser.add_argument('-miRNAFile', default='', help='FASTA format of miRNA '\
    'sequences')
parser.add_argument('-tarPred', nargs='?', const='H', help='Mode of target '\
    'prediction. H for heuristic. E for exhaustive. H is default if no mode '\
    'is specified')
parser.add_argument('-tarScore', nargs='?', const='S', help='Scoring mode '\
    'for target prediction. S for seedless. N for normal. S is default if '\
    'no mode is specified')
parser.add_argument('-libs', nargs='*', default=[], help='List of PARE '\
    'library files in tag count format. Data can be converted into tag '\
    'count format using')
parser.add_argument('-tagLen', default=20, help='Minimum length of PARE tag, '\
    'tags longer than tagLen will be chopped to the specified length. 20 is '\
    'default')
parser.add_argument('--tag2FASTA', action='store_true', default=False, help=
    'Convert tag count file for PARE libraries to FASTA files for mapping')
parser.add_argument('--map2DD', action='store_true', default=False, help=
    'Map the PARE reads to feature set')
parser.add_argument('--validate', action='store_true', default=False, help=
    'Flag to perform the validation of the potential cleave sites from '\
    'miRferno')
parser.add_argument('--repeats', action='store_false', default=True, help=
    'Flag to include PARE reads from repetitive regions')
parser.add_argument('--noiseFilter', action='store_false', default=True,
    help='Flag to include all PARE validations with p-value of <=.5, '\
    'irrespective of the noise to signal ratio at cleave site and category '\
    'of PARE read')
parser.add_argument('-accel', default='Y', help='Y to use '\
    'balanced multiple process scheme or else specify the number of '\
    'processors to be used. Y is default')

### Developer Options ###
parser.add_argument('--generateFasta', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('--fileFrag', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('--indexStep', action='store_true', default=False,
    help=argparse.SUPPRESS)
parser.add_argument('-splitCutoff', default=20, help=argparse.SUPPRESS)
parser.add_argument('-maxHits', default=30, help=argparse.SUPPRESS)
parser.add_argument('--cat4Filter', action='store_false', default=True,
    help=argparse.SUPPRESS)

args = parser.parse_args()

## Various checks for dependencies within command line arguments
# If either gff or genome file is given without the other and featureFile
# is not given, exit.
if(((args.gffFile and not args.genomeFile) or (args.genomeFile and not
        args.gffFile)) and (not args.featureFile)):
    print("gffFile and genomeFile must both be given to create feature set")
    exit()

# If the user input both a genome and feature file, exit as both cannot be
# supplied for proper execution
if(args.genomeFile and args.featureFile):
    print("genomeFile and featureFile cannot both be supplied for execution")
    exit()

# If gffFile and genomeFile are given turn on extraction, frag and index steps
# must be set on
if(args.gffFile and args.genomeFile):
    args.generateFasta = True
    args.fileFrag = True
    args.indexStep = True

# If featureFile is given, frag and index steps must be set on
if(args.featureFile):
    # If featureFile is given and gffFile is given, give a warning letting
    # user know the gffFile will be ignored and the input fasta file may
    # have been intended as a genomeFile
    if(args.gffFile):
        print("Warning: You have input a gffFile but input a FASTA file as "\
        "the featureFile. If you intended for this to be used in conjunction "\
        "with the gff file to create a feature file, please press 'ctrl+c' "\
        "to cancel the execution and rerun with the FASTA file under the "\
        "argument 'genomeFile'. If this is in fact the feature file, allow "\
        "sPARTA to continue its execution.")
        time.sleep(10)
    args.fileFrag = True
    args.indexStep = True

# If indexStep is on and tarPred is off, turn tarPred and tarScore on
if(args.indexStep):
    if(not args.tarPred):
        args.tarPred = 'H'
    if(not args.tarScore):
        args.tarScore = 'S'

# If tarPred is on, then tarScore will default to S
if(args.tarPred and not args.tarScore):
    args.tarScore = 'S'

# If tarPred is on, then miRNAFile must be provided
if(args.tarPred and not args.miRNAFile):
    print("miRNA file must be given to perform target prediction")
    exit()

# If tag2FASTA is on, turn map2DD on
if(args.tag2FASTA and not args.map2DD):
    args.map2DD = True

# If tag2FASTA is on, then libraries must be defined
if(args.tag2FASTA and not args.libs):
    print("libs must be assigned to perform tag2FASTA")
    exit()

# If validate is on, then libraries must be input
if(args.validate and not args.libs):
    print("At least one library must be given to perfor the validate")
    exit()

# genomeFeature must be an integer
args.genomeFeature = int(args.genomeFeature)

################### STEPS #########################

###################################################
#### MPPP FUNCTIONS ###########################

#### Extract coordinates from GFF file -- Will update script with latest version of scripts
def extractFeatures(genomeFile,gffFile):
    
    fh_in = open(gffFile,'r')
    fh_in.readline() ## GFF version
    gffRead = fh_in.readlines()
    genome_info = [] ## List to hold all coordinates for which fasta will be fetched
    for i in gffRead:
        ent = i.strip('\n').split('\t')
        #print (ent)
        if ent[2] == 'gene':
            chrID = ent[0][3:]
            strand = ent[6].translate(str.maketrans("+-","WC"))
            geneName = ent[8].split(';')[0].split('=')[1]
            geneType = 'gene'
            #print(chrID,strand,geneName,ent[3],ent[4],geneType)
            genome_info.append((chrID,strand,geneName,int(ent[3]),int(ent[4]),geneType))
        
    #genome_info_inter = genome_info ###This list will also hold intergenics
    genome_info_sorted = sorted(genome_info, key=operator.itemgetter(0,1,3))
    genome_info.sort(key=lambda k:(k[0],k[1],k[3])) ### Sorting on basis of chr_id, strand and gene start for coord calulation and reduce GetFASTA time
    genome_info_inter = genome_info
    alist = []###list maintained to check if first gene on chromosome and strand shows up than intergenic is just the start of gene
    for i in range(0, int(len(genome_info))-1):###maybe for i in range(0, len(genome_info -1))
        #print (i)
        gene1 = (genome_info[i])
        gene2 = (genome_info[i+1])
        gene_type = 'inter' ###set to integentic by default
        #print(gene1,gene2)
        
        ##Remove/skip redundant genes with same start and end......What about Overlapping genes????
        if gene1[3] == gene2[3] and gene1[4] == gene2[4]:
            ##gene is same/overlapping consider next gene
            pass
        
        else:
            ##Calculate coordinates of intergenic regions
            if tuple(gene1[0:2]) not in alist:##Only chr_id and strand is checked. This is first gene on chromosome and strand intergenic region is from position1
                print ('Caching gene coords for chromosome: %s and strand: %s\n' % (gene1[0], gene1[1]))
                alist.append((gene1[0:2]))
                inter_start1 = 1
                inter_end1 = int(gene1[3])-1###1 nt before start of Gene1 a.k.a the first gene on chromosome in this case
                ##As two genes are read together, the upstream intergenic region gor gene2 must be calculated in same step
                inter_start2 = int(gene1[4])+1##From end of first gene of chromosome
                inter_end2 = int(gene2[3])-1###Till start of second gene
                
                if gene1[1] == 'w': ##The gene is on positive strand so upstream
                    inter_name1 = ('%s_up' % (gene1[2]))
                    inter_name2 = ('%s_up' % gene2[2])

                else: ##Its on negative strand
                    inter_name1 = ('%s_down' % (gene1[2]))
                    inter_name2 = ('%s_up' % gene1[2])
                genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                genome_info_inter.append((gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type

            else:
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]:###If chr_id and strands are equal than find intergenic. These are gene on same chromosme and strand
                    inter_start = int(gene1[4])+1###End of Gene 1
                    inter_end = int(gene2[3])-1 ###1 nt before start of gene 2
                    if gene2[1] == 'w': ##Positive strand
                        inter_name = ('%s_up' % (gene2[2]))
                    else:## reverse strand
                        inter_name = ('%s_up' % (gene1[2]))
                    genome_info_inter.append((gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type))
                
                else: ###That means gene1 is at end of one chromosome and gene 2 is begining of chromosome so we have to extract intergenic at end of one chromosome
                    inter_start = int(gene1[4])+1###End of gene1
                    #inter_end = chromo_dict[gene1[0]]###End of chromosome searched using chromosome id of gene1 from chromosome dictionary
                    #inter_end = '-'### Extract till end of chromosome
                    inter_end = '-'
                    if gene1[1] == 'w':##Positive strand end of chromosme
                        inter_name = ('%s_down' % (gene1[2]))
                    else: ##Negative strand first intergenic of chromosme
                        inter_name = ('%s_up' % (gene1[2]))                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type))##Chr_id, strand
    
    #print("\n\nThis is length of genome info inter append:%s" % (len(genome_info_inter)))
    ###Sort the list after adding intergenic regions on on basis of chr_id and strand that is essential while caching chromosme during slicing sequences
    genome_info_inter_sort = sorted(genome_info_inter, key=operator.itemgetter(0,1))
    
    ###Filter list to remove unwanted types like miRNA,tRNA,rRNA,snoRNA,snRNA, short or no intergenic
    
    gene_coords_file = './coords'####To check wheter coords are printed in chr_id and strand sorted or not
    coords_out = open(gene_coords_file, 'w')
    coords = []## List that will hold genes to fetch, this removes unecessary RNAs and also fix miRNA double entry i.e as gene and miRNA        
    
    if args.genomeFeature == 2: ## Both gene and inter
        for ent in genome_info_inter_sort:
            print(ent)
            if ent[4] == '-': ## End of chromosome
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                
            elif int(ent[4])-int(ent[3]) > 25:###If there is no intergenic regon b/w genes or too short than filter
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
    
    else:
        if args.genomeFeature == 0:
            genomeFilter = 'gene'
        elif args.genomeFeature == 1:
            genomeFilter = 'inter'
        for ent in genome_info_inter_sort:
            if (ent[5] == genomeFilter):
                #print(ent)
                if ent[4] == '-': ## End of chromosome
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                    
                elif int(ent[4])-int(ent[3]) > 25:###If there is no intergenic regon b/w genes or too short than filter
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
    
    print ("Number of coords in 'coords' list: %s" % (len(coords)))
    coords_out.close()    
    fh_in.close()
    
    return coords

### With Inbuilt fragFASTA functionality - Removes pyfasta dependency - Break coords list and prepare fragmented genome - Lazy file loading i.e. only required chromosomes are loaded and flushed - 
def getFASTA1(genomeFile,coords):
    fastaOut = './genomic_seq.fa'
    fh_out = open(fastaOut, 'w')

    ## Read Genome file - Need lazy loading or chunk reading - Else will fail on large genomes
    fh_in = open(genomeFile, 'r')
    genomeFile = fh_in.read()
    genomeList = genomeFile.split('>')[1:] ##First block is empty
    chromoDict = {} ##Chrid,seq as a tuple
    for i in genomeList:
        chromoInfo = i.partition('\n') ## [before split, sep, after split]
        chrid = chromoInfo[0].split()[0] ##>Chr1 CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2007-12-20
        #chrSeq = chromoInfo[2].rstrip('\n') ## Not works - Why?
        chrSeq = chromoInfo[2].replace("\n", "")
        #print (chrSeq)
        chromoDict[chrid] = [chrSeq]
        #print(chrid,chrSeq[0:50] )
    
    ### For every entry either gene or intergenic, use gene_info for just genes
    chromo_mem = []##used as memory to keep track of Chr+strand in memory, if its in here chromosome will not be read, untill its not here and appended
    for i in coords: 
        #print (i)
        gene = i[2]
        chr_id = i[0]
        strand = i[1]
        start = i[3]-1###Adjusted to pick from begining because when spliced it does not include starting position so starting from one poistion back
        end = i[4]
        #print('start:%s End:%s Chr:%s Strand:%s' % (start,end,chr_id,strand))
        
        ### Add chromosome and strand to list and fetch in memory, will be executed first time a unique chr_id and strand found
        ### as list is sorted on chr_is and strand
        if tuple(i[0:2]) not in chromo_mem: 
            chromo_mem.append(tuple(i[0:2]))   ##Append for first time
            print ("Reading chromosome:%s and strand: '%s' into memory to splice genes" % (i[0],i[1]) )
            chrKey = 'Chr' + i[0]
            chromo = str(chromoDict[chrKey])
            #print('Chromosome:',chromo)
            
            gene_seq = chromo[start:end]  ##Get first gene on that strand
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1]
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
         
        ###If chr_id is in chromo_mem i.e. its loaded in memory just splice your sequence      
        
        elif end == '-': ## Till end of chromosome
            print('Fetching gene %s to prepare FASTA file' % (gene))
            gene_seq = chromo[start:]##
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1]
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
            
        else:
            print('Fetching gene %s to prepare FASTA file' % (gene))
            gene_seq = chromo[start:end]##
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1]
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))

    time.sleep(10)
    fh_out.close()
    
    ## Intergenic splicing problem
    ### fasta file empty problem - Fixed
    
    return fastaOut

###To efficiently break fasta file you need to distribute size and complete sequences
def fragFASTA(FASTA):
    
    ####Purge older files
    shutil.rmtree('./genome', ignore_errors=True)
    os.mkdir('./genome')
    
    pattern = ".*\.[0-9].*\.fa" ## all_cdna.100.fa or all_cdna.1.fa
    print ("\n***Purging older files***")
    for file in os.listdir():
        if re.search(pattern,file):
            print (file,'is being deleted')
            os.remove(file)
    
    statInfo = os.stat(FASTA)
    filesize =round(statInfo.st_size/1048576,2)
    print('\n**Size of FASTA file: %sMB**' % (filesize))### Covert to MB
    
    if filesize <= args.splitCutoff: ## No need to split file too small
        fls = []
        fls.append(FASTA)
        print ('No fragmentation performed for file %s' % (fls))
        
    else: ###Split file
        ##Check if chromosome file or cDNA - Get Number of headers in file
        fh_in = open(FASTA, 'r')
        seq_count = fh_in.read().count('>')
        print('**Number of headers in file: %s**\n'% (seq_count))
        #if genome == 'N':
        if seq_count >= 30: ## Check file is chromosme based or contig/cDNA based 
            
            ##Calculate number of fragments to split in - Should be done automatically
            if filesize <= 3072:
                splitnum = str(args.maxHits)
            elif filesize > 3072 and filesize <= 5120:
                splitnum = str(round(args.maxHits*1.25))
            else:
                splitnum = str(round(args.maxHits*1.5))        
            
            print ("Size based fragmentation in process for '%s' file" % (FASTA))
            retcode = subprocess.call(["pyfasta","split", "-n", splitnum, FASTA])
            fls = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)] ## file list using regex
            #fls = glob.glob(r'%s.[0-9]{1-3}.fa' % (FASTA.split('.')[0])) ## fragmented file list ##
            print ('The fragments: %s' % (fls))
               
        
        else: ## chromosome file - Break file on number of fragments
            splitnum = str(seq_count) ## Number of files to split in, try to keep size close to 100MB
            if fragFasta == 'Y':
                print ("Header based fragmentation in process for '%s' file" % (FASTA))
                retcode = subprocess.call(["pyfasta","split", "-n", splitnum, FASTA])
            fls = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)]
            #fls = glob.glob(r'%s.[0-9]{1,3}.fa' % (FASTA.split('.')[0])) ## fragmented file list
            #os.chdir("../")
            
            print ('The fragments: %s' % (fls))
            
    ##Get back to working directory
    #os.chdir("../")
    #print(os.getcwd())   
    return fls

##Input miRNAs from given file else from server
def miRinput():
    miRs = [] ## List to hold miRNAs
    #fh = open(args.miRNAFile)
    
    miRNA_file_clean = CleanHeader(args.miRNAFile)## Not required when info is downloaded from server in case of miRPage script
    fh_miRNA = open(miRNA_file_clean, 'r')
    fh_out2 = open('miRinput_RevComp.fa', 'w')
    mir_base = fh_miRNA.read()
    mir_blocks= mir_base.split('>')
    for i in mir_blocks[1:]:
        #print (i)
        block = i.strip('\n')##Remove newline from end of entry
        ent = block.split('\n')##Use the newline between header and sequence to split
        #print (ent)
        #print ('%s,%s,%s,%s' % (ent[0],'None','None',ent[1]))
        miRs.append((ent[0],'None','None',ent[1]))## None at position 2 and 3 because so that name and seq at same index when data collected from table
        fh_out2.write('>%s\n%s\n' % (ent[0],ent[1].translate(str.maketrans("AUTGC","TAACG"))[::-1]))## make rev comp of miRNA so that it matches the target site in genome rather than mapping miRNA to genome - in target finder file make sure that miRNA is complemented again to main original seq but not direction
    fh_miRNA.close()
    mirTable = 'None'##To be filled in final table
    print ('Total number of miRNAs in given file: %s\n' % (len(miRs)))
    
    fh_out2.close()
        
    #for i in miRs:
    #    print (i)
        
    return miRs ##miRs holds the list of miRNA name and query where as miRtable holds flag -table name or local

#Attempt3 - Using Bowtie 2 - Goal is to use the same index made for Degradome mapping
def tarFind3(frag):

    file_out = './predicted/%s.targ' % (frag.rpartition('.')[0]) ## Result File

    ### Make or select index
    index = "./index/%s_index" % (frag) ##General index name - Many index files will me made with different extensions
    if args.indexStep:
        print('**Creating index of cDNA/genomic sequences:%s\n**' % (index))
        retcode = subprocess.call(["bowtie2-build", frag, index])

    else: ### Check for existing index
        if os.path.isfile('%s.1.bt2' % index): ## The index file have many parts - so check forced for only one type of extension
            retcode = 0
            print('**Found index of cDNA/genomic sequences:%s\n**' % (index))
        else:
            print('**Could not find index of cDNA/genomic sequences:%s\n**' % (index))
            sys.exit()

    if retcode == 0: ### Index creation sucessful or index already exists
        print ('Predicting targets for frag:%s using index:%s' % (frag,index))
        nthread2 = str(nthread)
        if args.tarPred == 'H': ## Heurustic
            intervalFunc = str("L,4,0.1")
            minScoreFunc = str("G,-20,-2")
            readGap = str("24,8")
            refGap = str("12,8")
            ### Changed -D 5 to 6, changed -R 1 to 2 | Jan 13 -D 6 -> -D 3
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 3","-R 2","-N 1","-L 8","-i","S,4,0.5","--rdg","24,8","--rfg","12,8","--min-score","G,-20,-2","--norc","--no-unal","--no-hd","-p",nthread2, "-x", index, "-f" ,"miRinput_RevComp.fa","-S", file_out])
        elif args.tarPred == 'E': ##Exhaustive
            print ("You chose 'Exhaustive mode' for target identification - Please be patient")
            intervalFunc = str("L,2,0.1")
            minScoreFunc = str("G,-20,-2")
            readGap = str("24,8")
            refGap = str("12,8")
            #### Jan 13 -D 7 -> -D 4 |
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 4","-R 2","-N 1","-L 6","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nthread2, "-f", index, "miRinput_RevComp.fa","-S", file_out])

        else:
            print ('''\nPlease choose correct target prediction mode - Heuristic (H) or Exhaustive (E)\n
                   System will exit now''')
            sys.exit()

    ### Check for proper completion of Target prediction
    if retcode2 == 0:## The bowtie mapping exit with status 0, all is well
                    print('\n miRNAs mapped to Fragment: %s' % (frag))
    else:
        print ("There is some problem with miRNA mapping '%s' to cDNA/genomic seq index" % (frag))
        print ("Script exiting.......")
        sys.exit()

###PArse the bowtie based mapping and generate score
def tarParse3(targComb):
    
    print ('\n**Target prediction results are being generated**')
    ## Input / Output file ######
    print ("File for parsing: '%s' in predicted folder\n" % (targComb))
    fh_in = open(targComb,'r')
    TarPred =  './predicted/%s.parsed.csv' % (targComb.rpartition('/')[-1]) ### Similar to parsed target finder format
    fh_out = open(TarPred,'w')
    fh_out.write('miRname,Target,BindSite,miRseq,tarSeq,Score,Mismatch,CIGAR\n')
    
    #### Regenerate Target sequence with all features #####
    acount = 0 ##Total number of interactions from predictions
    parseCount = 0 ## Total number of interactions scores and written to result file
    for i in fh_in:
        #print(i)
        acount += 1
        ent = i.strip('\n').split('\t')
        #print('\n%s\n' % ent)
        miRrevcomp = ent[9] ### miRNA complemented and reversed to map genome using bowtie. That is target sequence if mimatches and gaps are added
        miRrev = miRrevcomp.translate(str.maketrans("TACG","AUGC")) ## Re-translated to get miR but still in reverse orientation - OK      
        tarHash = list(miRrevcomp) ## Strings are immutable covert to list - To rebuilt a traget seq
        #tar = miRrev
        #print('Original read mapped i.e miRNA revcomp',miRrevcomp)
        
        ##gap/bulges - Identify gaps/bulges and modify miRNA read used for mapping to regenerate target -  Add gap to target seq first to make miR length comparable to target
        gapinfo = ent[5]
        gappos = re.split("[A-Z]",gapinfo) ## In python format - gap in target seq and bulge in miRNAseq
        gapNuc = re.findall("[A-Z]",gapinfo)
        posCount = 0
        for x,y in zip(gappos[:-1],gapNuc):## In gap pos list which is made of alphabet splitting there is always am empty value at end because string has alphabet at last
            #print(x,y)
            if y == 'I':
                #tarHash.insert(posCount,'-') ## Another method as below need to time which is fast
                tarHash[posCount] = '-' ###OK
                posCount += int(x)
            else:
                posCount += int(x)       
        #print('Target seq after gap manipulation: %s' % (''.join(tarHash)))
        #print('This is the mirna in complement',miRrev)
        
        ##Mismatches - Identify mismatches and modify miRNA read used for mapping to regenerate target
        misinfo = ent[-2].split(':')[-1] ## Reverse index because XS:i is optional column ## MD:Z:16C3 - these positions are from references - so if there is an insertion/bulge in miRNA i.e. gap that it should be added to these positions before editing miRNA to tar
        #print ('This is the mismatch info block:%s' % (misinfo))
        mispos = re.split("[A,T,G,C,N]",misinfo) ##Found N in one case so included, N confimed in sequence too, will be counted as mismatch
        misposCorrect = [int(x)+1 for x in mispos] ## add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
        misNuc = re.findall("[A,T,G,C,N]",misinfo) ## Found N in one case so included, N confimed in sequence too, will be counted as mismatch
        posCount = 0
        for x,y in zip(misposCorrect,misNuc):
            #print(x,y)
            posCount += x ## Covert bowtie pos to python format
            gaps = tarHash[:posCount-1].count('-') ## Convert bowtie pos to python format -  Can give problem if more than one gap - but more than one gap not allowed V07 modification
            #print ('Position of mismatch:%s' % (posCount))
            tarHash[posCount-1+gaps] = y

        tar = ''.join(tarHash).replace("T","U") ### target converted to RNA format, will help in catching wobbles ahead
        bindsite = '%s-%d' % (ent[3],int(ent[3])+(len(miRrev)-1))

        
        ### Calculate score #######
        gap = [] ## List  to hold gap pos
        mis = [] ## List to hold mismatch position
        wobble = [] ## List to hold Wobble pos
        nt_cnt = 1 ## Keep track of actual position,
        
    #print('miRNA: %s\n%s' % (miRrevcomp[::-1],miRrevcomp[::-1].replace("T","U") ))

        #for x,y in zip(miRrevcomp[::-1].replace("T","U"),tar[::-1]):## Orientation changed to read from 5' miRNA
        for x,y in zip(miRrevcomp[::-1].replace("T","U"),tar[::-1]):## Orientation changed to read from 5' miRNA
            #print(miRrev[::-1][nt_cnt-1],x,y)## Print miRNA, rev complemmnetry miRNA used for matching, target
            if x == '-' or y == '-':
                #print('gap')
                gap.append(nt_cnt)
                if y == '-':
                    nt_cnt+=1
                
            elif x == 'A' and y == 'G': ### If in reference its 'G' than miRNA should have 'U' i.e. T but this is revcomplememnt of miRNA so complement of 'U' is A - Tested OK - v08 modifcation
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == 'C' and y == 'U': ### If in reference its 'U' than miRNA should have 'G' but this is rev complememnt of miRNA so complement of 'G' is C - Tested OK - v08 modification
                #print ('wobble')
                wobble.append(nt_cnt)
                nt_cnt+=1
            elif x == y:
                #print('match')
                nt_cnt+=1
            else:
                #print('mismatch')
                mis.append(nt_cnt)
                nt_cnt+=1
                
        #print('MimatchList:%s | GapList = %s | WobbleList = %s' % (mis,gap,wobble)) ## Poistion of mismatch gap and wobble

        score = 0   ## Initialize
        #print (mis)
        
        if args.tarScore == 'S': ## Allowed 3 MM, 2 Wob, 1 Gap
            mis2 = list(mis)
            #if set([10,11]).issubset(mis): ## Works well but took 1 sec more than below in Rice timed test
            if 10 in mis and 11 in mis: ## Check for sunsequent mismatch at 10 and 11 if yes than strict penalty ## if set(['a','b']).issubset( ['b','a','foo','bar'] )
                score += 2.5
                #print('Removing 10')
                mis2.remove(10)
                #print ('Removing 11')
                mis2.remove(11) ## So that they are not counted again
                
            for i in mis2:
                    score += 1
            for i in gap:
                score += 1.5
            for i in wobble:
                if (i+1 in mis) or (i-1 in mis): ## Mismatches around wobble - Strong penalty
                    score += 1.5
                elif (i+1) in mis and (i-1 in mis): ## Mismatches on both sides - Stronger penalty
                    score += 2
                else:
                    score += 0.5
        else:
            ##Heuristic and Exhaustive
            for i in mis:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in gap:
                if i>= 2 and i<=13:
                    score += 2
                else:
                    score += 1
            for i in wobble:
                if i>= 2 and i<=13:
                    score += 1
                    #print ('Wobble pos:%s' % (i))
                else:
                    score += 0.5
                    #print ('Wobble pos:%s' % (i))
        ###################
            
        #print(ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo)## MiRname, Tarname, mirSeq,Taerseq,binding site
        fh_out.write('>%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo))
        parseCount  += 1
    
    
    print("Total number of interactions from 'miRferno':%s AND total interactions scored: %s" % (acount,parseCount))
    fh_in.close()
    fh_out.close()
    
    


    return TarPred

### Convert tag count files to fasta file without multiplying tag count to read - So every read is uniq
def tag2FASTA2(lib):
    print("'%s' tag count file being converted to FASTA format" % (lib))
    fh_in = open(lib,'r') ### Open tag count file
    fh_out = open('./PARE/%s_PARE_tags.fa' % (lib), 'w')##Naming file with lib_ids name, 0000 means all the libs
    tag_num = 1 ### for tag naming that will be used in FASTA format
    for tag in fh_in:##All the entries of the library
        #print(tag.strip('\n').split('\t'))
        ent = tag.strip('\n').split('\t')
        tag = ent[0]
        if len(tag) >= args.tagLen: ##Length of tag specified by the user
            fh_out.write('>%s\n%s\n' % (tag_num, tag[:args.tagLen]))
            tag_num += 1
        else:
            #print ('Length is not 20nt')
            pass
    fh_out.close()

##Map degradome to the transcripts
def mapdd2trans(anIndex):##Creates index on fly and map PARE tags to index  
    mismatch = str(0) ### mismatch value for bowtie mapping
    nthread2 = str(nthread)
    index = anIndex.rsplit('.', 2)[0]
    indexLoc = './index/%s' % index
    #for lib in libs:
    dd_file = ('./PARE/%s_PARE_tags.fa' % (templib))
    map_out = ('./dd_map/%s_%s_map' % (templib,index))
    print ('\n**The library %s is being mapped to transcriptome index file: %s**\n' % (dd_file,indexLoc))
    ### --very--fast settings does not effect number of reads mapped as we are not allowing any mismatches | --no
    #retcode2 = subprocess.call(["bowtie2", "-a","--norc","--no-hd", "-t","-p",nthread2, "-f", indexLoc, dd_file,"-S",map_out]) ### Default mode is end-end
    
    retcode2 = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread2, "-f", indexLoc, dd_file,"-S",map_out]) ###Score min sets limimum score to '0'. In ened to end mode that means no mismatch/gap
    ## Optimized | no-unal was the culprit | 
    #retcode2 = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "-t","-p",nthread2, "-f", indexLoc, dd_file,"-S",map_out]) ###Score min sets limimum score to '0'. In ened to end mode that means no mismatch/gap
    
    if retcode2 == 0:## The bowtie mapping exit with status 0, all is well
        print('\nDegradome from PARE lib: %s mapped to cDNA/Trascript file' % (templib)) 
    else:
        print ("There is some problem with mapping of PARE lib: %s to cDNA/genomic seq index" % (templib))
        print ("Script exiting.......")
        sys.exit()
        
##Combine Files given as list - Make it an independent function
def FileCombine():

    print('\n****************************************')
    #targ_fls = [file for file in os.listdir() if re.search(r'.*.[0-9].*.targ', file)] - Does not recognize unslpit fasta file
    #targ_fls = [file for file in os.listdir() if re.search(r'.*.targ', file)]
    targ_fls = [file for file in os.listdir('./predicted') if file.endswith ('.targ')]
    print ('Target files:',targ_fls)
    print ('\nCombining all the target prediction files for parsing and scoring\n')
    
    targComb = './predicted/All.targs'
    targ_out = open(targComb ,'w')
    
    for x in targ_fls:
        print (x)
        targfile = open('./predicted/%s' % (x), 'r')
        #targfile.readline()
        data = targfile.read()
        targfile.close()
        targ_out.write(data)
    
    targ_out.close()
        
    return targComb

##Module to clean headers of cDNA/transcript and/or miRNA input file
def CleanHeader(filename):
    #read file
    fh_in=open(filename, 'r')
    #write file
    out_file = ('%s_new_head.fa' % (filename))
    fh_out =open(out_file, 'w')
    
    print ('\nProcessing "%s" file to clean FASTA headers\n' % (filename))
    
    acount = 0 ## count the number of entries
    for i in fh_in:
        if re.match('>', i):
            header = i.split()##change here if separater is not white space
            new_head = header[0].split('|')[0]## Change here for number of fields you want strating from 0
            fh_out.write('%s\n' % new_head)
            acount+=1
    #        print(i)
    #        print(new_head)
        else:
            fh_out.write('%s' % i)
        
    fh_in.close()
    fh_out.close()
    return out_file

    print('The fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file, acount))

def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = round((args.accel/int(nthread))+1) ## 1 added so as to avoid 0 processor being allocated in serial mode
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)
    
def PPmultiple(module,alist1,alist2):
    start = time.time()
    npool = Pool(int(args.accel))
    npool.map(lambda args: module(*args), alist2)

def PPResults(module,alist):
    npool = Pool(int(args.accel))    
    res = npool.map_async(module, alist)
    results = (res.get())
    npool.close()
    return results
    
def feed(queue, parlist):
    print ('Feeder function started')
    for par in parlist:
        print ('Echo from Feeder: %s' % (par))
        queue.put(par)
    print ('**Feeder finished queing**')

def calc(queueIn, queueOut):
    print ('Worker function started')
    while True:
        try:
            par = queueIn.get(block = False)
            print ('Echo from Worker \n Dealing with:', par)
            res = function(par)
            queueOut.put((par,res))
        except:
            break
    print ('**Worker finished **')

def write(queue, fname):
    print ('Writer function started')
    fhandle = open(fname, "w")
    while True:
        
        try:
            par, res = queue.get(block = False)
            print >>fhandle, par, res
        except:
            break
    fhandle.close()

def readFile(filename):
    """Read a file into memory

    Args:
        filename: Name of the file to be read
    Returns:
        Variable wholeFile containing the entire file

    """

    f = open(filename, 'r')
    wholeFile = f.readlines()
    f.close()

    # Strip any newline characters from the input files
    for i in range(len(wholeFile)) :
        wholeFile[i] = wholeFile[i].replace('\r', '')
        wholeFile[i] = wholeFile[i].replace('\n', '')
    
    return(wholeFile)

def createDensityDataStructure(PAGeFile):
    """Create data structure for PAGe input

    Args:
        filename: wholeFile as output from readFile function
    Returns:
        Dictionary of the PAGe file. Format is the target gene as the key
        with a two dimension array as its items. The array will contain the
        position, raw read counts, HNA, # unique reads and catergory of PARE
        reads, respectively.
        List of the category proportions

    """

    PAGeDict = {}
    # List containing the category proportions
    categoryList = []

    for i in range(len(PAGeFile)):
        # Target line should lead with a '>' character. If this is seen,
        # a new target has been found and so a new key should be entered
        # into the dictionary
        if(PAGeFile[i][0] == '>'):
            targetGene = PAGeFile[i][1:]
            PAGeDict[targetGene] = {}

        # Ignore the following size lines leading with # unless the word
        # Category is seen after
        elif(PAGeFile[i][0] == '#'):
            try:
                # If the read line has the format containing the various
                # probability proportions, append the proportions to the 
                # categoryArray
                if(PAGeFile[i].split('_')[1].split('=')[0] == 'fraction'):
                    categoryList.append(float(entry.split('_')[1].
                        split('=')[1]))
            except:
                pass
        else:
            # Create a dictionary of dictionaries with the positions
            # being secondary keys
            position = PAGeFile[i].split('\t')[0]
            item = PAGeFile[i].split('\t')[1:]
            PAGeDict[targetGene][position] = item
#            # Create a dictionary of lists containing each position
#            toAppend = []
#            for j in range(len(PAGeFile[i].split('\t'))):
#                toAppend.append(PAGeFile[i].split('\t')[j])
#            PAGeDict[targetGene].append(toAppend)

    return(PAGeDict, categoryList)

def createTargetFinderDataStructure(targetFinderFile):
    """Create data structure for targetFinder input

    Args:
        targetFinderFile: wholeFile as output from readFile function
    Returns:
        2D Array of the entire targetFinder input file

    """

    targetFinderList = []

    for i in range(len(targetFinderFile)):
        targetFinderList.append(targetFinderFile[i].split(','))

    # Sort the targetFinderList on miRNA name and targetFinder score
    targetFinderListSorted = sorted(targetFinderList, 
        key=operator.itemgetter(0,5))

    return targetFinderListSorted

def pValueCalculator(target, targetFinderList, proportion):
    """ Function to calculate p-value

    Args:
        proportion: Proportion for given category

    Returns:
        p-value

    """

    r = robjects.r
    miRNAName = target[0]
    score = target[5]

    # If score is X.5, include X and X.0 values in the n for pValue calc
    if(float(score) % 1):
        n = sum(x.count(miRNAName) for x in targetFinderList if x[5].split('.')[0] == score.split('.')[0]) ## RH Change to allow X.0 and X.5 as equal.
    else:
        n = sum(x.count(miRNAName) for x in targetFinderList if float(x[5]) ==
            float(score))
    pval = 1-(r.pbinom(0,n,proportion))[0]
    return pval

def validatedTargetsFinder(PAGeDict):
    """Perform the mapping. Take all entries from targetFinderList and
       identify if a target location matches to the 10th or 11th position
       of the miRNA.

    Args:
        PAGeDict: dictionary of PAGe file
        targetFinderList: list of target finder file
        categoryList: list of category proportions
    Returns:
        validated targets in form of targetFinder file with appended
        cleavage site and p-value.

    """ 
    validatedTargets = []
    for target in targetFinderList:
        gene = target[1]
        # Get the start location of the target
        end = int(target[2].split('-')[1])

        if(gene in PAGeDict.keys()):
            currDict = PAGeDict[gene]
            cleavageSite = []
            location = 0
            targetAbundances = []
            targetCategories = []
            targetLocations = []
            # Find if there is a PARE read at the 10th, 11th or 12th position
            # If there is, append the read abundance, the category and
            # cleavage location to their respective list
            if(str(end-9) in currDict):
                targetAbundances.append(currDict[str(end-9)][0])
                targetCategories.append(currDict[str(end-9)][1])
                targetLocations.append(end - 9)
            if(str(end-10) in currDict):
                targetAbundances.append(currDict[str(end-10)][0])
                targetCategories.append(currDict[str(end-10)][1])
                targetLocations.append(end - 10)
            if(str(end-11) in currDict):
                targetAbundances.append(currDict[str(end-11)][0])
                targetCategories.append(currDict[str(end-11)][1]) 
                targetLocations.append(end - 11)

            # If there is a PARE cleavage at any of the above positions,
            # find the best candidate for retainer.
            if(targetCategories):
                ## Debugging statement retained for quick analysis. Use with
                ## Debugging code below
                #print(targetAbundances, targetCategories, targetLocations)

                # If there is only one minimum category, use this target as
                # most probably cleavage location
                if(targetCategories.count(min(targetCategories)) == 1):
                    cleavageIndex = targetCategories.index(min(
                        targetCategories))
                    location = targetLocations[cleavageIndex]
                    cleavageSite = currDict[str(location)]

                # If there is more than one minimum category, we must filter
                # further to base our target of interest on greatest read
                elif(targetCategories.count(min(targetCategories)) > 1):
                    # Get all indices of the minimum category
                    cleavageIndices = [i for i, x in enumerate(
                        targetCategories) if x == min(targetCategories)]

                    # Get list of abundances with minimum categories
                    abundances = [targetAbundances[index] for index in 
                        cleavageIndices]
                    
                    # Cleavage index will be the read with the greatest
                    # abundance. If there is a still a tie, use the 
                    # index of the category with greatest read abundance.
                    # If there is a tie, the lowest index will always be used.
                    cleavageIndex = targetAbundances.index(max(abundances))
                    location = targetLocations[cleavageIndex]
                    cleavageSite = currDict[str(location)]
            
            if(location):
                ## Debugging statement retained in conjunction with above.
                ## Shows if cleavage is 10th, 11th or 12 position. (not
                ## coordinated with locations that output prior.)
                #print(cleavageSite, -(location-end) + 1)
                windowSum = 0
                toAppend = list(target)
                # The category score is the 2nd position.
                categoryScore = cleavageSite[1]
                # Calulate the p-value.
                pValue = pValueCalculator(target, targetFinderList, 
                    categoryList[int(categoryScore)])

                for i in range(location-5, location+6):
                    if(str(i) in currDict):
                        windowSum += int(currDict[str(i)][0])

                # Add one to the location because we need to account for the
                # 0th position of the index     
                toAppend.append(str(location)) ### Editing for biological understanding - original (str(location+1))
                # Add PARE abundance
                toAppend.append(cleavageSite[0])
                # Add sum of reads within 5 bp of cleavage site in each
                # direction.
                toAppend.append(str(windowSum))
                # Add ratio of abundance of cleavage site to sum within 5 bp 
                # of the cleavage site in each direction
                toAppend.append(str(float(int(cleavageSite[0])/windowSum)))
                # Add category at cleavage site
                toAppend.append(str(categoryScore))
                # Append the p-value to the toAppend list
                toAppend.append(str("%.6f" % pValue))
                validatedTargets.append(toAppend)
    
    return(validatedTargets)

def createPAGeIndex(bowtieFilename):
    """Create data structure for targetFinder input. Then, performs the
       mapping. Take all entries from tagCountFile and identify all tags
       that map to 
    Args:
         bowtieFilename: name of map file
    Returns:
        Dictionary with the proper information for dd_PAGe output. Key
        is gene name and the value is a list with tuples as values. The first
        element in the tuple is the location of the hit, and the second
        element is the abundance for that sequence.
        A list of all abundance values are also returned
    """

    # Initialize dictionary for bowtie entries
    bowtieDict = {}
    
    # Append the directory to the bowtieFilename
    bowtieFilename = 'dd_map/' + bowtieFilename
    bowtieFile = readFile(bowtieFilename)

    if args.repeats:
        # Loop through the entire bowtie files to load these into the dictionary
        for entry in bowtieFile:
            if(entry.split('\t')[4] != '255'):
                pass
            else:
                # Store gene, location and sequence from each entry into variables
                gene = entry.split('\t')[2]
                location = entry.split('\t')[3]
                sequence = entry.split('\t')[9]

                # Try to identify if the sequence is a key in the bowtie dictionary
                try:
                    bowtieDict[sequence]
                # If not, create one with a value of an empty dictionary
                except:
                    bowtieDict[sequence] = {}

                # Try to identify if gene is a key in the nested dictionary
                try:
                    bowtieDict[sequence][gene]

                # If not, create one with a value of an empty list
                except:
                    bowtieDict[sequence][gene] = []

                # Append the location of the sequence in the gene to the list
                bowtieDict[sequence][gene].append(location)        

    else:
        for entry in bowtieFile:
            # Store gene, location and sequence from each entry into variables
            gene = entry.split('\t')[2]
            location = entry.split('\t')[3]
            sequence = entry.split('\t')[9]

            # Try to identify if the sequence is a key in the bowtie dictionary
            try:
                bowtieDict[sequence]
            # If not, create one with a value of an empty dictionary
            except:
                bowtieDict[sequence] = {}

            # Try to identify if gene is a key in the nested dictionary
            try:
                bowtieDict[sequence][gene]

            # If not, create one with a value of an empty list
            except:
                bowtieDict[sequence][gene] = []

            # Append the location of the sequence in the gene to the list
            bowtieDict[sequence][gene].append(location)

    PAGeDict = {}
    allHits = []

    # Loop through all tags and their hit counts
    for entry in tagCountFile:
        # sequence and hits used multiple times, so store as variables.
        sequence = entry.split('\t')[0]
        hits = entry.split('\t')[1]

        # If the tag exists in the bowtie file, assign the hits to the 
        # data structure. Otherwise, a pass must be made
        try:
            # Loop through all gene locations
            for key in bowtieDict[sequence].keys():
                # Test to see if key exists in PAGeDict yet. If not, create
                # it
                try:
                    PAGeDict[key]
                except:
                    PAGeDict[key] = {}

                # Loop through all locations of the gene in which the tag maps
                # to
                for location in bowtieDict[sequence][key]:
                    # Create a dictionary off of the gene with hits as its key
                    PAGeDict[key][location] = int(hits)
                   
                    # Append the hits to the hits if hits > 2
                    if(int(hits) > 2):
                        allHits.append(int(hits))

        except:
            pass

    return(PAGeDict, allHits)

def unambiguousBaseCounter(transcriptomeFile, tagLen):
    """Get the counts of ambiguous bases in the transcriptome file as well
       as counts of ambiguous bases that are within the ends of the
       transcriptome +/- the tagLen.

    Args:
        transcriptomeFile: wholeFile of the trnascritome fasta file
        tagLen: Number of bases that an N is allowed to be away from the ends
            of the gene in order to be counted

    Returns:
        Total number of ambiguous bases and ambiguous bases tagLen-bp away from
        the ends of the gene.

    """

    baseCounts, baseCountsOffTagLen = 0, 0

    # Start at line 1, reading only sequences and jump 2 lines to read
    # the next sequence until end of file
    for i in range(1,len(transcriptomeFile),2):
        currentLine = transcriptomeFile[i]
        baseCounts += len(currentLine) - currentLine.count('N')
        baseCountsOffTagLen += (len(currentLine) - 2 * tagLen) - currentLine[
            tagLen:len(currentLine)-tagLen].count('N')

    return baseCounts, baseCountsOffTagLen

def writeDensityFile(PAGeDict, mode, allHits, baseCounts, baseCountsOffTagLen,
    outputFile, transcriptomeFilename, library):
    """Write validated targets to an output file

    Args:
        PAGeDict: Dictionary of genes with target locations and hits for
            each tag mapped to those locations.
        mode: 0 (genic) or 1 (intergenic)
        allHits: List of all abundance values
        baseCounts: Total number of unambiguous bases
        baseCountsOffTagLen: Total number of unambiguous bases tagLen-bp away from
            the ends of the gene.
        outputFile: file to output PAGe inforation
        transcriptomeFilename: Name of transcriptome file
        library: Name of library being analyzed

    """
    # Create data structure to keep track of the counts of each category
    categoryCounts = [0,0,0,0,0]
    categoryList = []
    f = open(outputFile,'w')
    #g = open('output/%s_categoryCounts.csv' % library,'w')
    # Get the count of the total genes
    numGenes = len(PAGeDict.keys())

    if(mode == 1):
        globalMedian = numpy.median(allHits)
        seventyFivePercentile = stats.scoreatpercentile(allHits, 75)
        ninetyPercentile = stats.scoreatpercentile(allHits, 90)
        print('median = %s\nseventyFivePercentile = %s\nninetyPercentile = %s'
            % (globalMedian, seventyFivePercentile, ninetyPercentile))
        #g.write('gene,cat0,cat1,cat2,cat3\n')

    else:
        #g.write('gene,cat0/1,cat2,cat3\n')
        pass
        

    # Sort genes so that genes are in alphabetical order
    for gene in sorted(PAGeDict.keys()):
        catTwoAbun = []
        catThreeAbun = []
        # Genic tracks category 0 and 1 as one
        if(mode == 0):
            catZeroOneAbun = []
        # Intergenic separates categories 0 and 1
        else:
            catZeroAbun = []
            catOneAbun = []
       
        f.write('>%s\n' % str(gene))
        #g.write(str(gene) + ',')
        if(mode == 0):
            geneHits = []
            multHitsFlag = 0

            # Store all of the abundance values for each location
            for hits in PAGeDict[gene].values():
                geneHits.append(hits)

            # Calculate median and max on gene
            median = numpy.median(geneHits)
            maxHit = max(geneHits)
            if(len([i for i, x in enumerate(geneHits) if x == maxHit]) > 1):
                multHitsFlag = 1

            # Sort all locations in which a tag maps in that gene so they can
            # be listed in order. Must use key as int because in dictionary,
            # locations are stored as strings
            for location in sorted(PAGeDict[gene].keys(),key=int):
                hits = PAGeDict[gene][location]

                # Calculate category
                if(hits == 1):
                    category = '4'
                    categoryCounts[4] += 1
                    PAGeDict[gene][location] = (hits, 4)
                elif(hits <= median):
                    category = '3'
                    categoryCounts[3] += 1
                    PAGeDict[gene][location] = (hits, 3)
                    catThreeAbun.append(hits)
                elif(hits > median and hits != maxHit):
                    category = '2'
                    categoryCounts[2] += 1
                    PAGeDict[gene][location] = (hits, 2)
                    catTwoAbun.append(hits)
                elif(hits > median and multHitsFlag):
                    category = '1'
                    categoryCounts[1] += 1
                    PAGeDict[gene][location] = (hits, 1)
                    catZeroOneAbun.append(hits)
                else:
                    category = '0'
                    categoryCounts[0] += 1
                    catZeroOneAbun.append(hits)
                    PAGeDict[gene][location] = (hits, 0)
                f.write('%s\t%s\t%s\n' % (str(location), str(hits), category))

            #g.write(str(max(catZeroOneAbun) if catZeroOneAbun else 0) + ',' +
                #str(max(catTwoAbun) if catTwoAbun else 0) + ',' +
                #str(max(catThreeAbun) if catThreeAbun else 0) + '\n')

        elif(mode == 1):
            # Sort all locations in which a tag maps in that gene so they can
            # be listed in order.
            for location in sorted(PAGeDict[gene].keys(), key=int):
                hits = PAGeDict[gene][location]

                # Calculate category
                if(hits <= 2):
                    category = '4'
                    categoryCounts[4] += 1
                    PAGeDict[gene][location] = (hits, 4)
                elif(hits <= globalMedian):
                    category = '3'
                    categoryCounts[3] += 1
                    PAGeDict[gene][location] = (hits, 3)
                    catThreeAbun.append(hits)
                elif(hits > globalMedian and hits <= seventyFivePercentile):
                    category = '2'
                    categoryCounts[2] += 1
                    PAGeDict[gene][location] = (hits, 2)
                    catTwoAbun.append(hits)
                elif(hits > seventyFivePercentile and
                        hits <= ninetyPercentile):
                    category = '1'
                    categoryCounts[1] += 1
                    PAGeDict[gene][location] = (hits, 1)
                    catOneAbun.append(hits)
                else:
                    category = '0'
                    categoryCounts[0] += 1
                    PAGeDict[gene][location] = (hits, 0)
                    catZeroAbun.append(hits)

                f.write('%s\t%s\t%s\n' % (str(location), str(hits), category))

            #g.write(str(max(catZeroAbun) if catZeroAbun else 0) + ',' + 
                #str(max(catOneAbun) if catOneAbun else 0) + ',' +
                #str(max(catTwoAbun) if catTwoAbun else 0) + ',' + 
                #str(max(catThreeAbun) if catThreeAbun else 0) + '\n')

    f.write('# Transcriptome=%s\n' % transcriptomeFilename)
    f.write('# Genes=%s\n' % numGenes)
    f.write('# Uncorrected non-ambiguous bases=%s\n' % baseCounts)
    f.write('# Eligible bases for degradome-derived 5 prime ends=%s\n' %
        baseCountsOffTagLen)
    for i in range(len(categoryCounts)):
        f.write('# Category %s_bases=%s\n' % (i, categoryCounts[i]))
    for i in range(len(categoryCounts)):
        categoryList.append(categoryCounts[i] / baseCountsOffTagLen)
        f.write('# Category %s_fraction=%s\n' % (i, categoryCounts[i] / 
            baseCountsOffTagLen))
    
    f.close()
    #g.close()
    return(categoryList)

def writeValidatedTargetsFile(header, validatedTargets, outputFile):
    """Write validated targets to an output file

    Args:
        header: the stripped header from the target file
        validatedTargets: list of the validated targets
        outputFile: file to output validated targets to

    """

    # Get indicies of various elements
    pValueIndex = len(validatedTargets[0]) - 1
    categoryIndex = len(validatedTargets[0]) - 2
    windowRatioIndex = len(validatedTargets[0]) - 3
    
    # Code to perform noise corrected pValue
    for i in range(len(validatedTargets)):
        validatedTargets[i].append("%.6f" % (float(validatedTargets[i][pValueIndex])/float(validatedTargets[i][windowRatioIndex])))

    # get corrected p-value index
    correctedPValueIndex = len(validatedTargets[0]) - 1

    f = open(outputFile,'w')

    f.write(header + ',cleavage position,PARE reads,10 nt window abundance,'\
        'PARE reads/window abundance,category,p-value,noise corrected p '\
        'value\n')

    # noiseFilter requires that p value be <= .25 and window ratio be
    # >= .25 for a tag to be considered
    if(args.noiseFilter):
        for target in validatedTargets:
            # Include any target with pvalue <= .25 and with window ratio
            # >= .25
            if((float(target[correctedPValueIndex]) <= .25) and (float(target[
                    windowRatioIndex]) >= .25)):
                # If category 4s are not to be included and the category
                # of the current category is 4, then omit it.
                if(args.cat4Filter and (int(target[categoryIndex]) == 4)):
                    pass
                else:
                    for j in range(len(target)):
                        f.write(str(target[j]))
                        # If at last element, write new line
                        if(j == correctedPValueIndex):
                            f.write('\n')
                        # Otherwise, write a comma
                        else:
                            f.write(',')

    # If noise filter is not on, then include any target with p value < .5 
    else:
        for target in validatedTargets:
            # Include any target with pvalue < .5
            if(float(target[correctedPValueIndex]) < .5):
                for j in range(len(target)):
                    f.write(str(target[j]))
                    # If at last element, write new line
                    if(j == correctedPValueIndex):
                        f.write('\n')
                    # Otherwise, write a comma
                    else:
                        f.write(',')

#### MAIN FUNCTION ################################
def main():
    if args.generateFasta:
        coords = extractFeatures(args.genomeFile,args.gffFile) ## Extracts Coords from GFF3
        fastaOut = getFASTA1(args.genomeFile,coords) ##Creates FASTA file
        print('This is the extracted file: %s' % (fastaOut))
    # If featureFile is given
    else:
        print("\nThe input FASTA file is considered 'as is' for analysis\n")
        fastaOut = args.featureFile ### Make it better
    
    ### FRAGMENTATION ###################
    ###Script Timer
    runLog = 'runtime_%s' % datetime.datetime.now().strftime("%m_%d_%H_%M")
    fh_run = open(runLog, 'w')
    print('tarPred: %s | tarScore: %s | Uniq filter: %s' % (args.tarPred,args.tarScore,args.repeats))
    fh_run.write('tarPred:%s | tarScore: %s | Uniq filter:%s\nLibs:%s\nGenomeFile:%s | GenomeFeature:%s' % (args.tarPred,args.tarScore,args.repeats, ','.join(args.libs),args.genomeFile,args.genomeFeature))
    fh_run.write ('\nLibs: %s' % (','.join(args.libs)))
    FragStart = time.time()
    
    if args.fileFrag:
        start = time.time()###time start
        fragList = fragFASTA(fastaOut)##Fasta is a list of fragmented files
        end = time.time()
        print ('fileFrag time: %s' % (round(end-start,2)))
    else:
        fragList = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)] ## Need better handiling - What is main file name is ATH_12345.fa like from phytozome
        print ('The fragments: %s' % (fragList))
        
    FragEnd = time.time()
    print ('\n\nThe script run time is %s\n\n' % (round(FragEnd-FragStart,2)))
    fh_run.write('Fragmentation time is : %s\n' % (round(FragEnd-FragStart,2)))

    #####################################
    
    ## TARGET PREDICTION ###################
    miRs = miRinput()
    
    TFStart = time.time()
    
    ## Remove results from previous run
    if args.indexStep:
        shutil.rmtree('./index', ignore_errors=True)
        os.mkdir('./index')
    if args.tarPred and args.tarScore:
        shutil.rmtree('./predicted', ignore_errors=True)
        os.mkdir('./predicted')
        print('\nFragments to be indexed and used for TP: %s' % (fragList))
    
        start = time.time()###time start
        ### Serial mode - Test and Trouble shooting purpose
        #for i in fragList:
        #    tarFind3(i)
        ## Parallel mode
        PP(tarFind3,fragList)
        end = time.time()
        print ('Target Prediction time: %s' % (round(end-start,2)))
        
        
        targComb = FileCombine()
        #targComb = 'All.targs' ## Test - open line above when real
        
        start = time.time()###time start
        predTargets = tarParse3(targComb)
        end = time.time()
    
        #print ('Target Prediction time: %s' % (round(end-start,2)))
        
    elif not args.tarPred and args.tarScore:
        targComb = FileCombine()
        #targComb = 'All.targs' ## Test - open line above when real
        
        start = time.time()###time start
        predTargets = tarParse3(targComb)
        end = time.time()
        
        #print ('Target Scoring time: %s' % (round(end-start,2)))
    
    else: ## Target prediction is OFF
        print("!!Target prediction is OFF - Files in 'predicted' folder might be old!!")
        predTargets = './predicted/All.targs.parsed.csv'
    
    ###Timer
    TFEnd = time.time()
    print ('\n\nThe script run time is %s\n\n' % (round(TFEnd-TFStart,2)))
    fh_run.write('Target prediction time is : %s\n' % (round(TFEnd-TFStart,2)))
    #########################################
    
    ## PARE PROCESS AND MAP #################
    PAREStart = time.time()
    
    if args.tag2FASTA:
        shutil.rmtree('./PARE',ignore_errors=True)
        os.mkdir('./PARE')
        PP(tag2FASTA2,args.libs)
    
    indexFls = [file for file in os.listdir('./index') if file.endswith ('index.1.bt2')]
    print ('These are index files: ',indexFls)
    
    if args.map2DD:
        shutil.rmtree('./dd_map',ignore_errors=True)
        os.mkdir('./dd_map')
        global templib
        for templib in args.libs:
            ## Serial -Test
            #for dd in indexFls:
                #print ('Lib:%s mapped to index: %s' % (templib,dd))
                #mapdd2trans(dd)
            
            ### Parallel
            PP(mapdd2trans,indexFls)
    
    ##Timer
    PAREEnd = time.time()
    print ('\n\nPARE processing run time is %s\n\n' % (round(PAREEnd-PAREStart,2)))
    fh_run.write('PARE processing and mapping time is : %s\n' % (round(PAREEnd-PAREStart,2)))
    ## INDEXER ##############################
    ###Timer AK - will be written to run log file
    PredStart = time.time()
    
    if args.validate:
        predTargets = './predicted/All.targs.parsed.csv'
        shutil.rmtree('./PAGe',ignore_errors=True)
        os.mkdir('./PAGe')
        shutil.rmtree('./output', ignore_errors=True)
        os.mkdir('./output')
        validatedTargetsFilename = './output/MPPPValidated.csv'
    
        # Read the transcriptome and target finder file into memory
        transcriptomeFile = readFile(fastaOut)
        targetFinderFile = readFile(predTargets)
        # Strip the header from the targetFinder File, but save it first
        header = targetFinderFile[0]
        del(targetFinderFile[0])
    
        # Store the targetFinder file into memory as a 2D Array.
        global targetFinderList
        targetFinderList = createTargetFinderDataStructure(
            targetFinderFile)
    
        # Analyze transcriptome file for unambiguous bases
        baseCounts, baseCountsOffTagLen = unambiguousBaseCounter(
            transcriptomeFile, args.tagLen)
    
        for tagCountFilename in args.libs:
            PAGeIndexDict = {}
            # Variable holding all hits > 2
            allHits = []
            PAGeIndexList = []
            global tagCountFile
            tagCountFile = readFile(tagCountFilename)
            library = tagCountFilename.split('_')[0]
            PAGeOutputFilename = './PAGe/%s_PAGe' % library
            validatedTargetsFilename = './output/%s_validated' % library
    
            bowtieFiles = [file for file in os.listdir('dd_map') if file.startswith('%s' % library)]
            print("Creating PAGe Index dictionary for lib %s" % library)
            PAGeStart = time.time()
            
            # If parallel version, results have to be merged together
            PAGeIndexAndHits = PPResults(createPAGeIndex, bowtieFiles)
            # Merge the PAGeIndex and allHits computed on separate cores
            for element in PAGeIndexAndHits:
                PAGeIndexDict.update(element[0])
                allHits.extend(element[1])
                PAGeIndexList.append(element[0])

            PAGeEnd = time.time()
            print("PAGe Indexing took %.2f seconds" % (PAGeEnd - PAGeStart))
            fh_run.write("PAGe Indexing took: %.3f seconds\n" % (PAGeEnd-PAGeStart))
    
            # Write the PAGe file
            print("Writing PAGeIndex file...")
            PAGeWriteStart = time.time()
            global categoryList
            categoryList = writeDensityFile(PAGeIndexDict, args.genomeFeature,
                allHits, baseCounts, baseCountsOffTagLen, PAGeOutputFilename,
                fastaOut, library)
            PAGeWriteEnd = time.time()
            print("File written. Process took %.2f seconds" % (PAGeWriteEnd - PAGeWriteStart))
            fh_run.write("PAGe index file written. Process took %.2f seconds\n" % (PAGeWriteEnd - PAGeWriteStart))
    
            # Find the validated targets
            print("Finding the validated targets")
            validatedTargetsStart = time.time()
            
            # Run validatedTargetsFinder in parallel. Takes every other element
            # of PAGeIndexAndHits array for only PAGeIndex dictionaries
            validatedTargetsList = PPResults(validatedTargetsFinder, PAGeIndexList)
            validatedTargets = []
            for validatedTarget in validatedTargetsList:
                validatedTargets.extend(validatedTarget)

            validatedTargetsEnd = time.time()
            print("All validated targets found in %.2f seconds" % (validatedTargetsEnd - validatedTargetsStart))
            fh_run.write("All validated targets found in %.2f seconds\n" % (validatedTargetsEnd - validatedTargetsStart))
            
            # Write the validated targets to a file
            # Write the validated targets to a file
            if(validatedTargets):
                print("Writing validated targets")
                writeValidatedTargetsFile(header, validatedTargets,
                    validatedTargetsFilename)
    ##Script timer
    PredEnd = time.time()
    print ('\n\nIndexing and Prediction run time is %s\n\n' % (round(PredEnd-PredStart,2)))
    fh_run.write('Indexing and Prediction run time is : %s seconds \n' % (round(PredEnd-PredStart,2)))
    fh_run.write('Script run time is : %s\n' % (round(PredEnd-FragStart,2)))
    fh_run.close()

#### RUN ##########################################

if __name__ == '__main__':
    
    if args.accel == 'Y':
        args.accel = int(multiprocessing.cpu_count()*0.85)
    else:
        args.accel = int(args.accel)

    
    start = time.time()
    main()
    end = time.time()
    print ("The complete 'sPARTA' run time is %s seconds" % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting..\n')
    sys.exit()
