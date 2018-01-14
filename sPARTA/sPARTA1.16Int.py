#!/usr/local/bin/python3
## sPARTA-Int v1.09 - miRferno reached End of support -- Atul

#### PYTHON FUNCTIONS ##############################
import sys,os,re,time,timeit,csv,glob, string
import subprocess, multiprocessing
import shutil
import itertools as it
import operator
from multiprocessing import Process, Queue, Pool
from operator import itemgetter
import pyfasta
import mysql.connector as sql
import rpy2.robjects as robjects
from scipy import stats
import bisect
import numpy
import datetime
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
RStats = importr('stats')

#########################################################################
#### USER SETTINGS ######################################################

Local               = 'N'                               ## 'N' means internal version and 'Y' means external version i.e. local files used instead of miRNA and fasta files from server
GenomeDB            = 'RICE_MSU7_genome'                ## DB to get fastaOut file
PAREdb              = 'RICE_sbsQIFA_PARE'               ## Make sure that your Library is in the DB 

gffFile             = 'Athaliana_167_gene.gff3'         ## From Phytozome - Not used if Local == N
genomeFile          = 'genomic_seq.fa'                  ## From Phytozome - Not used if Local == N

miRNA_file          = 'No.fa'                           ## If not present in directory than miRcode (below) will be sued to fetch miRNAs from miRbase table
miRCode             = 'gma'                             ## Code of species to fetch miRNAs. Not used if miRNA file is present in directory
ver                 = '21'                              ## mirbase version to use - Not used if miRNA file is present in directory
miRTable            = 'mirBASE'                         ## DO NOT TOUCH - Not used if miRNA file is present in directory

genomeFeature       = 1                                 ## 0 for gene and 1 for inter; 2 for both
getPARELib          = 'N'                               ## PARE libraries are present on Disk 'N' or need to be fetched 'Y'
libs                =['3628','3629']                    ## Name of PARE libraries - MANDATORY

#### STEPS ######################
generateFasta       = 'Y'                               ## If Local = Y - Use GFF file and splice out genes
fileFrag            = 'N'                               ## Fragment file is 'Y' else use files split earlier - Local compatible step
##

indexStep           = 'N'
tarPred             = 'N'                               ## Target prediction required or not 'Y' and 'N' - Local compatible step
tarScore            = 'N'                               ## Target scroing - Local compatible step
##
tag2FASTAstep       = 'N'                               ## Convert PARE tag count file to FASTA file
map2DDstep          = 'N'                               ## Map PARE reads to degradome
predStep            = 'N'                               ## Reza functions
revMapstep          = 'N'                               ## Reverse map to match co-ordinates with our genome
uniqRes             = 'N'                               ## Combine results to a single file and filter out redundant entries

##### END USER SETTINGS #################################################
#########################################################################

##### Advanced Options ##################################################
TPscore             = 'S'                               ## 'S' - Seedless and 'N' for normal target scoring
TPmode              = 'E'                               ## 'H' for heuristic and 'E' for exhaustive
maxHits             = 30
dataserver          = 'raichu.dbi.udel.edu' 
destserver          = 'raichu.dbi.udel.edu'
##
tagLen              = 20                                ## Minimum PARE tag length
nthread             = 6                                 ## Need automatic calculation like nProc
splitCutoff         = 40                                ## cutoff for fragmenting FASTA file
tableUpdate         = 'Y'                               ## Update to specified table
extraAnno           = 'N'                               ## Fetched annotation from new tables - phased_loci and chr_inverted

##### Misc Options - Developer Only ####################################
bowtieFilterSwitch  = 'Y'                               ## Reza Added: Filter MAPQ column for Uniq (255) and non-uniq reads
nproc               = 'Y'                               ## Used by parallel processing
noiseFilter         = 'N'                               ## Remove noise from results i.e. corr p-val >= 0.25, PARE ratio <= 0.25 AND Cat 4 validations
cat4Show            = 'N'                               ## Filter out Category 4 - Developer option
standardCleave      = 'N'                               ## Turn off miRNA size-specific cleavage search. Cleavage checked on 9,10,11, nt
##### BE VERY CAUTIOUS ##########
TableWipe           = 'N'                               ## Wipe table before update - CAN WIPE OTHERS RESULTS - DO NOT TOUCH
#######################################################################
#                                                                     #
#                                                                     #
############# FUNCTIONS ###############################################

def ConnectToDB(server, infile):

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

## Extract coordinates from GFF file -- Will update script with latest version of scripts
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
    genome_info.sort(key=lambda k:(k[0],k[1],k[3])) ### Sorting on basis of chr_id, strand and gene start for coord calulation and reduce getFASTAServer time
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
                ##As two genes are read together, the upstream intergenic region for gene2 must be calculated in same step
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
    
    if genomeFeature == 2: ## Both gene and inter
        for ent in genome_info_inter_sort:
            print(ent)
            if ent[4] == '-': ## End of chromosome
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                
            elif int(ent[4])-int(ent[3]) > 25:###If there is no intergenic regon b/w genes or too short than filter
                coords.append(ent[0:])
                coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
    
    else:
        if genomeFeature == 0:
            genomeFilter = 'gene'
        elif genomeFeature == 1:
            genomeFilter = 'inter'
        
        for ent in genome_info_inter_sort:
            if (ent[5] == genomeFilter):
                #print(ent)
                if ent[4] == '-': ## End of chromosome
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
                    
                elif int(ent[4])-int(ent[3]) > 25:### If there is no intergenic region b/w genes or too short than filter
                    coords.append(ent[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (ent[0:]))
    
    print ("Number of coords in 'coords' list: %s" % (len(coords)))
    coords_out.close()    
    fh_in.close()
    
    return coords

## This module interacts with genome DB to take gene coords, calculate intergenic co-ordinates as well as generate a FASTA file for input
def extractCoords(con,db):
    cur= con.cursor()
    
    ###Filteration of RNA types not done at this stage becasue we need all the genes to get intergenic regions
    ###after filteration few genes will be missing and counted as intergenic region
    ###########################################################TEST LINE#################################################################
    #cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))###extra % escapes the % in query
    #####################################################################################################################################
    ##If we remove miRNA here than only gene entry will be removed but its sequence will be covered into intergenic
    if genomeFeature == 1: ##Intergenic
        print("Fetching intergenic coords - If script stucks here for more then a minute, then cancel - ctrl+c and rerun ")
        cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type NOT LIKE 'mirna' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))### Extra % escapes the % in query
    elif genomeFeature == 0: ##Protein Coding
        print("Fetching gene coords  - If script stucks here for more then a minute, then cancel - ctrl+c and rerun ")
        cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))### Extra % escapes the % in query

    genome_info = cur.fetchall() ### List with information of all the genes in genome

    ##Check if list is empty
    if not genome_info:
        print ('^^^Gene Coords query returned with an empty list..Exiting^^^')
        sys.exit()
    ####(1, 'c','AT1G01020', 5928, 8737, protein_coding)
    #print (genome_info)
    
    ##Find length of chromosomes to calculate intergenics
    cur.execute('SELECT chr_id, length FROM %s.chromosome_master' % (db))
    chromo_len = cur.fetchall()
    chromo_dict = dict(chromo_len)  ### Made a dict so that chromosome numer could be searched to get chromosome length
    # print ('These are the chromosomes: %s and their length' % (chromo_dict))
    
    genome_info_inter = genome_info ###This list will also hold intergenics
    
    ####GET INTERGENIC REGIONS AND APPEND TO THE LIST WITH GENE COORDS
    alist = []###list maintained to check if first gene on chromosome and strand shows up than intergenic is just the start of gene
    #for gene1, gene2 from izip(*[iter(genome_info)]*2):
    #for i,j in pairwise (genome_info):
    #for gene1, gene2 in it.izip(genome_info[1:], genome_info):
    #for i in range(0, int(len(genome_info))-1):###maybe for i in range(0, len(genome_info -1))
    

    for i in range(0, int(len(genome_info))+1): ## May 20-14 - Modifed from above after extensive trouble shooting - Now the last entry is read and both up and Down calculated
        #print (i)
        gene1 = (genome_info[i])
        gene2 = (genome_info[i+1])
        gene_type = 'inter' ###set to intergenic by default
        #print(gene1,gene2)
        
        ##Remove/skip redundant genes with same start and end......What about Overlapping genes????
        if gene1[3] == gene2[3] and gene1[4] == gene2[4]:
            ##gene is same/overlapping consider next gene
            pass
        
        else:
            ##Calculate coordinates of intergenic regions
            ##Is this is first gene on chromosome and strand intergenic region is from position1 - Only chr_id and strand is checked. 
            if tuple(gene1[0:2]) not in alist:
                print ('\n------Caching gene coords for chromosome: %s and strand: %s------\n' % (gene1[0], gene1[1]))
                #print ('Gene1:%s\nGene2:%s' % (gene1,gene2))
                alist.append((gene1[0:2]))
                
                inter_start1 = 1
                inter_end1 = gene1[3]-1###1 nt before start of Gene1 a.k.a the first gene on chromosome in this case
                ##As two genes are read together, the upstream intergenic region gor gene2 must be calculated in same step
                
                ## If both the genes belong to same chr and strand i.e. chromosome has atleast two genes
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]: 
                    inter_start2 = gene1[4]+1##From end of first gene of chromosome
                    inter_end2 = gene2[3]-1###Till start of second gene
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene2[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])

                    #if inter_name1 == inter_name2:
                    #print ('\nLoop1 - First Gene and more than one gene on this chr and strand')
                    #print (gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type)
                    #print (gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type)
                    
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                    genome_info_inter.append((gene2[0],gene2[1],inter_name2,inter_start2,inter_end2,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                
                ## If the first two genes are not from same strand i.e. there is only one gene on this chromosme and strand- code added after Aspragus scaffolds had just one gene
                else: ## intergenic from end of chromosme/scaffold
                    inter_start2 = gene1[4]+1##From end of first gene of chromosome
                    inter_end2 = chromo_dict[gene1[0]]###Till end of chromosome
                    
                    if gene1[1] == 'w': ##The gene is on positive strand so upstream
                        inter_name1 = ('%s_up' % (gene1[2]))
                        inter_name2 = ('%s_down' % gene1[2])
                    else: ##Its on negative strand
                        inter_name1 = ('%s_down' % (gene1[2]))
                        inter_name2 = ('%s_up' % gene1[2])
                    #if inter_name1 == inter_name2:
                    
                    #print ('\nLoop2 - First gene on this chromosme and strand but also the only one')
                    #print (gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type)
                    #print (gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type)
                    
                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name1,inter_start1,inter_end1,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene type
                    genome_info_inter.append((gene1[0],gene1[1],inter_name2,inter_start2,inter_end2,gene_type))##Chr_id_strand,intergeneic name, inter start, inter end, gene typ
                
            else:
                if gene1[0] == gene2[0] and gene1[1] == gene2[1]:###If chr_id and strands are equal than find intergenic. These are gene on same chromosme and strand
                    inter_start = gene1[4]+1###End of Gene 1
                    inter_end = gene2[3]-1 ###1 nt before start of gene 2
                    if gene2[1] == 'w': ##Positive strand
                        inter_name = ('%s_up' % (gene2[2]))
                    else:## reverse strand
                        inter_name = ('%s_up' % (gene1[2]))
                    #print ('\nLoop3 - Not the first gene on chr and strand')
                    #print (gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type)
                    genome_info_inter.append((gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type))
                
                else: ###That means gene1 is at end of one chromosome and gene 2 is begining of chromosome so we have to extract intergenic at end of one chromosome
                    inter_start = gene1[4]+1###End of gene1
                    inter_end = chromo_dict[gene1[0]]###End of chromosome searched using chromosome id of gene1 from chromosome dictionary
                    if gene1[1] == 'w':##Positive strand end of chromosme
                        inter_name = ('%s_down' % (gene1[2]))
                    else: ##Negative strand first intergenic of chromosme
                        inter_name = ('%s_up' % (gene1[2]))
                        
                    #print ('\nLoop4 - Not the first gene on chromosome and Strand AND the last gene on chromosome')
                    #print (gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)
                    #
                    genome_info_inter.append((gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type)) ## Chr_id, strand


    ## Additional check for scaffolded genomes, if there are no genes in a scffold it's whole seqeunce will be fetched as intergenic
    if genomeFeature == 1:
        for i in chromo_dict.keys():
            alen = chromo_dict[i]
            # print("Chr:%s | Length:%s" % (i,alen))
            if tuple((i,'c')) in alist:
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_c_all' % (i))
                genome_info_inter.append((i,'c',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

            if tuple((i,'w')) in alist:
                pass
            else:
                # print("Get the tuple")
                inter_name = ('chr%s_w_all' % (i))
                genome_info_inter.append((i,'w',inter_name,1,alen,'inter')) ## Chr_id, strand, name, start, stop, length

    
    
    ###Sort the list after adding intergenic regions on on basis of chr_id and strand that is essential while caching chromosme during slicing sequences
    genome_info_inter_sort = sorted(genome_info_inter, key=operator.itemgetter(0,1))
    #print(genome_info_inter_sort)
    ###Write all cooords for troubleshooting
    all_coords_out = open('all_coords', 'w')
    for i in genome_info_inter_sort:
        all_coords_out.write('%s,%s,%s,%s,%s,%s\n' % (i[0:]))
    all_coords_out.close()
            
    
    ###Filter list to remove unwanted types like miRNA,tRNA,rRNA,snoRNA,snRNA, short or no intergenic
    gene_coords_file = './gene_coords'####To check wheter coords are printed in chr_id and strand sorted or not
    coords_out = open(gene_coords_file, 'w')
    gene_coords = []## List that will hold genes to fetch, this removes unecessary RNAs and also fix miRNA double entry i.e as gene and miRNA
    for entry in genome_info_inter_sort:
        #print(entry[5])
        ###If RNA type is removed here than that region is not included in analysis but if RNA is removed in mySQL query than only gene is removed and region becomes intergenic
        #if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA'): ##Replace this IF with tuple check i.e miRNA in tuple
        if genomeFeature == 1: ## Inter
            if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA' or entry [5] == 'protein-coding' or entry [5] == 'protein_coding' or entry [5] == 'misc_RNA'):
                pass
            else:
                if entry[4]-entry[3] > 25:###If there is no intergenic regon b/w genes or too short than filter
                    #gene_coords.append(entry[:5])
                    #coords_out.write('%s,%s,%s,%s,%s\n' % (entry[0:5]))
                    gene_coords.append(entry[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (entry[0:]))
        else: ## Protein coding
            if (entry[5] == 'miRNA' or entry[5] == 'tRNA' or entry[5] == 'rRNA' or entry [5] == 'snoRNA' or entry [5] == 'inter'):
                pass
            else:
                if entry[4]-entry[3] > 25:###If there is no intergenic regon b/w genes or too short than filter
                    gene_coords.append(entry[0:])
                    coords_out.write('%s,%s,%s,%s,%s,%s\n' % (entry[0:]))
            
            
    coords_out.close()
    #print(gene_coords)
    return gene_coords ### A list of selected gene coords

## With Inbuilt fragFASTA functionality - Removes pyfasta dependency - Break coords list and prepare fragmented genome - Lazy file loading i.e. only required chromosomes are loaded and flushed - 
def getFASTALocal(genomeFile,coords):

    # print("Saving genome to 'genomic_seq.fa")
    # sys.exit()
    fastaOut = './genomic_seq.fa'
    fh_out = open(fastaOut, 'w')

    ## Read Genome file - Need lazy loading or chunk reading - Else will fail on large genomes
    fh_in = open(genomeFile, 'r')
    genomeFile = fh_in.read()
    genomeList = genomeFile.split('>')[1:] ## First block is empty
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
            print('Chromosome:',chromo)
            
            gene_seq = chromo[start:end]  ##Get first gene on that strand
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATGC"))
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
         
        ###If chr_id is in chromo_mem i.e. its loaded in memory just splice your sequence      
        
        elif end == '-': ## Till end of chromosome
            print('Fetching gene %s to prepare FASTA file' % (gene))
            gene_seq = chromo[start:]##
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATGC"))
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
            
        else:
            print('Fetching gene %s to prepare FASTA file ' % (gene))
            gene_seq = chromo[start:end]##
            #print(gene_seq)
            if strand == 'C':
                gene_seq_rev = gene_seq[::-1].translate(str.maketrans("TAGC","ATGC"))
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))

    time.sleep(10)
    fh_out.close()
    
    ## Intergenic splicing problem
    ### fasta file empty problem - Fixed
    
    return fastaOut

## This module gets FASTA sequence of supplied coords from genome DB
def getFASTAServer(con,db,gene_coords):
    ####Get the FASTA chromosome and strand wise
    ####NOTE the coords list must be sorted on chr_id and strand as FASTA seq of chr_id and strand is loaded when first instance  encountered. so to elimanate repeated loading sort them
    cur= con.cursor()
    genome_file = './genomic_seq.fa'
    fh_out = open(genome_file, 'w')
    #test_out = open('test_out', 'w')
    chromo_mem = []##used as memory to keep track of Chr+strand in memory, if its in here chromosome will not be read, untill its not here and appended
    for i in gene_coords: ###For every entry either gene or intergenic, use gene_info for just genes
        #print (i)
        chr_id = i[0]
        strand = i[1]
        gene = i[2]
        start = i[3]-1###Adjusted to pick from begining because when spliced it does not include starting position so starting from one poistion back
        end = i[4]
        #print('%s,%s,%s,%s,%s\n' % (chr_id,strand,gene,start,end))
        
        cur = con.cursor()
        if tuple(i[0:2]) not in chromo_mem: ### Add chromosome and strand to list and fetch in memory, will be executed first time a unique chr_id and strand found
            ###as list is sorted on chr_is and strand
            chromo_mem.append(tuple(i[0:2]))###Than append
            #print(chromo_mem)
            
            print ("Reading chromosome:%s and strand: '%s' into memory to splice genes" % (i[0],i[1]) )
            cur.execute("SELECT chromosome FROM %s.chrom_sequence where chr_id = %s and strand = '%s'" % (db, chr_id, strand))###Feed into the memory
            chromo = cur.fetchall()
            # print ("\nHere is the 5' of chromosome: %s\n"% (chromo[0][0][0:100])) ### To get the string out of brackets            
            # sys.exit()

            print('Fetching gene %s to prepare FASTA file' % (gene))
            gene_seq = chromo[0][0][start:end]#Get first gene on that strand
            if strand == 'c':
                gene_seq_rev = gene_seq[::-1]
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
                
            #print (gene_seq)           
                
        else: ## If chr_id is in chromo i.e. its loaded in memory just splice your sequence
            print('Fetching gene %s to prepare FASTA file' % (gene))
            gene_seq = chromo[0][0][start:end]##
            if strand == 'c':
                gene_seq_rev = gene_seq[::-1]
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
            #print (gene_seq)
            #fh_out.write('>%s\n%s\n' % (gene,gene_seq))
    time.sleep(10)
    fh_out.close()
    
    return genome_file

## To efficiently break fasta file you need to distribute size and complete sequences
def fragFASTA(FASTA):
    
    ##

    shutil.rmtree('./genome', ignore_errors=True)
    os.mkdir('./genome')
    
    pattern = ".*\.[0-9].*\.fa" #
    print ("\n***Purging older files***")
    for file in os.listdir():
        if re.search(pattern,file):
            print (file,'is being deleted')
            os.remove(file)
    
    statInfo = os.stat(FASTA)
    filesize =round(statInfo.st_size/1048576,2)
    print('\n**Size of FASTA file: %sMB**' % (filesize))##
    
    if filesize <= splitCutoff: ## 
        fls = []
        fls.append(FASTA)
        print ('No fragmentation performed for file %s' % (fls))
        
    else:
        fh_in = open(FASTA, 'r')
        seq_count = fh_in.read().count('>')
        print('**Number of headers in file: %s**\n'% (seq_count))
        #if genome == 'N':
        if seq_count >= 30: #
            
            if filesize <= 3072:
                splitnum = str(maxHits)
            elif filesize > 3072 and filesize <= 5120:
                splitnum = str(round(maxHits*1.25))
            else:
                splitnum = str(round(maxHits*1.5))        
            
            print ("Size based fragmentation in process for '%s' file" % (FASTA))
            retcode = subprocess.call(["pyfasta","split", "-n", splitnum, FASTA])
            fls = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)] ## file list using regex
            #fls = glob.glob(r'%s.[0-9]{1-3}.fa' % (FASTA.split('.')[0])) ## fragmented file list ##
            #print ('The fragments: %s' % (fls))
               
        
        else: ## 
            splitnum = str(seq_count) ## 
            if fragFasta == 'Y':
                print ("Header based fragmentation in process for '%s' file" % (FASTA))
                retcode = subprocess.call(["pyfasta","split", "-n", splitnum, FASTA])
            fls = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)]

            
            #print ('The fragments: %s' % (fls))
    
    memFile = open('frag.mem','w')
    memFile.write("fasta=%s\n" % (FASTA))
    memFile.write("genomeFeature=%s\n" % (str(genomeFeature)))
    memFile.write("size=%sMB\n" % (filesize))
    memFile.write("frags=%s" % (','.join(fls)))
    memFile.close()
             
    return fls

## Input miRNAs from given file else from server
def miRinput(con):
    miRs = [] ## List to hold miRNAs
    #fh = open(miRNA_file)
    
    #try:## If user specified some local miRNA file in config
    if os.path.exists(miRNA_file):
        fh = open(miRNA_file, 'r')##Test if file exists or not
        miRNA_file_clean = CleanHeader(miRNA_file)## Not required when info is downloaded from server in case of miRPage script
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
        miRFile = miRNA_file_clean
        print ('Total number of miRNAs in given file: %s\n' % (len(miRs)))
        
        fh_out2.close()
        
    #except NameError:## NO Local file specified by User - Fetch  from table
    else:
        print ('\nFetching miRs from the miRNA table')
        miRFile = 'miRinput.fa'
        fh_out = open(miRFile, 'w')## File with miRNAs for troupleshooting
        fh_out2 = open('miRinput_RevComp.fa', 'w')
        
        cur = con.cursor()
    ######################TEST SECTION - COMMENT OUT WHEN NOT TESTING###########################
        #cur.execute("SELECT mirna_name, mirna_precursor_name, score, mature_sequence FROM master.mirna_filter_result where score > 1 group by mirna_name limit 24")
    #####################REAL WORLD - UNCOMMENT WHEN NOT-TESTING#############
        #cur.execute('SELECT mirna_name, mirna_precursor_name, score, mature_sequence FROM master.mirna_filter_result where score > 1 group by mirna_name')##For kevins table
        cur.execute("SELECT mir_name, accession, organism, mir_seq FROM mir_master.mirBASE WHERE mir_name like '%s%%' AND version = %s AND mir_len between 21 and 22" % (miRCode,ver))
    #####################Uncomment when testing is done#########################################
        miRs = cur.fetchall()
        print ('Total number of unique miRs found in miR table: %s' % (len(miRs)))
        
        ##Extract miRname and miR sequence for PARE validation and output in FASTA format as read by cleaveland
        #fh_out = open('miRinput.fa', 'w')
        mirTable = miRTable
        print('****miR info cached****\n')
    
        for miR in miRs:
            mirname = str(miR[0])##Converted to string from unicode, a unicode mirname example: u'meyers-miR544-5p'
            mirseq = str(miR[3]) ##Converted to string from unicode, a unicode mirname example: u'TGAAGATGAAGAAGATGAAGAAGA'
            #print (mirname, mirseq)
            fh_out.write('>%s\n%s\n' % (mirname,mirseq))
            fh_out2.write('>%s\n%s\n' % (mirname,mirseq.translate(str.maketrans("AUGC","TACG"))[::-1]))## make rev comp of miRNA so that it matches the target site in genome rather than mapping miRNA to genome - in target finder file make sure that miRNA is complemented again to main original seq but not direction
                
        fh_out.close()
        fh_out2.close
        
        
    #for i in miRs:
    #    print (i)
        
    return miRs,miRTable,miRFile##miRs holds the list of miRNA name and query where as miRtable holds flag -table name or local

## Attempt3 - Using Bowtie 2 - Goal is to use the same index made for Degradome mapping 
def tarFind3(frag):
    
    file_out = './predicted/%s.targ' % (frag.rpartition('.')[0]) ## Result File
    
    ### Make or select index
    index = "./index/%s_index" % (frag)
    if indexStep == 'Y':
        print('**Creating index of cDNA/genomic sequences:%s\n**' % (index))
        retcode = subprocess.call(["bowtie2-build", frag, index])

    else: ### Check for existing index
        if os.path.isfile(index):
            retcode = 0
        else: 
            print('**Found index of cDNA/genomic sequences:%s\n**' % (index))
            sys.exit()
            
    if retcode == 0: ### Index creation sucessful or index already exists
        print ('Predicting targets for frag:%s using index:%s' % (frag,index))
        nthread2 = str(nthread)
        if TPmode == 'H': ## Heurustic
            print ("You chose 'Heuristic mode' for target identification")
            intervalFunc = str("L,4,0.1")
            minScoreFunc = str("L,-24,-0.5") ###~34.5 - 35 i.e 34 - Stable v1.08
            
            # minScoreFunc = str("L,-2,-0.5") ## 21-24nt -> 16-19nt match -> 32 - 36 score -> penalties can reuce it to zero-2 = -34 cutoff - Test v1.09 - failed
            # readGap = str("22,14") ## Bulge not allowed in gap - Mar-11-15
            
 
            refGap= str("10,8") ## Bulge in miR  + 3MM, 2 seprate or consequite bulges filtered out later - - updated-Mar-23-15 - Can be (10,12) to gain speed - it will effect 1 gap and MM but make 2 gaps impossible
            readGap = str("10,6") ## Bulge in tar + 3MM, 2bulge in tar + 2MM and no bulge 6MM updated-Mar-23-15,
            misPen = str("6,2") ## Mismatch penalty, w/o phred score mx is used i.e. first one - can be changed to (5,2) to improve sensistivity - It will only effect 1 gap scenario and 2 gaps cases will be as is
            matScore = str("0") ## Match score
            ### Changed -D 5 to 6, changed -R 1 to 2 | Jan 13 -D 6 -> -D 3
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 3","-R 2","-N 1","-L 8","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--mp",misPen,"--ma",matScore,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nthread2, "-x", index, "-f" ,"miRinput_RevComp.fa","-S", file_out])
        
        elif TPmode == 'E': ##Exhaustive
            print ("You chose 'Exhaustive mode' for target identification - Please be patient")
            intervalFunc = str("L,2,0.1")
            minScoreFunc = str("L,-24,-0.5") ### ~34.5 - 35 i.e 34 - Stable v1.08
            
            # minScoreFunc = str("L,-2,-0.5") ## 21-24nt -> 16-19nt match -> 32 - 36 score -> penalties can reuce it to zero-2 = -34 cutoff - Test v1.09 - Failed
            # readGap = str("22,14") ## Bulge not allowed in gap - Mar-11-15
                        
            refGap= str("10,8") ## Bulge in miR  + 3MM, 2 seprate or consequite bulges filtered out later - updated-Mar-23-15
            readGap = str("10,4") ## Bulge in tar + 4MM, 2bulge in tar + 3MM and no bulge 6MM - updated-Mar-23-15,
            misPen = str("5,2") ## Mismatch penalty, w/o phred score mx is used i.e. first one
            matScore = str("0") ## Match score

            #### Jan 13 -D 7 -> -D 4 | 
            retcode2 = subprocess.call(["bowtie2","-a","--end-to-end","-D 4","-R 2","-N 1","-L 6","-i",intervalFunc,"--rdg",readGap,"--rfg",refGap,"--mp",misPen,"--ma",matScore,"--min-score",minScoreFunc,"--norc","--no-hd","--no-unal","-p",nthread2, "-x", index,"-f", "miRinput_RevComp.fa","-S", file_out])
        
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

## Parse the bowtie based mapping and generate score - VERY COMPLEX
def tarParse3(targComb):

    ''' Modifying this function is worst nightmare of life - Needs cleaning
    cutoffs w/o bulge or gap - 5MM + 1 wobble 
    bulge in miRNA - 1 bulge + 3MM
    bulge in reference - 1bulge+4mm or 2 consequite bulges +3MM'''
    
    print ('\n**Target prediction results are being generated**')
    ## Input / Output file ######
    print ("File for parsing: '%s' in predicted folder\n" % (targComb))
    fh_in = open(targComb,'r')
    TarPred =  './predicted/%s.parsed.csv' % (targComb.rpartition('/')[-1]) ### Similar to parsed target finder format
    fh_out = open(TarPred,'w')
    fh_out.write('miRname,Target,BindSite,miRseq,tarSeq,Score,Mismatch,CIGAR\n')
    
    #acount = 0
    #for i in fh_in:
    #    print ('Check lines:%s' % (i))
    #    acount +=1
    #print ('Total lines read:',acount)
    #sys.exit()
    
    #### Regenerate Target sequence with all features #####
    acount = 0 ##Total number of interactions from predictions
    parseCount = 0 ## Total number of interactions scores and written to result file
    for i in fh_in:
        # print("\nEntry:",i.strip("\n"))
        acount += 1
        ent = i.strip('\n').split('\t')
        #print('\n%s\n' % ent)
        miRrevcomp = ent[9]                 ## miRNA complemented and reversed to map genome using bowtie. That is target sequence if mismatches and gaps are added
        tarHash = list(miRrevcomp)          ## Strings are immutable convert to list - To rebuilt a traget seq
        # print("\nTarHash",tarHash)
        
        miRrev = miRrevcomp.translate(str.maketrans("TACG","AUGC")) ## Re-translated to get miR but still in reverse orientation - OK
        mirHash = list(miRrev)
        
        #print('Original read mapped i.e miRNA revcomp',miRrevcomp)
        
        ## Gap and Bulges (with reference to miRNA) - Identify gaps and bulges and modify miRNA read used for mapping to regenerate target as well as miRNA
        ## Add gap to target sequence  first to make miR length comparable to target
        gapinfo = ent[5]
        gappos = re.split("[A-Z]",gapinfo)   ## In python format - gap in target seq and bulge in miRNAseq
        gapNuc = re.findall("[A-Z]",gapinfo)
        # print("gappos:",gappos,"| gapNuc:",gapNuc)
        
        ###########################################################################################
        ## SECTION - A - FIND GAPS AND BULGES AND ADD INDICATORS TO MIRNA OR TARGET SEQUENCES
        ###########################################################################################

        posCount = 0
        ## At this point both mirHash and tarHash has same length and perfect complementrity as tarHash is essentially reverse complemented miRNA used for matching
        for x,y in zip(gappos[:-1],gapNuc): ## In gap pos list which is made of alphabet splitting there is always am empty value at end because string has alphabet at last
            # print(x,y)
            if y == 'I':                    ## There was an insertion in miRNA and gap in reference and bulge in miRNA
                for i in range(int(x)):     ## For as many as bulges in miRNA - like 11M 2I 11M              
                    tarHash[posCount] = '-' ## Replace existing nucleotide (from revcomp miRNA) with a gap
                posCount += int(x)          ## This only effects consequitve bulges in miR, if there are multiple insertions like 3I, then "3" needs to be added once and not in every iterneration
            
            elif y == 'D':                  ## There was a deletion in miRNA i.e. gap in miRNA and bulge in reference - In this case length of both miRNA and target will increase
                for i in range(int(x)):             ## For as many as gaps in miRNA                
                    mirHash.insert(posCount,'-')    ## When counted in python insertion will occur after posCount value -  Tested OK
                    tarHash.insert(posCount,'^')    ## Add bulge markers in target sequence as well - Tested OK
                posCount += int(x)                  ## This only effects consequitve gaps in miR, if there are multiple insertions like 3I, then "3" needs to be added once and not in every iterneration
            
            else:
                posCount += int(x)
        
        # print('Target %s seq after manipulation: %s' % (ent[2],''.join(tarHash))) ## Has '-' at gap and '^' at extra nucleotide position (i.e. gap in MiR) -OK
        # print('miRNA %s after after maipulations: %s' % (ent[0],''.join(mirHash))) ## Has '-' at gap pos -OK
        
        #########################################################################################
        ## SECTION -B - GET CORRECT POSITIONS FOR MISMATCHES,GAPS AND BULGES 
        ## AND REGENERATE TARGET BY INSERTING CORRECT NUCLETIDES AT EDITS AND BULGES
        #########################################################################################
        
        ## Mismatches - Identify mismatches and modify miRNA read used for mapping to regenerate target
        misinfoBlock = ent[-2].split(':')[-1] ## Reverse index because XS:i is optional column ## MD:Z:16C3 - these positions are from references - so if there is an insertion/bulge in miRNA i.e. gap that it should be added to these positions before editing miRNA to tar
        # print ('This is the mismatch info block:%s' % (misinfoBlock))
        
        ## Deletion (gap) in miRNAS i.e y= D - which has been added but replace the bulge '^' in target with actual sequence
        ## Should work if there is dletion in miRNA and deletion in traget i.e. two bulges one in miRNA and one in target - NO MM possible
        if '^' in misinfoBlock: 
            misinfo = misinfoBlock.replace('^','')          ## 11^A13 - Here miRNA was 24nt but misinfo shows 25nt as 11+A+13 - Replace the '^' inserted in target with 'A'
            mispos = re.split("[A,T,G,C,N]",misinfo)        ## Found N in one case so included, N confimed in sequence too, will be counted as mismatch
            # print('Mismatch info:%s | Mismatch pos:%s'%(misinfo,mispos))       

            ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
            misposCorrect = []                              ## List hold corrected positions, because edit is a nucleotide next to integers in misPos
            # posIndex = -1                                   ## Index of position for misposCorrect, -1 because after first addition to list it will be incremented to 0
            for x in mispos:
                # print(x)
                ## Assumption: There will be no empty entry in mispos at the begining, because there will be a position to indicate comsequtive edits, others are handled here
                if x:
                    misposCorrect.append(int(x)+1)          ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
                    # posIndex +=1
                else: ## If 'x' is empty like in case of three consequentive gaps in miRNA - Mismatch info:11GTA8 | Mismatch pos:['11', '', '', '8']
                    # y = misposCorrect[posIndex]             ## Get the last corrected position, add one to get position for empty entry
                    # misposCorrect.append(int(y)+1)          ## In consequitve edits this is position just next to last one
                    # posIndex +=1
                    misposCorrect.append(int(0)+1)

            misNuc = re.findall("[A,T,G,C,N]",misinfo)      ## Found N in one case so included, N confirmed in sequence too, will be counted as mismatch
            # print('Misafter:',mispos,'Mispos', misposCorrect,'MisNuc',misNuc)

            ## Replace the nucleotides at bulges(^) and mismatches in tarHash to give actual targets
            ## And also add nucleotide to target (replace ^) if gap in miRNA
            #print('Unedited target:%s-%s' % (''.join(tarHash),len(''.join(tarHash))))
            posCount = 0
            for x,y in zip(misposCorrect,misNuc):
                posCount += x                               ## Convert bowtie positions to python format
                # print(x,y,posCount)
                ## Account for gap before replacing the nucleotide with that in target
                gaps = tarHash[:posCount-1-1].count('-')      ## -1 to convert to python, -1 because - count at positions in target before the current mismatch/bulge position
                tarHash[posCount-1+gaps] = y                  ## Replaced the bulge or mimatch with nucleotide in target - OK
        
        else:   ## Normal i.e y = I - With insertion(bulge) in miRNA and gap in target - OK - What id there is a bulge in miRNA???
                ## In this case miRNA already had inserted nucleltides and target has been added '-' in section-A. Just replace mimatches at correct postions of target
            
            misinfo = misinfoBlock
            mispos = re.split("[A,T,G,C,N]",misinfo)        ## Found N in one case so included, N confimed in sequence too, will be counted as mismatch   
            misposCorrect = [int(x)+1 for x in mispos]      ## Add one to every position to get position where mismatch occured instead of position after which mismatch occured - This is an index and not position
            misNuc = re.findall("[A,T,G,C,N]",misinfo)      ## Found N in one case so included, N confirmed in sequence too, will be counted as mismatch
            # print('Misafter:',mispos,'Mispos', misposCorrect,'MisNuc',misNuc)
            
            ## Count for gaps, since they are added by us and replace MM nucleotides to give actual targets
            posCount = 0
            for x,y in zip(misposCorrect,misNuc):
                posCount += x                                 ## Keep adding the positions, as these are cumulative - MD:Z:2T12C4 - Misafter: ['2', '12', '4'] Mispos [3, 13, 5] MisNuc ['T', 'C'] - Tested OK
                # print(x,y,posCount)
                ## Account for gaps before replacing the nucleotide with that in target
                ## Can give problem if more than one gap? But more than one gap not allowed V07 modification?
                gaps = tarHash[:posCount-1-1].count('-')      ## -1 to convert to python,-1 because - count at positions in target before the current mismatch/bulge position
                tarHash[posCount-1+gaps] = y                  ## TESTED - OK

        tar = ''.join(tarHash).replace("T","U")             ## Target converted to RNA format, will help in catching wobbles ahead
        mir = ''.join(mirHash)
        bindsite = '%s-%d' % (ent[3],int(ent[3])+(len(miRrev)-1))
        
        # print ("Target:%s-%s | miRNA:%s-%s" % (tar,len(tar),mir,len(mir)))
        
        ###################################################################################
        ## SECTION-C - GET POSITIONAL INFORMATION ON GAPS, MISMATCHES, WOBBLES AND BULGES 
        ## AND CHOOSE VALID INTERACTIONS 
        ###################################################################################

        mirGap = [] ## List to hold gaps in miRNA
        tarGap = [] ## List  to hold gap in targets
        mis = []    ## List to hold mismatch position
        wobble = [] ## List to hold Wobble pos
        # print('miRNA: %s\n%s' % (miRrevcomp[::-1],miRrevcomp[::-1].replace("T","U") ))

        ## Read from mapping file -> miRrevComp -> uncomplement -> miRrev -> morhash-> miR 
        valid = 1 ## Validity flag [0] - Invalid and [1] - valid, if more then 1 bulges in miR or 2bulges in tar or bulges in 9,10,11 = invalid
        nt_cnt = 1 ## Keep track of actual position,
        for x,y in zip(mir.translate(str.maketrans("AUGC","UACG",))[::-1],tar[::-1]): ## Orientation changed to read from 5' miRNA - OK
            
            #print(miRrev[::-1][nt_cnt-1],x,y)## Print miRNA, rev complemmnetry miRNA used for matching, target
            # if x == '-' or y == '-' or x == '^' or y == '^':
            #     #print('gap')
            #     gap.append(nt_cnt)
            #     if y == '-':
            #         nt_cnt+=1

            if x == '-' or x == '^': ## Don't think '^' would be here, since it has been replaced with nucleotide
                #print('miRNA gap')
                mirGap.append(nt_cnt)
                if nt_cnt == 9 or nt_cnt == 10 or nt_cnt == 11:
                    # print("@miRNA has gap/bulge in 9th and 10th position")
                    valid = 0 ## This is invalid as it has bulge in miRNA at 9th and 10th pos, this interaction will be skipped

            elif y == '-' or y == '^':
                #print('target gap')
                tarGap.append(nt_cnt)
                if y == '-': ## Don't think '^' would be here, since it has been replaced with nucleotide
                    nt_cnt+=1
                    if nt_cnt == 9 or nt_cnt == 10 or nt_cnt == 11:
                        # print("@target has gap/bulge in 9th and 10th position")
                        valid = 0 ## This is invalid as it has bulge in miRNA at 9th and 10th pos, this interaction will be skipped

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
        # print('MismatchList:%s | miRGapList = %s | tarGapList = %s | WobbleList = %s' % (mis, mirGap, tarGap, wobble)) ## Poistion of mismatch gap and wobble

        ## Check if there are more then mpermitted gaps in miRNA and target
        if len(mirGap) > 2:
            # print("@miRNA %s has more then two gaps" % (ent[0]))
            valid = 0
        elif len(tarGap) > 1:
            # print("@target has more then one gap")
            valid = 0
        elif len(mis) > 5:
            # print("More then 5 mismatches not allowed")
            valid = 0
        elif len(mis) + len(wobble) > 6: 
            # print("Six edits are not allowed - It's too much")
            valid = 0
        else:
            pass

        ## Decide to report this interaction if its valid
        ## Validity flag, if more then 1 bulges in miR or 2bulges in tar or bulges in 9,10,11 = invalid
        if valid == 0: 
            # print("Skipping this entry - It's biologically invalid\n")
            continue ## Go to main for loop
        else:
            ## Go for scoring
            pass

        gap = mirGap+tarGap ## Combine gaps or scoring
        
        #####################################################################################
        ## SECTION-D - SCORE THE INTERACTIONS 
        #####################################################################################

        score = 0   ## Initialize
        #print (mis)
        if TPscore == 'S': ## Allowed 3 MM, 2 Wob, 1 Gap
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
        
        ## Correctly output mismatches - If there is no mismatch then misinfo is just the length of match like - 15G2A1T1 (if there is mismatch) and 21(if no mismatch)
        ## Since the second ouput which originally is part of mapping file confussing in mismatch column - take mismatch info from 'mis' list
        
        if mis or gap: ## If list of mismathes has positions of mismatches, then output the block from mapping file
            mismatches = misinfo
        else:
            mismatches = 0

        #print(ent[0],ent[2],bindsite,miRrev,tar,score,misinfo,gapinfo)## MiRname, Tarname, mirSeq,Taerseq,binding site
        fh_out.write('>%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[2],bindsite,mir.replace("U","T"),tar.replace("U","T"),score,mismatches,gapinfo))
        parseCount  += 1
    
    
    print("Total number of interactions from 'miRferno':%s AND total interactions scored: %s" % (acount,parseCount))
    fh_in.close()
    fh_out.close()

    return TarPred

## Get PARE reads from the server - Can we paralleize this module?
def TagAbundanceFile(con,PAREdb,libs):

        for alib in libs:##For all the libraries
            print ('\nPreparing PARE reads file for library %s ' % (alib))
        #print (lib[0])
            #print ('Caching tag and count information from server for PARE alib %s' % (alib[0]) )
            cur = con.cursor()
            cur.execute("SELECT tag, norm from %s.run_master WHERE lib_id = %s AND (hits between 1 and 20)" % (PAREdb,alib))
            lib_info = cur.fetchall()
            print('These are the tags:',lib_info[:10])
            
            fh_out = open('%s' % (alib), 'w')##Naming file with lib_ids name
            for ent in lib_info:## All the PARE tags in a library     
                fh_out.write('%s\t%s\n' % (ent[0],ent[1]))
            
            fh_out.close()

## Convert tag count files to fasta file without multiplying tag count to read - So every read is uniq
def tag2FASTA2(lib):
    print("'%s' tag count file being converted to FASTA format" % (lib))
    fh_in = open(lib,'r') ### Open tag count file
    fh_out = open('./PARE/%s_PARE_tags.fa' % (lib), 'w')##Naming file with lib_ids name, 0000 means all the libs
    tag_num = 1 ### for tag naming that will be used in FASTA format
    for tag in fh_in:##All the entries of the library
        #print(tag.strip('\n').split('\t'))
        ent = tag.strip('\n').split('\t')
        tag = ent[0]
        if len(tag) >= tagLen: ##Length of tag specified by the user
            fh_out.write('>%s\n%s\n' % (tag_num, tag[:tagLen]))
            tag_num += 1
        else:
            #print ('Length is not 20nt')
            pass
    fh_out.close()

## Map degradome to the transcripts
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
    #retcode2 = subprocess.call(["/home/bioinfo/tools/bowtie2.1/bowtie2", "-a","--norc","--no-hd", "-t","-p",nthread2, "-f", indexLoc, dd_file,"-S",map_out]) ### Default mode is end-end
    retcode2 = subprocess.call(["bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "--no-unal", "-t","-p",nthread2, "-x", indexLoc,"-f", dd_file,"-S",map_out]) ###Score min sets limimum score to '0'. In ened to end mode that means no mismatch/gap
    
    ## Optimized | no-unal was the culprit | 
    #retcode2 = subprocess.call(["/home/bioinfo/tools/bowtie2.1/bowtie2", "-a", "--end-to-end", "-D 1", "-R 1", "-N 0", "-L 20", "-i L,0,1","--score-min L,0,0","--norc","--no-head", "-t","-p",nthread2, "-f", indexLoc, dd_file,"-S",map_out]) ###Score min sets limimum score to '0'. In ened to end mode that means no mismatch/gap
    
    if retcode2 == 0:## The bowtie mapping exit with status 0, all is well
        print('\nDegradome from PARE lib: %s mapped to cDNA/Trascript file' % (templib)) 
    else:
        print ("There is some problem with mapping of PARE lib: %s to cDNA/genomic seq index" % (templib))
        print ("Script exiting.......")
        sys.exit()

## Combine Files given as list - Make it an independent function
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

## Module to clean headers of cDNA/transcript and/or miRNA input file
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

## 
def PPBak(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(module, alist)

##
def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = round((nproc/int(nthread))+1) ## 1 added so as to avoid 0 processor being allocated in serial mode
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

##  
def PPmultiple(module,alist1,alist2):
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(lambda args: module(*args), alist2)

##
def PPResults(module,alist):
    npool = Pool(int(nproc))    
    res = npool.map_async(module, alist)
    if(res.get() == []):
        print("YEP!")
    results = (res.get())
    npool.close()           ### Added by Reza to close hanging 
    return results

##  
def feed(queue, parlist):
    print ('Feeder function started')
    for par in parlist:
        print ('Echo from Feeder: %s' % (par))
        queue.put(par)
    print ('**Feeder finished queing**')

##
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

##
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

## Parallelized version of RevMapCoord - introduced in version 2.04 - REGEX dropped and Dictionary introduced v2.05
def RevMapCoord(ent): ####DOES THIS MODULE WORKS IF USER USES MULTIPLE LIBRARIES - Parallelize
    
    ## Create a dictionary from list of coords to be searched later
    ## Gene_coords structure: 1, 'c','AT1G01020', 5928, 8737, protein_coding   
    
    print (ent)
    gene_name = ent[1]
    cleave_site = int(ent[8])
    bind_site = ent[2].split('-')
    
    ## Reverse map co-ordinates ##########################################################
    if Local == 'N':
        print ('\n**Reverse mapping of Co-ordinates will be performed - Web analysis**')
        ###Check whether gene is from reverse or posiive strand by memebr ship test on dictionary
        if gene_name in coord_dict_wat:
            print ('Entry: %s in positive strand: %s' % (ent[0:4],coord_dict_wat[gene_name]))
            geno_start = coord_dict_wat[gene_name][1]###Use for reverse mapping of postive genes
    
            #print('Looking for chr_id')
            chr_id = coord_dict_wat[gene_name][0]
            #print('chr_id found')
            strand = 'w' ## AVlilable in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_wat[gene_name][2] ## Gene type
            new_cleave_site = (int(geno_start)-1)+int(cleave_site)###1 is reduced to give correct positions
            new_bind_site_start = (int(geno_start)-1)+int(bind_site[0])
            new_bind_site_end = (int(geno_start)-1)+int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
        else:
            print ('Entry: %s in reverse strand: %s' % (ent[0:4],coord_dict_crick[gene_name]))
            geno_end = coord_dict_crick[gene_name][1]###use for reverse mapping of negative genes
            #print('Looking for chr_id')
            chr_id = coord_dict_crick[gene_name][0]
            #print('chr_id found')
            strand = 'c' ## Available in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_crick[gene_name][2] ##Gene type
            new_cleave_site = (int(geno_end)+1)-int(cleave_site)###1 is added to give correct positions
            new_bind_site_end = (int(geno_end)+1)-int(bind_site[0])###As the sequence was reversed before TF and CL, their binding start and end direction has also changed - Verified-OK
            new_bind_site_start = (int(geno_end)+1)-int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
   
    else: ## 'Y' i.e Local analysis
        print ('**No Reverse mapping of Co-ordinates will be performed - Local analysis**')

    ## Get window abundances and annotations ###############################################
    
    small_win = 0###Small window abundance for this validated gene
    large_win = 0###Large window abundance for this validated gene
    print ('*Acquiring PARE abundances in small and large window............*')

    if Local == 'N':
        ## Fix the gene name on reverse strand and than you can use gene name as it is by removing _up/_down
        cur = con.cursor() ## Connect to data server with table
        gene_query = gene_name
        print ('Gene name:', gene_query,'Large win start:',int(new_cleave_site)-5,'Cleave site:',new_cleave_site,'Large win end:',int(new_cleave_site)+5)
        
        ## Making just single query for region of large window, can also make two queries for large and small window and get sum(norm_sum) in this case no need to calculate small win and large win abundances as below
        ## Single query because sometimes server is busy or slow
        if gtype == 'inter':## Gene name not used
            cur.execute("SELECT position,norm_sum FROM %s.tag_position where chr_id = %s AND strand = '%s' AND position between %s and %s" % (PAREdb,chr_id,strand,int(new_cleave_site)-5,int(new_cleave_site)+5)) ## No gene name
            window = cur.fetchall()
            #print(window,'\n')
            genefunc = 'NA'
        
        else: ## Gene name used
            cur.execute("SELECT position,norm_sum FROM %s.tag_position where gene = '%s' AND chr_id = %s AND strand = '%s' AND position between %s and %s" % (PAREdb,gene_query,chr_id,strand,int(new_cleave_site)-5,int(new_cleave_site)+5)) ## Clear table before writing
            window = cur.fetchall()
            
            ## Gene Annotation
            cur.execute("SELECT gene,title FROM %s.gene_master where gene = '%s'" % (GenomeDB,gene_query))
            info = cur.fetchall()
            # print("Anno:",info)

            if info[0][1] != None:
                # print (info[0][1])
                genefunc = info[0][1].replace(",",";")
            elif info[0][1] == None:
                # print ('Title not found for: %s' % (gene_query))
                genefunc = 'NA'
                pass
            else: ## Not found
                # print ('Title not found for: %s' % (gene_query))
                genefunc = 'NA'
                pass


        if extraAnno == 'Y':
            ## Phased annotation - need phased_loci table in genome DB
            cur.execute("SELECT cluster_id,len FROM ASPARAGUS_UGA1_genome.phased_loci where chr_id = %s and (%s between start and end)" % (chr_id,int(new_cleave_site)))
            info2 = cur.fetchall()
            # print('Phasing Anno:',info2)
            if info2:
                phasing = '-'.join(str(i) for i in info2[0])
                # print('Phasing Anno:',phasing)
            # elif info2[0][1] == None:
            #     phasing = 'No'
            else:
                phasing = 'No'

            ## Inverted repeat annotation - Need chr_inverted table in genome DB
            cur.execute("SELECT * FROM ASPARAGUS_UGA1_genome.chr_inverted where chr_id = %s and ( (%s between start1 and end1)  or (%s between start2 and end2) )" % (chr_id,int(new_cleave_site),int(new_cleave_site)))
            info3 = cur.fetchall()
            # print("Inverted anno:",info3)

            if info3:
                inverted = 'IR'
            # elif info3[0][1] == None:
            #     inverted = 'No'
            else:
                inverted = 'No'

            
        ## Compute small window and large window using the data collected above
        for entry in window:
            if int(entry[0]) <= int(new_cleave_site+1) and int(entry[0]) >= int(new_cleave_site-1):## Range of Small Window, Every 5' mapping on specific gene is read and if 5'end is mapped within small window than added
                small_win += (int(entry[1]))
    
            if int(entry[0]) <= int(new_cleave_site+5) and int(entry[0]) >= int(new_cleave_site-5):## Range of Large window, GET LARGE AND SMALL WINDOW IN SAME LOOP
                large_win += (int(entry[1]))
        print('Small Window:',small_win, "| Large Window:",large_win)

    
    ##else: ## Local
    ##    print ('Looking up gene in density file............')
    ##    try:
    ##        found = den_dict[gene_name]
    ##        print ('Gene found in density file\n')
    ##        for i in found:###In gene summarization record first two entries are gene_name and size of gene, last entry is a blank line when record is fetched (not in file)
    ##            locus = i.split()##position in summarization entry along with raw read, repeat normalized read and unique reads mapped, and category
    ##
    ##            
    ##            if locus[0].isdigit() and int(locus[0]) <= int(cleave_site+3): ###Fixes BUG 14-3-13, found one line in summarized file as ['#', 'Transcriptome=./genomic_seq.fa'] instead of like ['12101480', '2', '1', '0', '2'] AND Reads file till large window limit+ 1 extra
    ##                #print ('Entry :',locus, 'Cleave Site :',cleave_site)
    ##                if int(locus[0]) <= int(cleave_site+1) and int(locus[0]) >= int(cleave_site-1):####Range of Small Window###Every 5' mapping on specific gene is read and if 5'end is mapped within small window than added
    ##                    #print ('Using %s for small window' % locus[0])
    ##                    small_win += (int(locus[1]))
    ##                #else: ###if the 5' end is not within the window - pass
    ##                #    pass
    ##                #
    ##                if int(locus[0]) <= int(cleave_site+7) and int(locus[0]) >= int(cleave_site-7):###Range of Large window###GET LARGE AND SMALL WINDOW IN SAME LOOP
    ##                    #print ('Using %s for large window' % locus[0])
    ##                    large_win += (int(locus[1]))
    ##                
    ##                #else: ###if the 5' end is not within the window - pass
    ##                #    pass
    ##            else: ## That means present position is more than extent of large window
    ##                #print ('\nWe have passed the large window range so breaking\n')
    ##                break
    ##        
    ##    #else: ###gene_name doesn't have a record in sumamrization file - tags in some lib mapped to gene but in '0_density' tags were missng and not mapped to any such gene
    ##    ##So Gene is missing from normalized tags summarization file a.k.a 0_density
    ##    except KeyError:
    ##        print ('\n ** Warning - There is no entry for %s found in PARE summarization**\n' % (gene_name))
        
    ##################################################################################################################
        
    ## BugFix -  13-3-13 - In lib based validation, target is validated just raw reads mapped by cleaveland and sometimes there is no unique and normalized read mapped at corresponding site in '0_density' file so skip those
    if small_win != 0:
        if Local == 'N':
            if extraAnno == 'Y':
                rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' % (','.join(ent),str(chr_id), strand, new_bind_site, new_cleave_site, small_win, large_win, round(small_win/large_win,2),genefunc,phasing,inverted))
                print ("Server:",rev_mapped_entry)
            else:
                rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s,%s' % (','.join(ent),str(chr_id), strand, new_bind_site, new_cleave_site, small_win, large_win, round(small_win/large_win,2),genefunc))
                print ("Server:",rev_mapped_entry)


        else: ## Local analysis has no Reverse mapping 
            rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%f\n' % (ent[0][1:],ent[1],ent[2], ent[3],ent[4],ent[6],ent[7],ent[5],small_win,large_win,round(small_win/large_win,2)))
            print ("Local:",rev_mapped_entry)
    
    else:
        ##If the first result entry has 0 abundance for small window than no 'rev_mapped entry' variable and nothing to return, therefore added 'Rev_mapped_entry' with ratio set to 0
        #rev_mapped_entry = (('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0][1:],gene_name,new_bind_site, new_cleave_site,ent[4],ent[6],ent[7],ent[5],small_win,large_win, 0)))
        ##Or rev_mapped_entry just mentions the bug that caused blank entry which can be filtered later
        print('E13-3-13')
        rev_mapped_entry = ('E13-3-13')
        pass
    # print ("mapping success")
    return rev_mapped_entry

## Generate a single results file by combining lib-wise results
def resultUniq():
    
    """Read revmapped files, remove header, combine files
    parse on pval,Uniq on miRname, tarname and cleavesite"""
    
    if Local == 'Y':
        revmapped_fls = [file for file in os.listdir('./output') if re.search(r'validated', file)]
        print ('Validated files:',revmapped_fls)
        print ('\nCombining all the validated files for uniq\n')
        revmappedComb = './output/ALL_validated.csv'

    else :        
        revmapped_fls = [file for file in os.listdir('./output') if re.search(r'revmapped.csv', file)]
        print ('revmapped files:',revmapped_fls)
        print ('\nCombining all the revmapped files for uniq and table update\n')
        revmappedComb = './output/ALL_revmapped.csv'
    

    fh_out = open(revmappedComb ,'w')
    
    ##Combine files
    for x in revmapped_fls:
        print (x)
        revamppedfile = open('./output/%s' % (x), 'r')
        header = revamppedfile.readline() ## Use later
        data = revamppedfile.read()
        revamppedfile.close()
        fh_out.write(data)
    fh_out.close()
    
    ## Read and parse All result file:
    fh_in=open(revmappedComb, 'r') 
    parsed_in = [line.strip('\n').split(',') for line in fh_in]
    parsed_in.sort(key=itemgetter(14))##Sorted on corrected p-value beofre removing redundant so that first entry stored is best among other same interactions across library
    
    
    if Local == 'Y':
        uniqRevmapped = './output/ALL_validated_Uniq.csv'
    else:
        uniqRevmapped = './output/ALL_revmapped_Uniq.csv'

    fh_output2=open(uniqRevmapped, 'w')
    
    
    added_keys=set()## A set to store first 3 elements from input file: miRNA-DNA, chr# and cleavage site and than use it to compare further entries in file
    parsed_out_count=0## To keep count of unique entries
    
    print('\nRemoving redundant entries aka more than one instance of an miRNA (PARE Validated)')
    
    if Local == 'Y':
        fh_output2.write('miRNA,Target,localBindSite,mirSEq,tarSeq,Score,Mismatch,CIGAR,cleavage position,PARE reads,winAbundance(10nt),localWinRatio,category,p-value,corrected p-val,Chr\n')
        for ent in parsed_in:
            genename = ent[1].split('.')[0]## To avoid different variations of same gene to be counted as uniq
            lookup=tuple((ent[0],genename,ent[1],ent[8]))## miR name + Target Gene + position of cleavage on gene, earlier miR sequence was used instead of miR name but that removes miR belonging to same family as they have same sequence
            if lookup not in added_keys:
                fh_output2.write('%s\n' % (','.join(ent)))
                parsed_out_count+=1
                added_keys.add(lookup)## Once a new entry is found it is recorded so as to compare and neglect further entries
            else:
                pass
    else:
        fh_output2.write('miRNA,Target,localBindSite,mirSEq,tarSeq,Score,Mismatch,CIGAR,cleavage position,PARE reads,winAbundance(10nt),localWinRatio,category,p-value,corrected p-val,Chr,Strand,BindSite,CleaveSite,SmallWin,LargeWin,Ratio,Func\n')
        for ent in parsed_in:
            genename = ent[1].split('.')[0]## To avoid different variations of same gene to be counted as uniq
            lookup=tuple((ent[0],genename,ent[16],ent[18]))## miR name + Target Gene + Strand + position of cleavage on gene, earlier miR sequence was used instead of miR name but that removes miR belonging to same family as they have same sequence
            if lookup not in added_keys:
                fh_output2.write('%s\n' % (','.join(ent)))
                parsed_out_count+=1
                added_keys.add(lookup)## Once a new entry is found it is recorded so as to compare and neglect further entries
            else:
                pass

 
    fh_output2.close()
    
    return uniqRevmapped

## Module that uploads data to already made table in a DB
def TableUpload(con3,uniqRevmapped,mir_tab_name,miRFile):###Con is connection and res_upload is file from last module which has genomic coords and need to be upload
    ##Upload scoring_input_extend to table - please make sure that columns are in file are in order of the columns in the table
    
    cur = con3.cursor()##Connect to destination server with table
    #res_file = 'scoring_input_extend_upload'
    
    try:## If user specified some local miRNA file in config
        fh = open(miRFile)##Test if file exists or not
        mir_tab_name = 'None'
    except NameError:
        pass
    
    if TableWipe is 'Y': ## If the Global setting is 'Y' than table will be dropped - Be careful   
        print ('\n**Clearing the table before updating.....**')
        cur.execute("TRUNCATE TABLE mir_master.mirPARE_v3")## Clear table before writing
        #con2.commit()
        print ('**Table cleared successfully, update in process**\n\n')
        
    ##Current implementation of mysql.connector does not support instant upload by local file - see ConnectToDB() module for implementation
    ##Original query - LOAD DATA LOCAL INFILE "./scoring_input_extend" INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ',';
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    ##So fill table on row by row basis
    ##Check for miRNA table variable if not than Local file used else - public/priv
    
    
    add_row = "INSERT INTO mir_master.mirPARE_v3(mir_name,target_name,chr_id,strand,binding_site,cleavage_site,target_score,cat,mir_seq,target_seq,p_value,small_win,large_win,win_ratio,pval_corr,mismatch,CIGAR,PAREdb,miRdb_type) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
    
    #test_data =[["Test1", "AT2G28350.1", "1329-1349", "1340", "2", "ACCGUAUGUCCCUCGGUCCGU", "AGGAAUACAGGGAGCCAGGCA", "0.000314"],["ath-miR160b-5p", "Test", "1329-1349", "1340", "2", "Test", "Test", "0.000314"]]
    print (PAREdb,mir_tab_name,uniqRevmapped)
    fh_in2 = open(uniqRevmapped, 'r')
    fh_in2.readline() ##Waste header
    
    for row in fh_in2:
        res_entry = row.split(',')
        print(res_entry[0][1:], res_entry[1], res_entry[15], res_entry[16], res_entry[17], res_entry[18], res_entry[5], res_entry[12],res_entry[3],res_entry[4],res_entry[13],res_entry[19],res_entry[20],res_entry[21],res_entry[14],res_entry[6],res_entry[7],PAREdb,mir_tab_name)
        res_upload = (res_entry[0][1:], res_entry[1], res_entry[15], res_entry[16], res_entry[17], res_entry[18], res_entry[5], res_entry[12],res_entry[3],res_entry[4],res_entry[13],res_entry[19],res_entry[20],res_entry[21],res_entry[14],res_entry[6],res_entry[7],PAREdb,mir_tab_name)
        cur.execute(add_row,res_upload)        
        con3.commit()
        
    cur.close()

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

    cleaveLocations = {}
    cleaveLocations[21] = [9, 10, 11]
    cleaveLocations[22] = [9, 10, 11]
    cleaveLocations[23] = [9, 10, 11, 12]
    cleaveLocations[24] = [9, 10, 11, 12]
    cleaveStandard = [9, 10, 11]

    validatedTargets = []
    for target in targetFinderList:
        gene = target[1]
        # Get the start location of the target
        end = int(target[2].split('-')[1])
        # Get the length of the mIRNA
        length = int(len(target[3]))

        if(gene in PAGeDict.keys()):
            currDict = PAGeDict[gene]
            cleavageSite = []
            location = 0
            targetAbundances = []
            targetCategories = []
            targetLocations = []
            #
            #
            #

            # If the length of the miRNA is in our dictionary and the cleave
            # location exists in the dictionary, add it to the target lists
            if(standardCleave == 'N' and length in cleaveLocations.keys()):
                for cleaveLocation in cleaveLocations[length]:
                    if(str(end-cleaveLocation) in currDict):
                        targetAbundances.append(currDict[str(end -
                            cleaveLocation)][0])
                        targetCategories.append(currDict[str(end -
                            cleaveLocation)][1])
                        targetLocations.append(end - cleaveLocation)

            # If the length of the miRNA is not in our dictionary, we will 
            # just investigate the 10th 11th and 12th positions 
            else:
                for cleaveLocation in cleaveStandard:
                    if(str(end-cleaveLocation) in currDict):
                        targetAbundances.append(currDict[str(end -
                            cleaveLocation)][0])
                        targetCategories.append(currDict[str(end -
                            cleaveLocation)][1])
                        targetLocations.append(end - cleaveLocation)
                

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

    # Print("Making PAGe dictionary %s" % bowtieFilename)
    # Initialize dictionary for bowtie entries
    bowtieDict = {}
    
    # Append the directory to the bowtieFilename
    bowtieFilename = 'dd_map/' + bowtieFilename
    bowtieFile = readFile(bowtieFilename)

    if(bowtieFilterSwitch == 'Y'):
        # Loop through the entire bowtie files to load these into the dictionary
        for entry in bowtieFile:
            if(entry.split('\t')[4] != '255'):
                pass
            else:
                ## Store gene, location and sequence from each entry into variables
                gene = entry.split('\t')[2]
                location = entry.split('\t')[3]
                sequence = entry.split('\t')[9]

                ## Try to identify if the sequence is a key in the bowtie dictionary
                try:
                    bowtieDict[sequence]
                ## If not, create one with a value of an empty dictionary
                except:
                    bowtieDict[sequence] = {}

                ## Try to identify if gene is a key in the nested dictionary
                try:
                    bowtieDict[sequence][gene]

                ## If not, create one with a value of an empty list 
                except:
                    bowtieDict[sequence][gene] = []

                ## Append the location of the sequence in the gene to the list
                bowtieDict[sequence][gene].append(location)        

    else:
        for entry in bowtieFile:
            ## Store gene, location and sequence from each entry into variables
            gene = entry.split('\t')[2]
            location = entry.split('\t')[3]
            sequence = entry.split('\t')[9]

            ## Try to identify if the sequence is a key in the bowtie dictionary
            try:
                bowtieDict[sequence]
            ## If not, create one with a value of an empty dictionary
            except:
                bowtieDict[sequence] = {}

            ## Try to identify if gene is a key in the nested dictionary
            try:
                bowtieDict[sequence][gene]

            ## If not, create one with a value of an empty list
            except:
                bowtieDict[sequence][gene] = []

            ## Append the location of the sequence in the gene to the list
            bowtieDict[sequence][gene].append(location)

    PAGeDict = {}
    allHits = []

    ## Loop through all tags and their hit counts
    for entry in tagCountFile:
        ## sequence and hits used multiple times, so store as variables.
        sequence = entry.split('\t')[0]
        hits = entry.split('\t')[1]

        ## If the tag exists in the bowtie file, assign the hits to the 
        ## data structure. Otherwise, a pass must be made
        try:
            ## Loop through all gene locations
            for key in bowtieDict[sequence].keys():
                ## Test to see if key exists in PAGeDict yet. If not, create it
                # 
                try:
                    PAGeDict[key]
                except:
                    PAGeDict[key] = {}
                
                ## Loop through all locations of the gene in which the tag maps
                ## to
                for location in bowtieDict[sequence][key]:
                    # 
                    PAGeDict[key][location] = int(hits)
                   
                    # Append the hits to the hits if hits > 2
                    if(int(hits) > 2):
                        allHits.append(int(hits))

        except:
            pass

    return(PAGeDict, allHits)

def unambiguousBaseCounter(transcriptomeFilename, tagLen):
    """Get the counts of ambiguous bases in the transcriptome file as well
       as counts of ambiguous bases that are within the ends of the
       transcriptome +/- the tagLen.

    Args:
        transcriptomeFilename: Name of the trnascritome fasta file
        tagLen: Number of bases that an N is allowed to be away from the ends
            of the gene in order to be counted

    Returns:
        Total number of ambiguous bases and ambiguous bases tagLen-bp away from
        the ends of the gene.

    """
    baseCounts, baseCountsOffTagLen = 0, 0

    transcriptomeFile = readFile(transcriptomeFilename)
    for i in range(1,len(transcriptomeFile),2):
        currentLine = transcriptomeFile[i]
        baseCounts += len(currentLine) - currentLine.count('N')
        baseCountsOffTagLen += (len(currentLine) - 2 * tagLen) - currentLine[
            tagLen:len(currentLine)-tagLen].count('N')

    f_output = open('baseCounts.mem', 'w')
    f_output.write(str(baseCounts) + '\n' + str(baseCountsOffTagLen))
    f_output.close()

def writePAGeFile(PAGeDict, mode, allHits, baseCounts, baseCountsOffTagLen,
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
    # 
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

    # Sort genes so that genes are in alphabetical order
    for gene in sorted(PAGeDict.keys()):
        # Dictionary for each gene which records all locations
        # with the same abundance and categories for output
        hitsDict = {}
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
                #f.write('%s\t%s\t%s\n' % (str(location), str(hits), category))

                # If the abundance already exists, append the location to the
                # location list in the tuple. Otherwise, create the entry and
                # add the category
                try:
                    hitsDict[hits][0].append(location)
                except:
                    hitsDict[hits] = ([location], category)

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
                #f.write('%s\t%s\t%s\n' % (str(location), str(hits), category))

                # If the abundance already exists, append the location to the
                # location list in the tuple. Otherwise, create the entry and
                # add the category
                try:
                    hitsDict[hits][0].append(location)
                except:
                    hitsDict[hits] = ([location], category)

        for hits in sorted(hitsDict.keys(), key=int):
            f.write('%s\t' % hits)
            f.write(','.join(sorted(hitsDict[hits][0], key=int)))
            f.write('\t%s\n' % hitsDict[hits][1])
            #f.write('%s\t%s\t%s\n' % (hits, sorted(hitsDict[hits][0], key=int), hitsDict[hits][1]))
            #for i in range(len(hitsDict[hits][0])):
            #    loc = sorted(hitsDict[hits][0])[i]
            #    if(i == len(hitsDict[hits][0]) - 1):
            #        print(i)
            #        f.write('%s\t%s;\n' % (loc, hitsDict[hits][1]))
            #    else:
            #        f.write('%s,' % loc)

        hitsDict.clear()
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

    # Code to find the corrected p-value
    pValueList = []
    # Get indicies of various elements
    pValueIndex = len(validatedTargets[0]) - 1
    categoryIndex = len(validatedTargets[0]) - 2
    windowRatioIndex = len(validatedTargets[0]) - 3
    for target in validatedTargets:
        pValueList.append(target[pValueIndex])
    
    # Code to perform noise corrected pValue
    for i in range(len(validatedTargets)):
        validatedTargets[i].append(float(validatedTargets[i][pValueIndex]) / 
            float(validatedTargets[i][windowRatioIndex]))

    # get corrected p-value index
    correctedPValueIndex = len(validatedTargets[0]) - 1

    f = open(outputFile,'w')
    
    
    f.write(header + ',cleavage position,PARE reads,10 nt window abundance,'\
        'PARE reads/window abundance,category,p-value,corrected p value\n')

    # noiseFilter requires that p value be <= .25 and window ratio be
    # >= .25 for a tag to be considered
    if noiseFilter == 'Y':
        for target in validatedTargets:
            # Include any target with pvalue <= .25 and with window ratio
            # >= .25
            if((float(target[pValueIndex]) <= .25) and (float(target[
                    windowRatioIndex]) >= .25)):
                # If category 4s are not to be included and the category
                # of the current category is 4, then omit it.
                if(cat4Show == 'Y' and (int(target[categoryIndex]) == 4)):
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
            if(float(target[pValueIndex]) < .5):
                # If category 4s are not to be included and the category
                # of the current category is 4, then omit it.
                if(cat4Show == 'Y' and (int(target[categoryIndex]) == 4)):
                    pass
                else:
                    for j in range(len(target)):
                        f.write(str(target[j]))
                        # If at last eleent, write new line
                        if(j == correctedPValueIndex):
                            f.write('\n')
                        # Otherwise, write a comma
                        else:
                            f.write(',')

#### MAIN FUNCTION #############################

def main():
    
    global con
    
    # print ("Is atul the owner?")
    # sys.exit()

    ## FASTA GENERATION ###################
    if Local == 'Y':
        if generateFasta == 'Y':
            coords = extractFeatures(genomeFile,gffFile) ## Extracts Coords from GFF3
            fastaOut = getFASTALocal(genomeFile,coords) ##Creates FASTA file
            global tagLen ## Required later for tag2FASTA step as well
            unambiguousBaseCounter(fastaOut,tagLen)
        else:
            print("\nThe input FASTA file is considered 'as is' for analysis\n")
            fastaOut = genomeFile ### Make it better
    else:
        if generateFasta == 'Y':
            con = ConnectToDB(dataserver,0)
            coords = extractCoords(con,GenomeDB)###Get selected genes and all intergenic coords - all coords all also fetched too
            fastaOut = getFASTAServer(con,GenomeDB,coords)###Extracts FASTA of supplied coords
            unambiguousBaseCounter(fastaOut,tagLen)
        else:
            print("\nThe input FASTA file is considered 'as is' for analysis\n")
            con = ConnectToDB(dataserver,0)
            coords = extractCoords(con,GenomeDB) ##For reverse mapping
            fastaOut = genomeFile ### Make it better
        
    ########################################################################################################
    ### FRAGMENTATION ######################################################################################
    
    ## Script Timer
    runLog = 'runtime_%s' % datetime.datetime.now().strftime("%m_%d_%H_%M")
    fh_run = open(runLog, 'w')
    print('TPmode: %s | TPscore: %s | Uniq filter: %s\n' % (TPmode,TPscore,bowtieFilterSwitch))
    fh_run.write ('Libs: %s\n' % (','.join(libs)))
    fh_run.write('TPmode:%s | TPscore: %s | Uniq filter:%s\nLibs:%s\nGenomeFile:%s | GenomeFeature:%s\n' % (TPmode,TPscore,bowtieFilterSwitch, ','.join(libs),genomeFile,genomeFeature))
    FragStart = time.time()
    
    if fileFrag == 'Y':
        start = time.time()###time start
        fragList = fragFASTA(fastaOut)##
        end = time.time()
        print ('fileFrag time: %s' % (round(end-start,2)))
    else:
        fragMem = open('frag.mem','r')
        fragRead = fragMem.readlines()
        print("Frags from earlier processed file: '%s' will be used in this run" % (fragRead[0].strip('\n').split('=')[1]))
        for i in fragRead:
            akey,aval = i.strip('\n').split('=')
            if akey == "genomeFeature":
                if int(aval) != genomeFeature:
                    print('\n------ERROR-------')
                    print("The earlier processed genome or feature set belongs to genomeFeature: %s" % (aval))
                    print("Your current genomeFeature input: %s" % (str(genomeFeature)))
                    print("Please input either correct genomeFeature value or re-run the analysis with '-featurefile' or '-genomeFile and -gffFile'")
                    print("\nSystem will exit now")
                    sys.exit(0)
                else:
                    pass
                
            if akey =="frags":
                fragList =  i.strip('\n').split('=')[1].split(',')
                #print('fragList from fragMem:', fragList)
        fragMem.close()
        #fragList = [file for file in os.listdir() if re.search(r'.*\.[0-9].*\.fa', file)] # deprecated v1.03
        print ('The fragments: %s' % (fragList))
        
    FragEnd = time.time()
    print ('\n\nThe fragmentation time is %s\n\n' % (round(FragEnd-FragStart,2)))
    fh_run.write('Fragmentation time is : %s\n' % (round(FragEnd-FragStart,2)))
    
    ##########################################################################################################
    ### Prepare list of miRNAs from either local file if specifed or miRNA table
    con = ConnectToDB(dataserver, 0)
    miRs,mirTable,miRFile = miRinput(con)

    ############################################################################################################
    ## TARGET PREDICTION ########################################################################################
    TFStart = time.time()
    
    ## Remove results from previous run
    if indexStep =='Y':
        shutil.rmtree('./index', ignore_errors=True)
        os.mkdir('./index')
        
    if tarPred =='Y' and tarScore =='Y':
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
        
    elif tarPred == 'N' and tarScore == 'Y':
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
    print ('\n\nThe Target prediction time is %s\n\n' % (round(TFEnd-TFStart,2)))
    fh_run.write('Target prediction time is : %s\n' % (round(TFEnd-TFStart,2)))
    
    ##############################################################################################################
    ## PARE PROCESS AND MAP ######################################################################################
    PAREStart = time.time()
    
    if getPARELib == 'Y': ## Get the files from run master
        TagAbundanceFile(con,PAREdb,libs)
        
    
    if tag2FASTAstep == 'Y':
        shutil.rmtree('./PARE',ignore_errors=True)
        os.mkdir('./PARE')
        PP(tag2FASTA2,libs) ##Convert tag count file to FASTA before mapping
    
    indexFls = [file for file in os.listdir('./index') if file.endswith ('index.1.bt2')]
    print ('These are index files: ',indexFls)
    
    if map2DDstep == 'Y':
        shutil.rmtree('./dd_map',ignore_errors=True)
        os.mkdir('./dd_map')
        global templib
        for templib in libs:
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
    
    ############################################################################################################
    # Pred MODULES ##############################################################################################
    
    ###Timer AK - will be written to run log file
    PredStart = time.time()
    
    if predStep == 'Y':
        ####Inputs
        fastaOut = 'genomic_seq.fa' ## Calculate genome statistics --- GET RID OF EITHER ONE OF VARIABLE
        predTargets = 'predicted/All.targs.parsed.csv'
        shutil.rmtree('./PAGe',ignore_errors=True)
        os.mkdir('./PAGe')
        shutil.rmtree('./output', ignore_errors=True) ## AK Added 
        os.mkdir('./output') ## AK added
        validatedTargetsFilename = './output/MPPPValidated.csv' ##Final results
    
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
    
        baseCountsFile = readFile('baseCounts.mem')
        baseCounts = int(baseCountsFile[0])
        baseCountsOffTagLen = int(baseCountsFile[1])
    
        for tagCountFilename in libs:
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
            categoryList = writePAGeFile(PAGeIndexDict, genomeFeature,
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
            if(validatedTargets):
                writeValidatedTargetsFile(header, validatedTargets,
                    validatedTargetsFilename)
    ##Script timer
    PredEnd = time.time()
    print ('\n\nIndexing and Prediction run time is %s\n\n' % (round(PredEnd-PredStart,2)))
    fh_run.write('Indexing and Prediction run time is : %s\n' % (round(PredEnd-PredStart,2)))
    fh_run.write('Script run time is : %s\n' % (round(PredEnd-FragStart,2)))
    fh_run.close()
    

    ###################################################################################################
    ##### REV COORDS ##################################################################################
    if revMapstep == 'Y':
        
        print('\n***Entering RevMapCoord- parallel***\n')
        global coord_dict_wat, coord_dict_crick
        coord_dict_wat = {} ## Dictionary of genes at watson strand, +
        coord_dict_crick = {}###Dictionary of genes at crick strand, - strand
        #
        ####------------TEST PURPOSE ONLY---------------READ FROM FILE----####
        ###coords_file = open('./gene_coords', 'r')
        ###for line in coords_file:
        ###    i = line.split(',')
        ###    #print (i)
        ####--------------------------------------------------------------####
        ##
        ##if Local =='Y':
        ##    ####For Local analysis
        ##    dens_start = time.time()
        ##    print('**Reading Density file**')
        ##    fh_in2 = open(norm_density_file, 'r')
        ##    global denfile
        ##    denfile = fh_in2.read().split('>')[1:] ###read once use repeatdly, first is always empty when split by '>'
        ##    print ('Time took to read density file is %s' % (round(time.time()-dens_start,2)))
        ##    print ('**Density file read in sub-main function**\n')
        ##    global den_dict
        ##    den_dict = {}
        ##    for ent in denfile:
        ##        ent_splt = ent.split('\n')
        ##        den_dict[ent_splt[0]] = ent_splt[2:-1]
        ##    print('**Dictionary made for Density file**\n')
        ##    #####################
            
        if Local == 'N':
            global nproc
            nproc ='1' ## Need better handling
            for i in coords:###gene_coords is a list in script, also written out as file of same name
                #strand = i.split(',')[1] for file###TEST if reading from file
                strand = i[1]
                if strand == 'c':### if entry in reverse strand
                    atuple = (i[0],i[4],i[5])
                    coord_dict_crick[i[2]] = atuple###Gene name as key and chr_id,strand end and gene type as value
                else:
                    atuple = (i[0],i[3],i[5])
                    coord_dict_wat[i[2]] = atuple##Gene name as key and chr_id,strand end and gene type as value
            print ('**Strand dictionary made**')
        
        #####Read the scoring input extend file and change coordinates
        ResFls = [file for file in os.listdir('./output') if file.endswith ('_validated')]
        for afile in ResFls:
            fh_in = open('./output/%s' % afile, 'r')### PARE VALIDATED results
            header = fh_in.readline() ##Waste header
            ScoInpExt = [] ##list of final results or parlist
            for res in fh_in:
                res_strp = res.strip('\n')
                ent =res_strp.split(',')
                ScoInpExt.append(ent)
            print ('**List from ScoInpExt ready to feed to RevMapCoord function**')
            
            #Write results to file
            revmapRes = './output/%s_revmapped.csv' % (afile)
            fh_out = open(revmapRes, 'w')
            fh_out.write('%s,chr,strand,BindSite,CleaveSite,SmallWin,LargeWin,Ratio,Func,phasing,inverted\n' % (header.strip('\n')))
            
            ## # TEST- SINGLE CORE - For troubleshooting ----##
            ValidTarGeno = []
            for i in ScoInpExt:
               z = RevMapCoord(i)
               ValidTarGeno.append(z)
            
            ## PARALLEL MODE - Uncommnet test above for normal use  
            # print ('**Reverse mapping initiated**\n\n')
            # ValidTarGeno = PPResults(RevMapCoord, ScoInpExt) ## Results are in form of list
            # print ('**Reverse mapping complete for:%s\n\n\n' % (afile))

            for i in ValidTarGeno: ## Write Results from list to file
                if i != 'E13-3-13':##Filter for error 13-13-13 that is small window abundance = 0 and ratio calculation error
                    fh_out.write('%s\n' % (i))
                else:
                    print (i)
            fh_in.close()
            #fh_in2.close()
            fh_out.close()
            
    ####6. UPDATE TABLE ###########################################################################################################
    if uniqRes == 'Y' and tableUpdate == 'Y':
        uniqRevmapped = resultUniq()
        con3 = ConnectToDB(destserver,1)###The second input '1' is for future use in case of local_infile upload to update table
        TableUpload(con3,uniqRevmapped,miRTable,miRFile)

    elif uniqRes == 'N' and tableUpdate == 'Y':
        uniqRevmapped = './output/ALL_revmapped_Uniq.csv'
        con3 = ConnectToDB(destserver,1)###The second input '1' is for future use in case of local_infile upload to update table
        TableUpload(con3,uniqRevmapped,miRTable,miRFile)
    
    elif uniqRes == 'Y' and tableUpdate == 'N':
        uniqRevmapped = resultUniq()

    else:
        print ('**Table Update not selected so results not updated to mySQL table**')
    ###############################################################################################################################

#### RUN ##########################################

if __name__ == '__main__':
    
    if nproc == 'Y':
        nproc = int(multiprocessing.cpu_count()*0.80)
    else:
        nproc == int(nproc)
    start = time.time()
    main()
    end = time.time()
    print ('Complete run time is %s' % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting..\n')
    sys.exit()

###################################################

### LOGS ################################################

##MPPP v01 -> v02
## Target prediction if OFF the main target prediction file is hard assigned
## dd2Trans optimized for speed - map files contain unalinged reads
## Target prediction optimized for senstivity

###MPPPv03 -> MPPV03 Int
##FASTA seq fetching from dataserver
## Rev map coord file generation Added
## target prediction schema mannualy verifed - Jan 19 - OK

###MPPPv03Int -> v04Int
#### Fixed critical bug in rebel mode scoring - IF statement in OR and AND statement
##### Changed FileCombine() and fh_in/fh_out filenames changed in tarParse3
##### Only first few miRNA mappings were being scored - FIXED - There were 'break' in both rebel and normal mode

###MPPPv04 -> MPPPv05
#### Optimized intergenic category calculation
#### Added miRNA input miRNArev complement file generated if miRNA file is input
#### Run log added that captures time and stats to a single file

####v05 -> v06:
###Core balancing - Automatic fragment number calculation
#### PP module used by bowtie based processing had nproc value modified to avoid overflow of server - Calculated by considering threads value

###V06 -> v07
####Bowtie - Modified the bowtie scoring for refgap and readgap modified to allow  1Gap+2MM or 5MM - Normal2 mode added for both exhaustive and heuristic
####Signal Category added in out put by Reza
#### Pvalue calculation modifed to have cumulative number of targets only for same score

###v07 -> v08
## adjusted p-value added to results
## Predictions with target score in decimel 0.5.1.5,2.5 etc now included in 0-1,1-2,2-3 section for p-value calculations
## This is the final corrected version with last category system, category 1 not calculated problems fixed


###v08 -> v09
## Reza added 12th position to match for validated targets. The schema was changed such a way that we add remove cleave indexes

###v09 -> v095
### Missing some steps between here...
### Corrected mistake where corrected p-values were computed on a per core level rather than a full compiled list of p-values.
### Changed pValue calculation to only use the number of trials of X.0 values with X.0. X.5 values remain full count of X.5 and X.0
### All the revmapped entries from different libraries need to be uniqed before updation i.e miRNA name, targetname and Cleave site [Remove header, combine files, Uniq, TableUpdate (Comment out header skip line)]
### The uniqed set should be used for updation to Online resource - Till than uniq manually and upload manually

###v095 -> v099
## Fixed the co-ordinates issue for scaffolds that is duplicated intergenic entries were bing generated and last gene on scaffold was n ot being included - Fixed on Scaffold- need testing on assembled genome
## Multiple libraries fails due to no threads available fixed
## Integration of all changes done with reza

###V099 -> v01
## Before uniq and upload sorted on corrected pvalue instead of p-value
## Uniq was giving error because binding site column was selected instead if cleave site - Fixed
## miRNA input with a 'T' was predicting wrong targets because rev complemting miRNA expected 'U' - Fixed
### Fixed TableUpload required miRNA input file name, a variable added to miRNA input module that is fed to table update
### Better GFF parser, The result header is screwed

## v1.0 -> v1.01 - sPARTA published !!!!!!!!!!!!!!!!
##Added functionality to fetch cleaner set of tags from the PARE DB

## v1.01 -> 1.02
## Reverse fix in miRferno, and implemented frag memory files, reza fixed one of the bug
## Header for output fixed
## Bowtie 2.4 compatibility fixed

## v1.02 -> v1.03
## In tarParse - Added functionality to catch gap in miRNA - Critical modificiation
## Needs further testing to see that when there a gap in miRNA - the target is regenrated accurately
## that means whther misInfo accounts for gaps in miRNA (deletion) or do we need to add number of gaps to correctly replace nucletide in targets
## This can be achieved by allowing more gaps and than studying the predicted results

## v1.03 -> v1.07 [major update]
## Added new PAGe index structure from sPARTA
## Added new capability to find targets for longer seqeunces i.e 24nt - Taken from sPARTA external


## v1.07 -> 1.08 [unstable]
## 'R' rebel scoring changes to 'S' seedless
## 'temp' folder renames to 'output'
## A switch added to generate uniq output file
## readGap value changed in exhaustive mode to allow one bulge in miRNA
## Added miRNA with gaps symbol, if any to output file
## Checked modifications from v1.03 for stable implementation - Pain in ass it was -OK
## Added miRNA-size based cleavage site search for miRNAs between 18-25nt length
## tarPARSE3 is updated in this version and sPARTA external need to have this module updated from here.

## v1.08 -> v1.09 [Critical Update]
## miRferno smallRNA mapping paramenetrs changes to improve sensitity especially to easily accomodate bulges and gaps
## Added functionality in miRferno to capture multiple gaps and or bulges in miRNA/targets
## Added functionality to allow consequitive bulges in target
## Added functionalty to allow bulge/gap in both miR and target
## Added a new validity to check if miR-traget secondary structure is biolgically possible

## v1.09 -> v1.10
## Bowtie parameters optimised for heuristic to improve speed with sensitive target prediction introduces in earlier version
## More then 5 mismatches and 6 edits (mismatches and wobbles) are not allowed
## Added generation fo uniq results in local mode

## v1.10 -> v1.11
## Fixed PARE tag fetching with hits less then 20. The filter was not working before, also maize had different tag abundance system that is Abun =0 if it
## hits multiple places.

## v1.11 -> v1.12
## Changed modules names to avoid super confusion about modules that extract features/coords and getFASTA froms erver or local file

## v1.12 -> v1.13
## Added loop for table upload and fileuniq to automate picking of uniq file

## v1.13 -> v1.14
## Fixed window length, from 5nt to 10nt. Mysql query were collecting abundances for just 5nt (+/-2)
## Added annotations from test tables - phased_loci and inverted repeats, controlled by a switch, extraAnno == Y

## v1.14 -> v1.15
## fixed writing results when extraAnno == 'N', result file is generated without phasing and inverted annotations

## v1.15 -> v1.16
## Major fix fro scaffolded genomes, if scaffold does not have any gene it will be fetched and not missed

### MODIFICATIONS ########################################
## getFASTALocal : 1) Split the coords file into required pieces as in fragFASTA 2) Read the genome in chunks i.e. required chromosomes 3) Make fragments
## nthread  Balancing - Threads are assigned on basis of load or number
## Bonferoni correction in Reza predictor
## Target gene plot with cleavage site market in Reza function
## In ddmap step - instead of intializing a globalTemp -Two arguments can directly be passed: http://stackoverflow.com/questions/13490346/python-gevent-pass-multiple-arguments-to-function-called-via-pool-map
## Reoptimize Heuristic and Exhaustive mode - Optimized/DONE
## Check if percentiles include outliers for intergenic category calculations
## getFASTAServer and fragFASTA could be merged

## v1.14 - Added two new columns to revmapped file - phasing and inverted - Modify table upload
