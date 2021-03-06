#!/usr/local/bin/python3


##This script makes use of Cleaveland3 pipeline to predict PARE validated targets. This version of script is written for biologists as some significant changes for them are made like below:
##Originally final results were filtered on basis of miR seq+targetname+cleavage site but for biologist miRseq was replaced by miR name
##You may want to run CL2 2nd script which creates an excel sheet of PARE reads in small and large windows because there could be non-specificity in targets.
##MiRNAs from same family can have same sequence and thus targets so you have to filter by using CL2 Jixian 2nd script
#USAGE: ./NAME_OF_SCRIPT.py
##This script was written by Atul kakrana: kakrana@udel.edu

###NEED TOOL FOLDER IN YOUR HOME

##


import mysql.connector as sql
import re
import sys
import os
import subprocess
import csv 
import glob
import time
import shutil
import itertools as it
import operator
import multiprocessing
from multiprocessing import Process, Queue, Pool
from operator import itemgetter
import datetime

################################################################ @@ CONFIG @@ ######################################################
##1. Input cDNA/transciptome file name
trascriptome_file = 'genomic_seq.fa'##The FASTA header will be cleaned automatically
##2. A miRNA File specified bu user - If file input used than table input will be turned off automatically
#miRNA_file = 'miRinput.fa'##The header will be cleaned automatically
#miRCode = 'osa'
###Bowtie processors
nthread = 24

##3. PARE database names
PAREdb = 'RICE_drt_PARE'

##3B. Genome database name
GenomeDB = 'RICE_MSU7_genome'

##4. Cleaveland folder
##Change CleaveLand address if SOY_Gmax101_genomerequired in 'CleavelandAnalysis' module

##5. Index filename/address
##Please change index filename in 'mapdd2trans' module. Index should be made on cDNA file to be used in this study and must have clean header, use 'clean_fasta_header_v1' for that.

##6. Bowtie folder
##Change folder address if required in 'mapdd2trans' module

##7A. You may want to change length of PARE tags in 'CleavelandAnalysis', tag length starts from 18 and goes to some 28, recommended 20
tag_len = 20

##7B. Pvalue for cleavland analysis - in case whole db is used as one PARE Db relaxed p-value should be used else 0.02 is good
pvalue = '0.5'

##8. TargetFinder score cutoff (max allowed is 7)
cutoff = str(7)

##9. Server with PARE data
dataserver = 'raichu.dbi.udel.edu'####Value should be in ''

##10. Server to upload data a.k.a destination server with table to hold results of PARE validation
destserver = 'raichu.dbi.udel.edu'#####Value should be in ''

##11. PARE validation done library wise or for whole DB. a value of zero means DB wise and 1 means library wise
LibPref = 1

##12. Number of processors to use for Parallel processing instances, 'Y' will use default setting of 75% of cores, Number of cores can be specifed in integer format too for ex. '12' will use 12 cores
nproc = 'Y'

##13. Genome level (value = 0) or Gene level (value = 1). Though default but Genome level can take long time.
#level = '0'

##14 A. Table to use for miRNA input - miRBASE or Kevin's - Must be in mir_master - Not required if Local miRNA file used
##14 B. miRNA code as in mirBASE - The script will fetch miRNAs from mirBASE table in mir_master
miRTable = 'mirBASE'
ver = '20'
miRCode = 'osa'

##15.  Drop the Results table before uploading results - Options - 'Y' and 'N' -- Be careful selecting 'Y' will wipe out the table and you will lose all the data
TableUpdate = 'N'
TableWipe = 'N'
##TableName = 'kakrana_test'

##16. Make index or use existing (index must be in 'index' folder withih same directory)
make_index = 'N'

##15. Local analysis i.e. no reverse coord mapping to be done for web_viewer
global Local
Local = 'N'

##16. Perform target prediction - 'Y' for Yes of 'N' for no
global TargetPred
TargetPred = 'N'

##17. Get PARE data or folder already exists
global GetPARE
GetPARE = 'N'
## 18. Inter or genic
genomeFeature = 1 ## 0 for gene and 1 for inter; 2 for both

########################################################### END OF CONFIG ############################################################


#########AA###################################TT###############################################UU########################################################LL##################

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
    
##Module to connect to DB
def ConnectToDB(server, infile):
    
    ##infile values are '0' when you dont want to pulaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    
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

## This module interacts with genome DB to take gene coords, calculate intergenic co-ordinates as well as generate a FASTA file for input
def GetCoords(con,db):
    cur= con.cursor()
    
    ###Filteration of RNA types not done at this stage becasue we need all the genes to get intergenic regions
    ###after filteration few genes will be missing and counted as intergenic region
    ###########################################################TEST LINE#################################################################
    #cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))###extra % escapes the % in query
    #####################################################################################################################################
    ##If we remove miRNA here than only gene entry will be removed but its sequence will be covered into intergenic
    if genomeFeature == 1: ##Intergenic
        cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type NOT LIKE 'mirna' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))###extra % escapes the % in query
    elif genomeFeature == 0: ##Protein Coding
        cur.execute("SELECT chr_id,strand,gene,start,end,type FROM %s.gene_master WHERE type LIKE 'protein_coding' AND gene NOT LIKE 'BLAST%%' ORDER BY chr_id,strand,start" % (db))###extra % escapes the % in query
    genome_info = cur.fetchall()###List with information of all the genes in genome
    ##Check if list is empty
    if not genome_info:
        print ('^^^Gene Coords query returned with an empty list..Exiting^^^')
        sys.exit()
    ####(1, 'c','AT1G01020', 5928, 8737, protein_coding)
    #print (genome_info)
    
    ##Find length of chromosomes to calculate intergenics
    cur.execute('SELECT chr_id, length FROM %s.chromosome_master' % (db))
    chromo_len = cur.fetchall()
    chromo_dict = dict(chromo_len)###Made a dict so that chromosome numer could be searched to get chromosome length
    print ('These are the chromosomes: %s and their length' % (chromo_dict))
    
    genome_info_inter = genome_info ###This list will also hold intergenics
    
    ####GET INTERGENIC REGIONS AND APPEND TO THE LIST WITH GENE COORDS
    alist = []###list maintained to check if first gene on chromosome and strand shows up than intergenic is just the start of gene
    #for gene1, gene2 from izip(*[iter(genome_info)]*2):
    #for i,j in pairwise (genome_info):
    #for gene1, gene2 in it.izip(genome_info[1:], genome_info):
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
                inter_end1 = gene1[3]-1###1 nt before start of Gene1 a.k.a the first gene on chromosome in this case
                ##As two genes are read together, the upstream intergenic region gor gene2 must be calculated in same step
                inter_start2 = gene1[4]+1##From end of first gene of chromosome
                inter_end2 = gene2[3]-1###Till start of second gene
                
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
                    inter_start = gene1[4]+1###End of Gene 1
                    inter_end = gene2[3]-1 ###1 nt before start of gene 2
                    if gene2[1] == 'w': ##Positive strand
                        inter_name = ('%s_up' % (gene2[2]))
                    else:## reverse strand
                        inter_name = ('%s_up' % (gene1[2]))
                    genome_info_inter.append((gene2[0],gene2[1],inter_name,inter_start,inter_end,gene_type))
                
                else: ###That means gene1 is at end of one chromosome and gene 2 is begining of chromosome so we have to extract intergenic at end of one chromosome
                    inter_start = gene1[4]+1###End of gene1
                    inter_end = chromo_dict[gene1[0]]###End of chromosome searched using chromosome id of gene1 from chromosome dictionary
                    if gene1[1] == 'w':##Positive strand end of chromosme
                        inter_name = ('%s_down' % (gene1[2]))
                    else: ##Negative strand first intergenic of chromosme
                        inter_name = ('%s_up' % (gene1[2]))                    
                    genome_info_inter.append((gene1[0],gene1[1],inter_name,inter_start,inter_end,gene_type))##Chr_id, strand
    
    
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

###This module gets FASTA sequence of supplied coords
def GetFASTA(con,db,gene_coords):
        ####Get the FASTA chromosome and strand wise
    ####NOTE the coords list must be sorted on chr_id and strand as FASTA seq of chr_id and strand is loaded when first instance  encountered. so to elimanate repeated loading sort them
    cur= con.cursor()
    genome_file = './genomic_seq.fa'
    fh_out = open(genome_file, 'w')
    #test_out = open('test_out', 'w')
    chromo_mem = []##used as memory to keep track of Chr+strand in memory, if its in here chromosome will not be read, untill its not here and appended
    for i in gene_coords: ###For every entry either gene or intergenic, use gene_info for just genes
        #print (i)
        gene = i[2]
        chr_id = i[0]
        strand = i[1]
        start = i[3]-1###Adjusted to pick from begining because when spliced it does not include starting position so starting from one poistion back
        end = i[4]
        #test_out.write('%s,%s,%s,%s,%s\n' % (chr_id,strand,gene,start,end))
        cur = con.cursor()
        if tuple(i[0:2]) not in chromo_mem: ### Add chromosome and strand to list and fetch in memory, will be executed first time a unique chr_id and strand found
            ###as list is sorted on chr_is and strand
            chromo_mem.append(tuple(i[0:2]))###Than append
            #print(chromo_mem)
            
            print ("Reading chromosome:%s and strand: '%s' into memory to splice genes" % (i[0],i[1]) )
            cur.execute("SELECT chromosome FROM %s.chrom_sequence where chr_id = %s and strand = '%s'" % (db, chr_id, strand))###Feed into the memory
            chromo = cur.fetchall()
            #print (chromo[0][0])###To get the string out of brackets            
            gene_seq = chromo[0][0][start:end]#Get first gene on that strand
            if strand == 'c':
                gene_seq_rev = gene_seq[::-1]
                fh_out.write('>%s\n%s\n' % (gene,gene_seq_rev))
                
            else:
                fh_out.write('>%s\n%s\n' % (gene,gene_seq))
                
            #print (gene_seq)           
                
        else:###If chr_id is in chromo i.e. its loaded in memory just splice your sequence
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

##Module to fetch miR name and sequence from master.mirna_filter_results table on raichu - this table was populated  by Kevin
def miRinput(con):
    miRs = [] ## List to hold miRNAs
    #fh = open(miRNA_file)
    try:## If user specified some local miRNA file in config
        fh = open(miRNA_file)##Test if file exists or not
        
        miRNA_file_clean = CleanHeader(miRNA_file)## Not required when info is downloaded from server in case of miRPage script
        fh_miRNA = open(miRNA_file_clean, 'r')
        mir_base = fh_miRNA.read()
        mir_blocks= mir_base.split('>')
        for i in mir_blocks[1:]:
            #print (i)
            block = i.strip('\n')##Remove newline from end of entry
            ent = block.split('\n')##Use the newline between header and sequence to split
            #print (ent)
            print ('%s,%s,%s,%s' % (ent[0],'None','None',ent[1]))
            miRs.append((ent[0],'None','None',ent[1]))## None at position 2 and 3 because so that name and seq at same index when data collected from table
        fh_miRNA.close()
        mirTable = 'None'##To be filled in final table
        print ('***miRNAs fetched from the specified file***\n')
        
    except NameError:## NO Local file specified by User - Fetch  from table
        print ('\nFetching miRs from the miRNA table')
        fh_out = open('miRinput.fa', 'w')## File with miRNAs for troupleshooting
        cur = con.cursor()
    
    ######################TEST SECTION - COMMENT OUT WHEN NOT TESTING###########################
        #cur.execute("SELECT mirna_name, mirna_precursor_name, score, mature_sequence FROM master.mirna_filter_result where score > 1 group by mirna_name limit 24")
    
    #####################REAL WORLD - UNCOMMENT WHEN NOT-TESTING#############
        #cur.execute('SELECT mirna_name, mirna_precursor_name, score, mature_sequence FROM master.mirna_filter_result where score > 1 group by mirna_name')##For kevins table
        cur.execute("SELECT mir_name, accession, organism, mir_seq FROM mir_master.mirBASE WHERE mir_name like '%s%%' AND version = %s" % (miRCode,ver))
    #####################Uncomment when testing is done#########################################
    
        miRs = cur.fetchall()
        print ('Total number of unique miRs found in miR table: %s' % (len(miRs)))
        
        ##Extract miRname and miR sequence for PARE validation and output in FASTA format as read by cleaveland
        #fh_out = open('miRinput.fa', 'w')
        
        for miR in miRs:
            mirname = str(miR[0])##Converted to string from unicode, a unicode mirname example: u'meyers-miR544-5p'
            mirseq = str(miR[3]) ##Converted to string from unicode, a unicode mirname example: u'TGAAGATGAAGAAGATGAAGAAGA'
            #print (mirname, mirseq)
            fh_out.write('>%s\n%s\n' % (mirname,mirseq))
        fh_out.close()
        mirTable = miRTable
        print('****miR info cached****\n')
        
    for i in miRs:
        print (i)
        
    return miRs, miRTable##miRs holds the list of miRNA name and query where as miRtable holds flag -table name or local

###Module to predict potential targets using targetfinder --- Parallel version
###As miRs are downloaded on fly, target prediction needs to be done on fly too
def TargetFinder(miR):
    #print (miR)
    query = str(miR[0]) ##Converted to string from unicode, a unicode mirname example: u'meyers-miR544-5p'
    mirseq = str(miR[3]) ##Converted to string from unicode, a unicode mirname example: u'TGAAGATGAAGAAGATGAAGAAGA'
    print('\n\nPredicting targets of miR: %s with sequence: %s' % (query,mirseq))
    ###To write result for every miRNA to separate file
    filename = miR[0].split()[0]
#    print(filename)
    fh_out=open('./target_finder_results/%s_targets' % (filename), 'w')
    
    ###Must use clean header file for cDNA
    
    ###Proposed change use getstatusoutput (targetfinder) and take results which is 2nd item of returened tuple
    pipe = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/TargetFinder_1.6/targetfinder.pl", "-c", cutoff, "-s", mirseq , "-d", genome_file, "-q", query], stdout=subprocess.PIPE, universal_newlines=True)
    result = pipe.stdout.read()
    
    print (result)
    
    fh_out.write('%s' % (result))
    fh_out.close()

####Module to get PARE data
def TagAbundanceFile(con,db):
    ##Remove pre-existing library from previous run
    shutil.rmtree('./PARE',ignore_errors=True)
    os.mkdir('./PARE')
    
    ##Get the number of libraries and make a list of their lib_ids:
    # Set cur to be the cursor so we can execute a query
    cur = con.cursor()
    cur.execute('select distinct(lib_id) from %s.run_master' % (db))
    only_libs = cur.fetchall()
#    print (libs)
    print ('\nTotal number of PARE libraries found: %s\n' % (len(only_libs)))
    
####!!!!!!!!!!!!!!!!!!!!!!!!!!TEST SECTION - COMMENT OUT WHEN NOT TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ##USe only the first library for testing purpose
    #libs = libs[0:2]
####!!!!!!!!!!!!!!!!!!!!!!!!!!END OF TEST SECTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ###As per disucssion in early Feb-13 called by Kevin incluidng BLAKE, JIXIAN and ME, instead of individual libraries a single file with norm hits should be used
    ###User can use either single or multiple libraries
    
    all_libs = [(0000,)] ## An artificial list so thet we can use this in place of libs which is a real list created in case of individual libraries
    if LibPref == 0:
        
        print ('PARE validation will be done DB wise, as selected')
        fh_out = open('./PARE/%s_PARE_tags.fa' % (all_libs[0][0]), 'w')##Naming file with lib_ids name, 0000 means all the libs
        print ('Caching tag and count information from server for all PARE libs\n')
        cur.execute('select tag, norm from %s.run_master' % (db))
        lib_info = cur.fetchall()
        
        tag_num = 1 ### for tag naming that will be used in FASTA format
        for ent in lib_info:##All the entries of the library
            if len(ent[0]) == tag_len: ## Length of tag specified by the user              
                for count in range(ent[1]):##Number of times the tag_count file
                    fh_out.write('>%s\n%s\n' % (tag_num, ent[0]))
                    tag_num += 1
            else:
                #print ('Length is not 20nt')
                pass
        fh_out.close()
            
    else: ###if LibPref = 1, that means PARE VALIDATION done individually
        print ('PARE validation will be done Library wise, as selected\n')
        for lib in only_libs:##For all the libraries
    #        print (lib[0])
            print ('Caching tag and count information from server for PARE lib %s' % (lib[0]) )
            cur.execute('select tag, raw_value from %s.run_master where lib_id = %s' % (db,lib[0],))
            lib_info = cur.fetchall()
    #        print(lib_info[0])
            
            ##Write to file
            #print('Writing PARE Fasta file for lib: %s\n' % (lib [0]))
            fh_out = open('./PARE/%s_PARE_tags.fa' % (lib[0]), 'w')##Naming file with lib_ids name
            tag_num = 1 ### for tag naming that will be used in FASTA format
            for ent in lib_info:##All the entries of the library
                if len(ent[0]) == tag_len:##Length of tag specified by the user            
                    for count in range(ent[1]):##Number of times the tag_count file
                        fh_out.write('>%s\n%s\n' % (tag_num, ent[0]))
                        tag_num += 1
                else:
                    #print ('Length is not 20nt')
                    pass
        fh_out.close()
        
        print("**Caching Done - Let's Rock the server\n**")
        
        #for i in only_libs:
        #    all_libs.append(i)###prepare a master list will all different libs and a db level lib that contains data of all PARE TAGS, required for PARE abundance at small window and large window calculation
        #
    return only_libs, all_libs###

##Function that gets just the list names, required to run script in parts. Run only when TagAbundanceFile function is not used
def JustLibList(con,db):
    cur = con.cursor()
    cur.execute('select distinct(lib_id) from %s.run_master' % (db))
    only_libs = cur.fetchall()
#    print (libs)
    print ('\nTotal number of PARE libraries found: %s\n' % (len(only_libs)))
    
    ##All libs includes an additional combined library
    all_libs = [(0000,)]
    #for i in only_libs:
    #        all_libs.append(i)###prepare a master list will all different libs and a db level lib that contains data of all PARE TAGS, required for PARE abundance at small window and large window calculation
    
    return only_libs, all_libs###   

####Module to map degradome to transcriptome/geome - RUN AWAY SAFE
def mapdd2trans(libs):##Creates index on fly and map PARE tags to index
    ##remove pre-existing directories from previous run
    if make_index == 'Y':    
        shutil.rmtree('./index',ignore_errors=True)
        os.mkdir('./index')
        #######################################TEST##################################
        #retcode = 0
        #######################################Comment after testing and open retcode line below ##########
        ##Make index on fly for the downloaded genomic sequences
        print('**Creating index of cDNA/genomic sequences\n**')
        retcode = subprocess.call(["/data2/homes/kakrana/tools/bowtie-0.12.7/bowtie-build", genome_file, "./index/an_index"])
    
    else:
        print('*Using existing index files*\n')
        retcode = 0 ##Pass the signal that index exists and process should not halt in next step
    
    
    shutil.rmtree('./dd_map',ignore_errors=True)
    os.mkdir('./dd_map')
    
    maxhits = str(30) ### If more than this than remoive the read
    mismatch = str(1) ### mismatch value for bowtie mapping
    report_seq = str(10) ###number of sequences to report from bowtie mapping, if ubique mapping than only one loci reported else the specifed
    nproc2 = str(nthread)
    if retcode == 0: ###That means index was prepared without any exit status error than procees to PARE mapping to index
        for lib in libs:
            print ('\n**The library %s is being mapped to transcriptome index file**' % (lib[0]))
            dd_file = ('./PARE/%s_PARE_tags.fa' % (lib[0]))
            map_out = ('./dd_map/%s_map' % (lib[0]))
    #        fh_in = openopen('%s_PARE_tags.fa' % (lib[0]), 'r')
            retcode2 = subprocess.call(["/data2/homes/kakrana/tools/bowtie-0.12.7/bowtie", "-v", mismatch,"-m",maxhits,"--best", "--strata", "-k", report_seq, "-p",nproc2, "-f", "./index/an_index", dd_file, map_out])
            #pipe = subprocess.Popen(["/data2/homes/kakrana/tools/bowtie-0.12.7/bowtie", "-v", mismatch,"--best", "--strata", "-k", report_seq, "-f", "./index/an_index", dd_file, map_out])##@@@@@@CHANGE INDEX FILENAME HERE@@@@BOWTIE FOLDER CHANGE
            #time.sleep(20) ## Fixes issue 1 below
            if retcode2 == 0:## The bowtie mapping exit with status 0, all is well
                print('\nDegradome from PARE lib: %s mapped to cDNA/Trascript file' % (lib[0]))
                
            else:
                print ("There is some problem with mapping of PARE lib: %s to cDNA/genomic seq index" % (lib[0]))
                print ("Script exiting.......")
                sys.exit()
            
    else:
        print ("There seems to be an error in building index for the cDNA/genomeic sequences before PARE mapping to index\n")

####Module where cleavland scripts are run - Parallel version - need method for exit code - subprocess.call throws result on screen
###Original cleavland module broken into two different parts so that PARE tags from DB level need not to be used in this module
def CleavelandSumm(lib):##Both summarize density and Run final cleaveland script | filename is clean Header cDNA_file
  
    #for lib fed by PP module in main()
    map_file = ('./dd_map/%s_map' % (lib[0]))
    density_file = ('./dd_density/%s_density' % (lib[0]))
    fh_out = open(density_file, 'w')
    run_info_file = ('./dd_density/%s_density_run_info' % (lib[0]))
    print ('Summarizing degradome map for PARE lib:%s' % (lib[0]))
    ##Summarize Degradome
    pipe = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/CleaveLand3_map2dd.pl", "-d", map_file, "-t", genome_file, "-s", seq_len, "-f", "bowtie"], stdout=subprocess.PIPE, universal_newlines=True)##@@@@@@@@CHANGE CLEAVELAND FOLDER HERE
    results = pipe.stdout.read()
    fh_out.write('%s' % (results))
    time.sleep(20) ## after issue 1 fix, this is kept to let subprocess finish gracefully - Neet wait() or communicate() or something else here
    fh_out.close()##Added later - Comeback if encounter some problems

###Module where cleavland results are generated - Parallel version
def CleavelandResult(lib):
    ##For lib fed by PP module in main()
    cl_results = ('./cl_results/%s_cl_results' % (lib[0]))
    density_file = ('./dd_density/%s_density' % (lib[0]))
    fh_out2=open(cl_results, 'w')
    cl_results_info = ('./cl_results/%s_cl_results_info' % (lib[0]))
    print ('Running final cleaveland analysis step for PARE lib:%s\n' % (lib[0]))
    pipe2 = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/CleaveLand3_analysis.pl", "-d", density_file, "-t", "./target_finder_results/", "-p", pval], stdout=subprocess.PIPE, universal_newlines=True)##@@@@@@@@CHANGE CLEAVELAND FOLDER HERE
    results2 = pipe2.stdout.read()
    fh_out2.write('%s' % (results2))
    time.sleep(10)
    fh_out2.close()##Added later - Comeback if encounter some problems

####Module to parse cleaveland results - can be removed --
def CleavelandParse(libs):
       
    for lib in libs:
        cl_results = ('./cl_results/%s_cl_results' % (lib[0]))##filename you got from final cleaveland3 analysis
        #cl_results = ('./cl_results/%s_cl_results' % (lib))##Test line
        print('\nCleaveland result file for lib %s is being parsed' % (cl_results))
        #Open the Cleaveland 3 final results file using csv reader
        fh_in=open(cl_results, 'r')
        fh_in.readline()##Header line wasted
        csv_results=csv.reader(fh_in, delimiter='\t')   
    
        # Open the file to write parsed results
        cl_parsed = ('./cl_results/%s_cl_results_parsed.csv' % (lib[0]))
        #cl_parsed = ('./cl_results/%s_cl_results_parsed_test.csv' % (lib))##Test line
        fh_out=open(cl_parsed, 'w')
        #Write Header
        fh_out.write('miRNA name,target_gene,target_location,Cleavage_site,targetfinder_score, p-value,category\n')
    
     
        entry_count=0        
        for i in csv_results:
        #    print(i)    
            mirna=i[0]
            target=i[1]
            arange=i[3]
            cleavage = i[4]
            score=i[2]
            category = i[5]
            p_val=(i[6])
            fh_out.write('>%s,%s,%s,%s,%s,%s,%s\n' % (mirna,target,arange,cleavage,score,p_val[:8],category)) ### String can't be rounded therefore p-valu selected till 7th alphabet
            entry_count+=1           
        
        print('The total number of result entries were:', entry_count)
        print('The output file "%s" is generated for further analysis\n' % (cl_parsed))
    
        fh_in.close()
        fh_out.close() 
        
    return(libs)##Just like that, of no use

####Module to parse targetfinder results
def TargetFinderParse(TFfoldername): ###TF v2
    
    ##Remove pre-existing 'scoring' directory if any from previous run
    shutil.rmtree('./scoring', ignore_errors=True)
    
    ##This module combines all the result files from target finder and than parses the single result file to be used for scoring  
    os.mkdir('./scoring')
    ##Cat all files in target finder folder
    
    files_to_cat = TFfoldername+'/*' 
    fls = glob.glob(r'./target_finder_results/*')###NEED to include files to cat in this line instead of explicitly mentioning folder
    TF_res_comb = ('./scoring/TF_results_all')
    TF_out = open(TF_res_comb, 'w') 
    for x in fls:
        TFfile = open(x, 'r')
        data = TFfile.read()
        TFfile.close()
        TF_out.write(data)
    TF_out.close()
    
#    TF_out.write('%s' % (results))
#    TF_out.close()    
    
    ##Prepare clean target finder combined file for parsing
    
    TF_res_clean = './scoring/TF_results_all_cleaned'
    fh_out2 = open(TF_res_clean, 'w') ###Intermidiate file to write clean results
    fh_in2=open(TF_res_comb, 'r')##Can give errors due to file opening and closing
    No_res_re = re.compile('No results')
    for i in fh_in2:
        if re.search(No_res_re, i): ##checks for file saying : No results for ath420
            pass
        else:
            fh_out2.write(i)
    fh_in2.close()
    fh_out2.close()
    
    ##Open a files to be written in CSV format
    TF_parsed_file = './scoring/TF_res_parsed.csv'
    fh_output=open(TF_parsed_file, 'w') # File to compare entries with validated ones   
    ##Write header to outputfile
    fh_output.write('miRNA name,target_gene,target_location,mirna_seq,tar_seq,targetfinder_score\n')
    
    ##Open clean TF combined file for reading    
    fh_in=open(TF_res_clean, 'r')##Can give errors due to file opening and closing
    csv_in=csv.reader(fh_in, delimiter='\t')
    csv_table=[]

    ##Populate the list
    for row in csv_in:
        csv_table.append(row)
    
    #print(csv_table[0:12])        
    #print(csv_table[0:12])##+6 to get next block i.e 6 lines makes an entry
    
    entries=int(len(csv_table)/6)
    print('\nTotal number of entries in combined TargetFinder results file:',entries)


    #Extract Entry, miRNA, Target
    scr_re = re.compile('score=\d{1}\.{0,1}\d{0,1}')
    ran_re = re.compile('range=\d*-\d*')
    result_entries=0
    a_count=0###incremented to move the next block
    b_count=6
    for i in range(entries):
        
        score =''
        ent_loc=csv_table[a_count:b_count]##A miR and target entry complete
        info_block=ent_loc[0][0].split('|')##block containing names and other info
#        print('%s\n' % (info_block))
        
        mirna=info_block[0].split(',')[0].split('=')[1]
        target=info_block[0].split(',')[1].split('=')[1].strip(' ')
#        print(info_block[3])
        
        scr = re.search(scr_re, str(ent_loc[0]))## make string of list entry before regex search
        score = scr.group().split('=')[1]##scr.group()returns the matched value
        
        ran = re.search(ran_re, str(ent_loc[0]))## make string of list entry before regex search
        arange = ran.group().split('=')[1]  
        
        tar_seq=ent_loc[2][0].split()[2]
        mir_seq=ent_loc[4][0].split()[2]
        
    #    print(mirna,target,score,range, mir_seq, tar_seq)
    ##    break

        fh_output.write('>%s,%s,%s,%s,%s,%s\n' % (mirna,target,arange,mir_seq,tar_seq,score))
        a_count+=6
        b_count+=6
        result_entries+=1
        
    print('Number of entries in %s file: %s\n' % (TF_parsed_file, result_entries))
        
    # Close output file: parsed_out
    os.remove('./scoring/TF_results_all_cleaned')
    fh_output.close()
    fh_in.close()
    
    return TF_parsed_file

####Module that combines results from both targetfinder and cleavland analysis.
def MatchCLandTF(libs,TF_parsed_file):

    ##This module matches CL files (one from each PARE library) and TF file (single file) to generate files (one for each PARE library) with sequence information for scoring  
    
    ##Because there are multiple files from CL output, one for each PARE library from the database
    for lib in libs:
        
        ##Prepare output files and name according to library 
        matched_file = ('./scoring/%s_validated_out' % (lib[0]))## This file should be combined at the end of module to unify results from all library
        #matched_file = ('./scoring/%s_validated_out_test' % (lib))## test line
        fh_output=open(matched_file, 'w') 
        #csv_out=csv.writer(fh_output, delimiter=',')
        
        matched_supp = ('./scoring/%s_validated_out_supp' % (lib[0]))
        #matched_supp = ('./scoring/%s_validated_out_supp_test' % (lib))##test line
        fh_out_supplemental=open(matched_supp, 'w')
        csv_out_sup=csv.writer(fh_out_supplemental, delimiter=',')
            
        
        ### Matching starts
        validated_tuples = []
        
        ##Read cleaveland parsed file and write column 1-3 in validated tiples for matching with target finder file
        
        cl_parsed = ('./cl_results/%s_cl_results_parsed.csv' % (lib[0]))
        #cl_parsed = ('./cl_results/%s_cl_results_parsed_test.csv' % (lib))##Test line
        print ('\nThe file %s is being joined with TargetFinder Results' % (cl_parsed))
        with open(cl_parsed) as fh1:
            fh1.readline()
            csv_reader = csv.reader(fh1)
            for row in csv_reader:
                validated_tuples.append(tuple(row[0:3]))
        print('The total number of entries i.e tuples in %s are %s' % (cl_parsed, len(validated_tuples)))

        
        ##Read TF file and compare with the column 1-3 in 'Validates_tuples' list 
        TF_parsed_file = './scoring/TF_res_parsed.csv'
        TF_match_list = []####The function of list to hold matching target finder entries which will be used later to join TF and CL entries
        with open(TF_parsed_file) as fh2:
            fh2.readline()
            csv_reader = csv.reader(fh2)
            match_count=0
            for row in csv_reader:
                if tuple(row[0:3]) in validated_tuples:
                    #csv_out.writerow(row[0:6])###+ add cleavage information
                    TF_match_list.append(row[0:5])                    
                    #csv_out_sup.writerow(row[0:6])            
                    match_count+=1
                else:
                    pass
        print('Total number of entries joined with target finder final results:%s\n' % (match_count))
#        print('Please use the output files %s and supplemental file %s  for further analysis' % (matched_file,matched_supp ))

        with open(cl_parsed) as fh1:###Here we join TF and CL entries to obtain miR and tar sequence
            fh1.readline()
            csv_reader = csv.reader(fh1)
            for ent in csv_reader:
                key = tuple(ent[0:3])
                #print(ent[0],ent[1],ent[2],ent[3])
                for i in TF_match_list:
                    #print(ent[0],ent[1],ent[2],ent[3])
                    #print (i)
                    if tuple(i[0:3]) == key:
                        #miR = re.sub('-','',i[3])
                        #tar = re.sub('-','',i[4])
                        miR = i[3]
                        tar = i[4]
                        fh_output.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[1],ent[2],ent[3],ent[4],ent[5],miR,tar,[ent6]))
                    else:
                        pass
    return libs ##Just like that
    
####Module to combine the validated out results from all the libraries so that we can collectively remove redundant files in the next step
####This was separated from MATCHCLandTF module because it was producing empty combined file but after making separate module is working fine
def FileCombine():
    
    validated_comb = './scoring/validated_all_lib_comb'
    validated_out = open(validated_comb ,'w')
       
    fls = glob.glob(r'./scoring/*validated_out')
    print ('Combining all the matched aka validated files from CL3 results of different PARE libraries to single file before removing redundant entries\n')
    for x in fls:
        CLfile = open(x, 'r')
        data = CLfile.read()
        CLfile.close()
        validated_out.write(data)
    validated_out.close()
    
#    validated_out.write('%s' % (results))
#    validated_out.close()
#    
    return validated_comb

####Module that removes redundancy, separates results on basis of miRNA size, generates result for biologist and does other stuff 
def UniqScoringInp(validated_comb):
    
    ## 1. Remove redundant entries-with all the attributes same i.e DNA, chr#, cleavage site, entry name, tar, miRNA
    ## 2. Write different scoring inputs for final scoring steps i.e sco_inp_21, sco_inp_22 etc 
    
    ##Reading parsed output having redundant entries
    fh_in=open('./scoring/validated_all_lib_comb', 'r') 
    ##fh_in.readline()##Waste header line##There is no header
    #parsed_in=csv.reader(fh_in, delimiter=',')
    parsed_in = [line.strip('\n').split(',') for line in fh_in]
    ##One line is ['>osa-miR812a', 'LOC_Os07g14310_up', '33472-33493', '33484', '1.5', '0.440755', 'CAGGUUGCAAAUUGGCAGGCAG', 'GUCUAACGUUUGACCGUCCGUC']
    parsed_in.sort(key=itemgetter(5))##Sorted on p-value beofre removing redundant so that first entry stored is best among other same interactions across library
    #for i in lines:
    #    print (i)

    ## Output file
    sco_inp_ext = "./scoring/scoring_input_extend"
    fh_output=open(sco_inp_ext, 'w') 
    csv_out=csv.writer(fh_output, delimiter=',')
    
    fh_output2=open('./scoring/scoring_input', 'w') 
    #csv_out2=csv.writer(fh_output, delimiter=',')
    ##Method:3- To find unique entries from the 'Matched i.e validated entries from cleaveland output'
    
    added_keys=set()## A set to store first 3 elements from input file: miRNA-DNA, chr# and cleavage site and than use it to compare further entries in file
    
    parsed_out_count=0## To keep count of unique entries
    print('\nRemoving redundant entries aka more than one instance of an miRNA (PARE Validated)')
    for ent in parsed_in:
        print(ent[0],ent[1],ent[2])
        genename = ent[1].split('.')[0]##To avoid different variations of same gene to be counted as uniq
        lookup=tuple((ent[0],genename,ent[2]))##miR name + Target Gene+position of cleavage on gene, earlier miR sequence was used instead of miR name but that removes miR belonging to same family as they have same sequence
        if lookup not in added_keys:##That means such a tuple has not been recorded yet and is unique
            csv_out.writerow(ent)###Complete information written to scoring_input_extend file
            fh_output2.write('%s\t%s\t%s\t%s\t%s\n' % (ent[0],ent[5],ent[1],ent[3],ent[4]))##Information relevent to connect with Jixans pipeline script 2
            #csv_out.writerow([row[0]]+[row[3]]+[row[4]])##MiRNA name, miRNA seq and target seq
            parsed_out_count+=1
    #        outlist.append(row)
            added_keys.add(lookup)## Once a new entry is found it is recorded so as to compare and neglect further entries
        else:
            pass
            
    print('The number of unique entries found and will be used for scoring:', parsed_out_count )
    print('\n->>Two results file generated at this step "scoring_input" and "scoring_input_extend"<<-')
    
    #print(outlist)
    #print('The total number of entries after removing duplicates:', len(outlist))
    
    fh_in.close()
    fh_output.close()
    
    ###############  sub-module: Segregate miRNAs on the basis of size before actual scoring - eliminates the need to run count_mir_len
    
    
    print ('\nClassifying results on basis of miR size | Size is calculated after gap adjustment')
    print('Dont Panic! the results file generated in last step are retained\n')
    
    
    ##The output files corresponding to the length of miRNA
    fh1_out=open('./scoring/scoring_inp_21', 'w')
    csv1=csv.writer(fh1_out, delimiter=',')
    fh2_out=open('./scoring/scoring_inp_22', 'w')
    csv2=csv.writer(fh2_out, delimiter=',')
    fh3_out=open('./scoring/scoring_inp_23', 'w')
    csv3=csv.writer(fh3_out, delimiter=',')
    fh4_out=open('./scoring/scoring_inp_24', 'w')
    csv4=csv.writer(fh4_out, delimiter=',')
    fh5_out=open('./scoring/scoring_inp_other', 'w')
    csv5=csv.writer(fh4_out, delimiter=',')
    
    
    ##The input file to be read and filtered according to their length
    fh_in1=open('./scoring/scoring_input_extend', 'r')
    csv_in1=csv.reader(fh_in1, delimiter=',')
    
    ## Read the file line by line and do operation
    for i in csv_in1:
        a=i[6].count('-')##Count the bulges and gap    
        original_len=int(len(i[6])-a)
    #    print('Correct length:', original_len)
        if original_len == 21:
            csv1.writerow(i[0:3])
        elif original_len == 22:
            csv2.writerow(i[0:3])
        elif original_len == 23:
    #        print(i[2])
            csv3.writerow(i[0:3])
        elif original_len == 24:
            csv4.writerow(i[0:3])
        else:
            print('An miRNA %s: %s of length %s was found, it normal' % (i[0][1:],i[6],original_len))
            csv5.writerow(i[0:3])
            #print(i[0])
            pass
        
    print('\nPlease use the output file of required miRNA length from "scoring" folder for your analysis\n')
    
    ##Close all the files of this module
    fh1_out.close()
    fh2_out.close()
    fh3_out.close()
    fh4_out.close()
    fh_in1.close()
    
    return sco_inp_ext
    
####Parallelized version of RevMapCoord - introduced in version 2.04 - REGEX dropped and Dictionary introduced v2.05
def RevMapCoord(ent): ####DOES THIS MODULE WORKS IF USER USES MULTIPLE LIBRARIES - Parallelize
    ###create a dictionary from list of coords to be searched later
    ###gene_coords structure: 1, 'c','AT1G01020', 5928, 8737, protein_coding   
    #print (ent)
    gene_name = ent[1]
    cleave_site = int(ent[3])
    bind_site = ent[2].split('-')
    
    
    if Local == 'N':
        print ('**Reverse mapping of Co-ordinates will be performed - Web analysis**')
    ###Check whether gene is from reverse or posiive strand by memebr ship test on dictionary
        if gene_name in coord_dict_wat:
            print ('\nEntry: %s in positive strand: %s' % (ent[0:4],coord_dict_wat[gene_name]))
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
            print ('\nEntry: %s in reverse strand: %s' % (ent[0:4],coord_dict_crick[gene_name]))
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

    ###use cleavage site to find abundance
   
    small_win = 0###Small window abundance for this validated gene
    large_win = 0###Large window abundance for this validated gene
    print ('*Acquiring PARE abundances in small and large window............*')
    #found = den_dict[gene_name]
    #print ('Gene summary info\n%s' % (found))
    #if found: ### gene_name from sco inp extend results found in PARE summarization file and his record is fetched
    
    
    ##Fix the gene name on reverse strand and than you can use gene name as it is by removing _up/_down
    cur = con.cursor()##Connect to data server with table
    #gene_query = gene_name.split('_')[0]###For up and down sequnces
    gene_query = gene_name
    print ('Gene name:', gene_query,'Large win start:',int(new_cleave_site)-2,'Cleave site:',new_cleave_site,'Large win end:',int(new_cleave_site)+2)
    
    ##Making just single query for region of large window, can also make two queries for large and small window and get sum(norm_sum) in this case no need to calculate small win and large win abundances as below
    ##Single query because sometimes server is busy or slow
    if gtype == 'inter':##Gene name not used
        cur.execute("SELECT position,norm_sum FROM %s.tag_position where chr_id = %s AND strand = '%s' AND position between %s and %s" % (PAREdb,chr_id,strand,int(new_cleave_site)-2,int(new_cleave_site)+2))## No gene name
        window = cur.fetchall()
        #print(window,'\n')
    
    else:##Gene name used
        cur.execute("SELECT position,norm_sum FROM %s.tag_position where gene = '%s' AND chr_id = %s AND strand = '%s' AND position between %s and %s" % (PAREdb,gene_query,chr_id,strand,int(new_cleave_site)-2,int(new_cleave_site)+2))## Clear table before writing
        window = cur.fetchall()
        #print(window,'\n')
    
    for entry in window:
        #print(ent[0],ent[1])
        if int(entry[0]) <= int(new_cleave_site+1) and int(entry[0]) >= int(new_cleave_site-1):####Range of Small Window###Every 5' mapping on specific gene is read and if 5'end is mapped within small window than added
            #print ('Using %s for small window' % ent[0])
            small_win += (int(entry[1]))

        if int(entry[0]) <= int(new_cleave_site+5) and int(entry[0]) >= int(new_cleave_site-5):###Range of Large window###GET LARGE AND SMALL WINDOW IN SAME LOOP
            #print ('Using %s for large window' % ent[0])
            large_win += (int(entry[1]))    
    #
    ###Check for bug 13-3-13-In lib based validation, target is validated just raw reads mapped by cleaveland and sometimes there is no unique and normalized read mapped at corresponding site in '0_density' file so skip those
    if small_win != 0:
        #print( ('%s,%s,%s,%s,%d,%d,%f\n' % (new_cleave_site,ent[6],ent[7],ent[5],small_win,large_win,round(small_win/large_win,2))) )
        #rev_mapped_entry = ('%d,%d,%f\n' % (small_win,large_win,round(small_win/large_win,2))) 
        rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%f\n' % (ent[0][1:],gene_name,chr_id,strand,new_bind_site, new_cleave_site,ent[4],ent[6],ent[7],ent[5],ent[8],small_win,large_win,round(small_win/large_win,2)))
        #print (rev_mapped_entry)
    else:
        ##If the first result entry has 0 abundance for small window than no 'rev_mapped entry' variable and nothing to return, therefore added 'Rev_mapped_entry' with ratio set to 0
        #rev_mapped_entry = (('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0][1:],gene_name,new_bind_site, new_cleave_site,ent[4],ent[6],ent[7],ent[5],small_win,large_win, 0)))
        ##Or rev_mapped_entry just mentions the bug that caused blank entry which can be filtered later
        rev_mapped_entry = ('E13-3-13')
        pass

    #return sco_inp_extend_geno
    return rev_mapped_entry

####Module that uploads data to already made table in a DB
def TableUpload(con3,res_upload,mir_tab_name):###Con is connection and res_upload is file from last module which has genomic coords and need to be upload
    ##Upload scoring_input_extend to table - please make sure that columns are in file are in order of the columns in the table
    
    cur = con3.cursor()##Connect to destination server with table
    #res_file = 'scoring_input_extend_upload'
    
    try:## If user specified some local miRNA file in config
        fh = open(miRNA_file)##Test if file exists or not
        mir_tab_name = 'None'
    except NameError:
        pass
        
    
    
    if TableWipe is 'Y': ## If the Global setting is 'Y' than table will be dropped - Be careful   
        print ('\n**Clearing the table before updating.....**')
        cur.execute("TRUNCATE TABLE mir_master.mirPARE_v2")## Clear table before writing
        #con2.commit()
        print ('**\nTable cleared successfully, update in process**')
        
    ##Current implementation of mysql.connector does not support instant upload by local file - see ConnectToDB() module for implementation
    ##Original query - LOAD DATA LOCAL INFILE "./scoring_input_extend" INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ',';
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    ##So fill table on row by row basis
    ##Check for miRNA table variable if not than Local file used else - public/priv
    
    
    add_row = "INSERT INTO mir_master.mirPARE_v2(mir_name,target_name,chr_id,strand,binding_site,cleavage_site,target_score,mir_seq,target_seq,p_value,small_win,large_win,win_ratio,PAREdb,miRdb_type) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"
    
    #test_data =[["Test1", "AT2G28350.1", "1329-1349", "1340", "2", "ACCGUAUGUCCCUCGGUCCGU", "AGGAAUACAGGGAGCCAGGCA", "0.000314"],["ath-miR160b-5p", "Test", "1329-1349", "1340", "2", "Test", "Test", "0.000314"]]
    print (PAREdb,mir_tab_name )
    fh_in2 = open(res_upload, 'r')
    
    for row in fh_in2:
        res_entry = row.split(',')
        print(res_entry[0], res_entry[1], res_entry[2], res_entry[3], res_entry[4], res_entry[5], res_entry[6], res_entry[7],res_entry[8],res_entry[9],res_entry[10],res_entry[11],res_entry[12],PAREdb,mir_tab_name)
        res_upload = (res_entry[0], res_entry[1], res_entry[2], res_entry[3], res_entry[4], res_entry[5], res_entry[6], res_entry[7],res_entry[8],res_entry[9],res_entry[10],res_entry[11],res_entry[12],PAREdb,mir_tab_name)
        cur.execute(add_row,res_upload)        
        con3.commit()
        
    cur.close()

###Tiny PP module - Just PP no result fetching
def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()###time start
    npool = Pool(int(nproc))
    npool.map(module, alist)###Cannot give multiple parameters at this time

##Tiny PP module with results reported back as results
def PPResults(module,alist):
    npool = Pool(int(nproc))    
    res = npool.map_async(module, alist)###
    results = (res.get())###results in form of a list as returend by the module; either (parameter, result) or (result) format
    return results
    
##PP functions to queue jobs and get results - Check latest PP.py for better implementation 
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
            res = RevMapCoord(par)
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

###########################################################################- MAIN -####################################################################
def main(): ###Variables need to be initialized in the ----main-------
    
###1. CLEAN FASTA HEADERS - OBSOLETE -TEST ############################################################################
    global genome_file
    global con
    
    ### Begin ###############################
    con = ConnectToDB(dataserver,0)###The second input '0' is for future use in case of local_infile upload to update table
    miRs,mirTable = miRinput(con)###Fetch list of miRNAs from either local file if specifed or miRNA table
    ########################################################################################################################
    #
    ####2. GET DATA - PARE, miR[optional] and Sequence[optional] ##########################################################
    if Local=='Y': ##In local no gene coords required
        genome_file = CleanHeader(trascriptome_file)
    
    else:        
        gene_coords = GetCoords(con,GenomeDB)###Get selected genes and all intergenic coords - all coords all also fetched too
        genome_file = GetFASTA(con,GenomeDB,gene_coords)###Extracts FASTA of supplied coords
      
    ################## TargetFinder Parallel ################
    ##Timer
    runLog = 'runtime_%s' % datetime.datetime.now().strftime("%m_%d_%H_%M")
    fh_run = open(runLog, 'w')
    TFstart = time.time()###time start
    if TargetPred == 'Y':
        TFstart = time.time()###time start
        shutil.rmtree('./target_finder_results',ignore_errors=True) ##Delete previous ./targetfinder folder 
        os.mkdir('./target_finder_results') ##Open new directory for current run
        #for i in miRs:
        #    TargetFinder(i)
        
        PP(TargetFinder,miRs)###TargetFinder in parallel
    else:
        pass
    ########################################################
    TFend = time.time()
    TFtime = (round(TFend-TFstart,2))
    fh_run.write('Target Prediction: %s\n' % (TFtime))
    print ('Time taken for target prediction:',TFtime)
    #if GetPARE == 'Y': ## Fetch from MySQL DB  
    #    only_libs,all_libs=TagAbundanceFile(con,PAREdb)####PARE DATA - NEED TO BE PARALLELIZED
    #else: ##PARE data already exists
    #    only_libs,all_libs=JustLibList(con,PAREdb)##To run script in parts, do not run if TagAbundanceFile is being used
    
    #######################################################################################################################
    
    ######3.INDEX, MAPPING and CLEAVELAND ###################################################################################
    analysis_start = time.time()###time start
    #Combined PARE map will not be used to generate result if libwise result selected
    if LibPref == 0:###If DB level analysis than use this list for all later modules
        libs = all_libs
    else:####if library wise analysis use this module for all later modules
        libs = only_libs
        
    mapdd2trans(libs)##Prepare index and map PARE to Index
    ##########-----------CLeaveland analysis parallel----------###
    print ('\nStarting Cleaveland3 based analysis\n')
    shutil.rmtree('./dd_density',ignore_errors=True)##remove pre-existing directories from previous run
    os.mkdir('./dd_density')##Make new directories for this run
    PP(CleavelandSumm,libs)###all_libs has inflrmation of libs as well as DB level PARE abundance that will be summarized as well to be used in small and large window calculation by RevMapCoord module
    
    shutil.rmtree('./cl_results',ignore_errors=True)
    os.mkdir('./cl_results')
    PP(CleavelandResult,libs)
    analysis_end = time.time()
    analysis_time = (round(analysis_end-analysis_start,2))
    #########-------------------------------------------------####
    ######
    ##############################################################################################################################
    ######
    ######
    ##########4. PROCESSING RESULTS ###############################################################################################
    CleavelandParse(libs)
    TF_parsed_file = TargetFinderParse('./target_finder_results')
    libs = MatchCLandTF(libs,TF_parsed_file)
    validated_comb = FileCombine()
    sco_inp_ext = UniqScoringInp(validated_comb)
    ###########################################################################################################################
    ####
    ####
    #######5. REVERSE MAPPING PARALLEL ########################################################################################
    #
    global nproc
    nproc ='1' ## hack - Need better handling
    
    print('\n***Entering RevMapCoord- parallel***\n')
    global coord_dict_wat, coord_dict_crick
    coord_dict_wat = {} ## Dictionary of genes at watson strand, +
    coord_dict_crick = {}###Dictionary of genes at crick strand, - strand
    ###------------TEST PURPOSE ONLY---------------READ FROM FILE----####
    ##coords_file = open('./gene_coords', 'r')
    ##for line in coords_file:
    ##    i = line.split(',')
    ##    #print (i)
    ###--------------------------------------------------------------####
    for i in gene_coords:###gene_coords is a list in script, also written out as file of same name
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
    fh_in = open(sco_inp_ext, 'r')### PARE VALIDATED results
    ScoInpExt = [] ##list of final results or parlist
    for res in fh_in:
        res_strp = res.strip('\n')
        ent =res_strp.split(',')
        ScoInpExt.append(ent)
    print ('**List from ScoInpExt ready to feed to RevMapCoord function**')
    
    #Write results to file
    fh_out = open(sco_inp_extend_geno, 'w')
    fh_out.write('miRName,Target,Chr,Strand,BindSite,CleaveSite,TFScore,miRSeq,TarSeq,pVal,Category,SmallWin,LargeWin,WinRatio')
    #
    ###--- TEST- SINGLE CORE - For troubleshooting ----##
    ##for i in ScoInpExt:
    ##    RevMapCoord(i)
    
    #
    print ('**Reverse mapping initiated**\n\n')
    ValidTarGeno = PPResults(RevMapCoord, ScoInpExt)###Results are a in form of list
    for i in ValidTarGeno:##Write Results from list to file
        if i != 'E13-3-13':##Filter for error 13-13-13 that is small window abundance = 0 and ratio calculation error
            fh_out.write(i)
        else:
            print (i)
    fh_in.close()
    #fh_in2.close()
    fh_out.close()

    print('Analysis time: %s' % (analysis_time))
    fh_run.write('Analysis Time: %s\n' % (analysis_time))
    fh_run.close()
    ###############################################################################################################################
    
    
    ####6. UPDATE TABLE ###########################################################################################################
    if TableUpdate == 'Y':
        con3 = ConnectToDB(destserver,1)###The second input '1' is for future use in case of local_infile upload to update table
        TableUpload(con3,sco_inp_extend_geno,miRTable)
    else:
        print ('**Table Update not selected so results not updated to mySQL table**')
    ###############################################################################################################################

################################################- RUN -#####################################################


if __name__ == '__main__':
    
    start = time.time()###time start
    ##GLOBAL VARIABLES###        
    ######TEST PARAMETERS######
    genome_file = './genomic_seq.fa_new_head.fa'###Supplied by 'GetFASTA' function in REAL, Required by 'TargetFinder' as PP doesnt allow multiple arguments
    #norm_density_file = './dd_density/0_density' ## For 'RevMapCoord' to read norm PARE summarization file. The list all_libs if sorted than this will be all_libs[0] directly mentioned in 'RevMapCoord' function
    sco_inp_ext = './scoring/scoring_input_extend'## Supplied by 'UniqScoringInp' function in REAL
    #gene_coords = './gene_coords'### Supplied as list by 'gene_coords' function in REAL. Uncomment Test portion in 'RevMapCoord' to directly read from file
    #all_libs = [(0,),(574,)] ## Supplied by 'TagAbundanceFile' function in REAL
    
    
    only_libs = [(1153,),(1154,),(1155,),(1156,),] ## Supplied by 'TagAbundanceFile' module in REAl -------------- PARTH - name this as library number for two libraries [(123,),(456,)]
    
    
    sco_inp_extend_geno = './scoring/sco_inp_ext_geno'
    seq_len = str(tag_len)## For 'CleavelandSumm' function
    pval = str(pvalue)##For 'CleavelandResult' function


    ###Processors to use####
    if nproc == 'Y':###Use default 70% of processors
        nproc = int(multiprocessing.cpu_count()*0.85)
    else:##As mannually entered by the user
        nproc == int(nproc)
        ###############
    main()
    end = time.time()
    print ('The script run time is %s' % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting........\n')
    sys.exit()


###############################################- END -####################################################

'''
Created on Jun 2, 2012

@author: atul
'''

###Changelog since JAN-2013 - OBSOLETE - NEEDS UPDATION

##1. p-value was added to the scoring_input_extend file


###Bugfix since JAN-2013
##1. The matched CL and TF results were not getting combined properly in MatchCLandTF() module.
##FOr some reason code to combine files was generating empty result file it is therefore the code is moved to a new module to endure that it is being run
##2.

###ISSUES
##1. For some reason sometimes dd_map file generated in mapdd2trans(libs) module becomes unreadable to CleaveLand3_map2dd.pl script in CleavelandAnalysis() module
###As a result dd_density file remains empty and empty CL result file is generated - The error shows up in commandline but does not terminates the run.
##2. mysql.connector library for python3 does not yet support connection argument -local_infile = [0/1]. It is therefore result file cannot be directly uploaded to table
### For time being table update is done row by row in module TableUpload() module

###Proposed changes:
##1. module RevCoords prepares a list along with file and list is fed directlt to TableUpload
##2. cleavland module might need change like the mapdd2trans...or introduce process.wait() get return code and than process for cl3 2nd script

###V.2.08
##Long intergenic bug fixed - It was because coords were sorthed on just strand and chromosome
##3	w	Medtr3te099910_up	34793331	34916088	122757
##3	w	Medtr3te100000_up	34916299	34953292	36993
##So id intergenic is calculated ~ 122757 - WRONG - As there are many genes between them whose sequence and intergenic already included
##Fixed by adding 'start' as sort parameter

##From CL3v2.09 -> 2.10
#Bug in mirinput fixed that was skipping the file and downloading the miRNAs from mir table
#While making index the hardcoded genomic file name changed to script generated for the case where transcriptome file is specified by the user. The file is processed by cleaninh header
#and new name used instead of hardcoded

##Check OK - NO BUG

###From CLv2.10 -> 2.11 -> 2.12
#Added a check at gene coords function to see if the query return empty list - Prints error and sys exit
#Added table drop or Not before uploading results - switch at the top
#fixed miRNA db name if user does not use miRBASE or anyother table - it will be set to none
#REmoved '0' file system so no genome mapping to '0', no cleaveland analysis and no use in rev map coord
#Revmap coord doesn not use '0'file instead use SQL queries

###From v.2.12 -> 2.12A
#Large window size changed from 5nt to 11nt

###From v.2.12 -> 2.13 --> 3.02 (Regression from 3.01)
##Make index switch added
##Tragets sorted on p-value before removing redundant. This increases the number of good targets as one with low p-value is selcted
##Tested extensively at Trakan/test/CL3vs2/2.Cl2 --> Works perfectly - OK

##Proposed -> v3.03
##Filter empty gene entries while retrieving from Genome DB and also from CLeanFASTA module - Implement and Test both scenarios
##Added 'version' to fetch miRNAs from miRbase related to specied version
##Change mapping parameters to report more genes - Made - Tested on 2 miRNAs - No Error works OK
##Use PARE tags with length more than 20 - Test from 3.01
##Functionality to select 'intergenic' or 'protein'coding' from settings
###Max hits filter added  to dd_map
###Added category to the output

