#!/usr/local/bin/python3

## atulkakrana@gmail.com
## Revmaps internal miRferno results to genome co-ordinates


import sys,os,re,time,timeit,csv,glob,string,shutil
import subprocess, multiprocessing
from multiprocessing import Process, Queue, Pool
import operator
from operator import itemgetter
import mysql.connector as sql


#########################################################################
#### USER SETTINGS ######################################################

Local = 'N'                         ## 'N' means internal version and 'Y' means external version i.e. local files used instead of miRNA and fasta files from server
GenomeDB = 'MAIZE_AGPv2_genome'      ## DB to get fastaOut file
# PAREdb = 'MAIZE_pub_PARE'        ## Make sure that your Library is in the DB 
genomeFeature = 1                   ## 0 for gene and 1 for inter; 2 for both


nthread = 6                             ## Need automatic calculation like nProc
nproc = 'Y'                             ## Used by parallel processing
generateFasta = 'N'                     ## Functionality missing
dataserver = 'raichu.dbi.udel.edu' 

##########################################################################

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
    chromo_dict = dict(chromo_len)###Made a dict so that chromosome numer could be searched to get chromosome length
    #print ('These are the chromosomes: %s and their length' % (chromo_dict))
    
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
        gene_type = 'inter' ###set to integentic by default
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

def RevMapCoord(ent): ####DOES THIS MODULE WORKS IF USER USES MULTIPLE LIBRARIES - Parallelize
    
    ## Create a dictionary from list of coords to be searched later
    ## Gene_coords structure: 1, 'c','AT1G01020', 5928, 8737, protein_coding   
    
    print (ent)
    gene_name = ent[1]
    bind_site = ent[2].split('-')
    # cleave_site = int(ent[8])
    
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
            # new_cleave_site = (int(geno_start)-1)+int(cleave_site)###1 is reduced to give correct positions
            new_bind_site_start = (int(geno_start)-1)+int(bind_site[0])
            new_bind_site_end = (int(geno_start)-1)+int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
        else:
            print ('Entry: %s in reverse strand: %s' % (ent[0:4],coord_dict_crick[gene_name]))
            geno_end = coord_dict_crick[gene_name][1] ### Use for reverse mapping of negative genes
            #print('Looking for chr_id')
            chr_id = coord_dict_crick[gene_name][0]
            #print('chr_id found')
            strand = 'c' ## Available in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_crick[gene_name][2] ##Gene type
            # new_cleave_site = (int(geno_end)+1)-int(cleave_site)###1 is added to give correct positions
            new_bind_site_end = (int(geno_end)+1)-int(bind_site[0])###As the sequence was reversed before TF and CL, their binding start and end direction has also changed - Verified-OK
            new_bind_site_start = (int(geno_end)+1)-int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
   
    else: ## 'Y' i.e Local analysis
        print ('**No Reverse mapping of Co-ordinates will be performed - Local analysis**')

    rev_mapped_entry = ('%s,%s,%s,%s,%s' % (','.join(ent),str(chr_id),strand,new_bind_site_start,new_bind_site_end))

    ## Get window abundances and annotations ###############################################
    
    # small_win = 0 ##Small window abundance for this validated gene
    # large_win = 0 ##Large window abundance for this validated gene
    # print ('*Acquiring PARE abundances in small and large window............*')

    # if Local == 'N':
    #     ## Fix the gene name on reverse strand and than you can use gene name as it is by removing _up/_down
    #     cur = con.cursor() ## Connect to data server with table
    #     gene_query = gene_name
    #     print ('Gene name:', gene_query,'Large win start:',int(new_cleave_site)-5,'Cleave site:',new_cleave_site,'Large win end:',int(new_cleave_site)+5)
        
    #     ## Making just single query for region of large window, can also make two queries for large and small window and get sum(norm_sum) in this case no need to calculate small win and large win abundances as below
    #     ## Single query because sometimes server is busy or slow
    #     if gtype == 'inter':## Gene name not used
    #         cur.execute("SELECT position,norm_sum FROM %s.tag_position where chr_id = %s AND strand = '%s' AND position between %s and %s" % (PAREdb,chr_id,strand,int(new_cleave_site)-5,int(new_cleave_site)+5)) ## No gene name
    #         window = cur.fetchall()
    #         #print(window,'\n')
    #         genefunc = 'NA'
        
    #     else: ## Gene name used
    #         cur.execute("SELECT position,norm_sum FROM %s.tag_position where gene = '%s' AND chr_id = %s AND strand = '%s' AND position between %s and %s" % (PAREdb,gene_query,chr_id,strand,int(new_cleave_site)-5,int(new_cleave_site)+5)) ## Clear table before writing
    #         window = cur.fetchall()
            
    #         ## Gene Annotation
    #         cur.execute("SELECT gene,title FROM %s.gene_master where gene = '%s'" % (GenomeDB,gene_query))
    #         info = cur.fetchall()
    #         # print("Anno:",info)

    #         if info[0][1] != None:
    #             # print (info[0][1])
    #             genefunc = info[0][1].replace(",",";")
    #         elif info[0][1] == None:
    #             # print ('Title not found for: %s' % (gene_query))
    #             genefunc = 'NA'
    #             pass
    #         else: ## Not found
    #             # print ('Title not found for: %s' % (gene_query))
    #             genefunc = 'NA'
    #             pass


    #     if extraAnno == 'Y':
    #         ## Phased annotation - need phased_loci table in genome DB
    #         cur.execute("SELECT cluster_id,len FROM ASPARAGUS_UGA1_genome.phased_loci where chr_id = %s and (%s between start and end)" % (chr_id,int(new_cleave_site)))
    #         info2 = cur.fetchall()
    #         # print('Phasing Anno:',info2)
    #         if info2:
    #             phasing = '-'.join(str(i) for i in info2[0])
    #             # print('Phasing Anno:',phasing)
    #         # elif info2[0][1] == None:
    #         #     phasing = 'No'
    #         else:
    #             phasing = 'No'

    #         ## Inverted repeat annotation - Need chr_inverted table in genome DB
    #         cur.execute("SELECT * FROM ASPARAGUS_UGA1_genome.chr_inverted where chr_id = %s and ( (%s between start1 and end1)  or (%s between start2 and end2) )" % (chr_id,int(new_cleave_site),int(new_cleave_site)))
    #         info3 = cur.fetchall()
    #         # print("Inverted anno:",info3)

    #         if info3:
    #             inverted = 'IR'
    #         # elif info3[0][1] == None:
    #         #     inverted = 'No'
    #         else:
    #             inverted = 'No'

            
    #     ## Compute small window and large window using the data collected above
    #     for entry in window:
    #         if int(entry[0]) <= int(new_cleave_site+1) and int(entry[0]) >= int(new_cleave_site-1):## Range of Small Window, Every 5' mapping on specific gene is read and if 5'end is mapped within small window than added
    #             small_win += (int(entry[1]))
    
    #         if int(entry[0]) <= int(new_cleave_site+5) and int(entry[0]) >= int(new_cleave_site-5):## Range of Large window, GET LARGE AND SMALL WINDOW IN SAME LOOP
    #             large_win += (int(entry[1]))
    #     print('Small Window:',small_win, "| Large Window:",large_win)

    
    # ##else: ## Local
    # ##    print ('Looking up gene in density file............')
    # ##    try:
    # ##        found = den_dict[gene_name]
    # ##        print ('Gene found in density file\n')
    # ##        for i in found:###In gene summarization record first two entries are gene_name and size of gene, last entry is a blank line when record is fetched (not in file)
    # ##            locus = i.split()##position in summarization entry along with raw read, repeat normalized read and unique reads mapped, and category
    # ##
    # ##            
    # ##            if locus[0].isdigit() and int(locus[0]) <= int(cleave_site+3): ###Fixes BUG 14-3-13, found one line in summarized file as ['#', 'Transcriptome=./genomic_seq.fa'] instead of like ['12101480', '2', '1', '0', '2'] AND Reads file till large window limit+ 1 extra
    # ##                #print ('Entry :',locus, 'Cleave Site :',cleave_site)
    # ##                if int(locus[0]) <= int(cleave_site+1) and int(locus[0]) >= int(cleave_site-1):####Range of Small Window###Every 5' mapping on specific gene is read and if 5'end is mapped within small window than added
    # ##                    #print ('Using %s for small window' % locus[0])
    # ##                    small_win += (int(locus[1]))
    # ##                #else: ###if the 5' end is not within the window - pass
    # ##                #    pass
    # ##                #
    # ##                if int(locus[0]) <= int(cleave_site+7) and int(locus[0]) >= int(cleave_site-7):###Range of Large window###GET LARGE AND SMALL WINDOW IN SAME LOOP
    # ##                    #print ('Using %s for large window' % locus[0])
    # ##                    large_win += (int(locus[1]))
    # ##                
    # ##                #else: ###if the 5' end is not within the window - pass
    # ##                #    pass
    # ##            else: ## That means present position is more than extent of large window
    # ##                #print ('\nWe have passed the large window range so breaking\n')
    # ##                break
    # ##        
    # ##    #else: ###gene_name doesn't have a record in sumamrization file - tags in some lib mapped to gene but in '0_density' tags were missng and not mapped to any such gene
    # ##    ##So Gene is missing from normalized tags summarization file a.k.a 0_density
    # ##    except KeyError:
    # ##        print ('\n ** Warning - There is no entry for %s found in PARE summarization**\n' % (gene_name))
        
    # ##################################################################################################################
        
    # ## BugFix -  13-3-13 - In lib based validation, target is validated just raw reads mapped by cleaveland and sometimes there is no unique and normalized read mapped at corresponding site in '0_density' file so skip those
    # if small_win != 0:
    #     if Local == 'N':
    #         if extraAnno == 'Y':
    #             rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' % (','.join(ent),str(chr_id), strand, new_bind_site, new_cleave_site, small_win, large_win, round(small_win/large_win,2),genefunc,phasing,inverted))
    #             print ("Server:",rev_mapped_entry)
    #         else:
    #             rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s,%s' % (','.join(ent),str(chr_id), strand, new_bind_site, new_cleave_site, small_win, large_win, round(small_win/large_win,2),genefunc))
    #             print ("Server:",rev_mapped_entry)


    #     else: ## Local analysis has no Reverse mapping 
    #         rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%f\n' % (ent[0][1:],ent[1],ent[2], ent[3],ent[4],ent[6],ent[7],ent[5],small_win,large_win,round(small_win/large_win,2)))
    #         print ("Local:",rev_mapped_entry)
    
    # else:
    #     ##If the first result entry has 0 abundance for small window than no 'rev_mapped entry' variable and nothing to return, therefore added 'Rev_mapped_entry' with ratio set to 0
    #     #rev_mapped_entry = (('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0][1:],gene_name,new_bind_site, new_cleave_site,ent[4],ent[6],ent[7],ent[5],small_win,large_win, 0)))
    #     ##Or rev_mapped_entry just mentions the bug that caused blank entry which can be filtered later
    #     print('E13-3-13')
    #     rev_mapped_entry = ('E13-3-13')
    #     pass
    # # print ("mapping success")
    
    return rev_mapped_entry

def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    ##PP is being used for Bowtie mappings - This will avoid overflooding of processes to server
    nprocPP = round((nproc/int(nthread))+1) ## 1 added so as to avoid 0 processor being allocated in serial mode
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

def PPmultiple(module,alist1,alist2):
    start = time.time()
    npool = Pool(int(nproc))
    npool.map(lambda args: module(*args), alist2)

def PPResults(module,alist):
    npool = Pool(int(nproc))    
    res = npool.map_async(module, alist)
    if(res.get() == []):
        print("YEP!")
    results = (res.get())
    npool.close()           ### Added by Reza to close hanging 
    return results

def main():

    ### GET COORDS ####
    global con

    if Local == 'Y':
        if generateFasta == 'Y':
            coords = extractFeatures(genomeFile,gffFile) ## Extracts Coords from GFF3
            # fastaOut = getFASTALocal(genomeFile,coords) ##Creates FASTA file
            global tagLen ## Required later for tag2FASTA step as well
            # unambiguousBaseCounter(fastaOut,tagLen)
        else:
            print("\nThe input FASTA file is considered 'as is' for analysis\n")
            # fastaOut = genomeFile ### Make it better
    else:
        if generateFasta == 'Y':
            con = ConnectToDB(dataserver,0)
            coords = extractCoords(con,GenomeDB)###Get selected genes and all intergenic coords - all coords all also fetched too
            # fastaOut = getFASTAServer(con,GenomeDB,coords)###Extracts FASTA of supplied coords
            # unambiguousBaseCounter(fastaOut,tagLen)
        else:
            print("\nThe input FASTA file is considered 'as is' for analysis\n")
            con = ConnectToDB(dataserver,0)
            coords = extractCoords(con,GenomeDB) ##For reverse mapping
            # fastaOut = genomeFile ### Make it better
    

    ### REVERSE MAP ##########

    print('\n***Entering RevMapCoord- parallel***\n')
    global coord_dict_wat, coord_dict_crick
    coord_dict_wat = {} ## Dictionary of genes at watson strand, +
    coord_dict_crick = {}###Dictionary of genes at crick strand, - strand
    shutil.rmtree('./revMapped', ignore_errors=True) ## AK Added 
    os.mkdir('./revMapped') ## AK added

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
    
    #####Read the scoring input extend file and change coordinates ##########
    
    ResFls = [file for file in os.listdir('./') if file.endswith ('targs.parsed.csv')]
    
    for afile in ResFls:
        fh_in = open('./%s' % afile, 'r')### PARE VALIDATED results
        header = fh_in.readline() ## Waste header
        ScoInpExt = [] ## list of final results or parlist
        for res in fh_in:
            res_strp = res.strip('\n')
            ent =res_strp.split(',')
            ScoInpExt.append(ent)
        print ('**List from ScoInpExt ready to feed to RevMapCoord function**')
        
        ## Write results to file
        revmapRes = './revMapped/%s_revmapped.csv' % (afile)
        fh_out = open(revmapRes, 'w')
        fh_out.write('%s,Chr_id,Strand,Bind_start,Bind_site\n' % (header.strip('\n')))
        
        ## # TEST- SINGLE CORE - For troubleshooting ----##
        # ValidTarGeno = []
        # for i in ScoInpExt:
        #    z = RevMapCoord(i)
        #    ValidTarGeno.append(z)
        
        ## PARALLEL MODE - Uncommnet test above for normal use  
        print ('**Reverse mapping initiated**\n\n')
        ValidTarGeno = PPResults(RevMapCoord, ScoInpExt) ## Results are in form of list
        print ('**Reverse mapping complete for:%s\n\n\n' % (afile))

        for i in ValidTarGeno: ## Write Results from list to file
            if i != 'E13-3-13':##Filter for error 13-13-13 that is small window abundance = 0 and ratio calculation error
                fh_out.write('%s\n' % (i))
            else:
                print (i)
        fh_in.close()
        #fh_in2.close()
        fh_out.close()

### RUN ########################
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


## v01
## Written to revmap internal mirFErno predicted targets