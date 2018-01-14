#!/usr/local/bin/python3
#### Script to get sequence for the gene identifier from genome DB - written by kakrana@udel.edu

import sys
import mysql.connector as sql
from operator import itemgetter

#####
server      = 'raichu.ddpsc.org'
genomedb    = 'RICE_MSU7_genome'
idFile      = 'RICE_MSU7_genome_gene_master.ids'   ## This file can be generated for whole genome - mysql -u kakrana -p -e "select gene,start,end,chr_id,strand from CICER_LIS1_genome.gene_master" > CICER_LIS1_genome_Genes.txt 
                                                

#####
header      = 'Y' ## ID FILE has header or not
sep         = '\t' ## Specify seprator used in spreadsheed or input file
geneID      = '1' ## Column with gene ID in excel format
seqType     = '2' ## 1 = CDS , 2 = Transcript (Minus introns, only model 1 of gene used), 3 = Complete gene + buffer (specify), 4 = 5' promoter of gene as specifed by buff5 below (including 5' UTR)
buff5       = 1000 ## Only for complete gene i.e.seqType == 3 and 4
buff3       = 500 ## Only for complete gene i.e.seqType == 3
#####


def ConnectToDB(server):
    
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

def genes2CDS(con,idFile,genomedb):
    print(' mode selected - Exons only')
    
    cur = con.cursor()
    
    file_out = idFile.split('.')[0]+'.fa'
    fh_out = open(file_out,'w')
    
    fh_in = open(idFile,'r')
    if header == 'Y':
        fh_in.readline()

    for anent in fh_in:
        anent_splt = anent.split(sep)
        gene = anent_splt[int(geneID)-1]
        print ('\nExtracting sequence for ID: ',gene)
        
        ## Select the lowest avilable model for gene
        cur.execute ("SELECT min(abs(model)) FROM %s.gene_position WHERE gene = '%s'" % (genomedb,gene)) ## Select a model if several available for gene for sequence extraction, T0 entries of a gene have negative model so removed from here
        model = cur.fetchall()[0][0]
        
        cur = con.cursor()
        cur.execute("SELECT * FROM %s.gene_position WHERE gene = '%s' AND model = %s AND utr = '5'" % (genomedb,gene,'-' + str(model))) ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        utr5_info = cur.fetchall()
        cur.execute("SELECT * FROM %s.gene_position WHERE gene = '%s' AND model = %s AND utr = '3'" % (genomedb,gene,'-' + str(model))) ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        utr3_info = cur.fetchall()
        cur.execute("SELECT * FROM %s.gene_position WHERE gene = '%s' AND model = %s AND utr = 'E'" % (genomedb,gene,'-' + str(model))) ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        utrE_info = cur.fetchall()
        cur.execute("SELECT * FROM %s.gene_position WHERE gene = '%s' AND model = %s AND utr = 'N'" % (genomedb,gene,model)) ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        gene_info = cur.fetchall() ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        #print('Gene Info:',gene_info)
        strand = gene_info[0][7]
        chr_id = gene_info[0][6]
        ## First remove the exon that has become UTR nnow i.e. utr E

        if utrE_info:
            print('This is the extended UTR',utrE_info)
            gene_info2 =[] ## List without exons that match extended UTR
            exon_remove = [] ##List of conflicting exons to remove
            for i in utrE_info:
                exon_remove.append(i[2])
                for exon in gene_info:
                    if exon[2] not in exon_remove:
                        gene_info2.append(exon)
                    else:
                        pass
        
        else:
            gene_info2 = list(gene_info) ## copy original gene_info to gene_info2 to maintain same name below
            
        print('This is gene info2:',gene_info2)
                        

        ## Sort Exons for concatnating
        info_sorted = sorted(gene_info2, key=itemgetter(2)) ## Sorted on exon ids - N, 5 and 3
        #print ('This is gene info2 sorted as per exon numbers:',info_sorted)
        cds = '' ## Empty Seq
        
        ## First and last exons incude UTRs if present. So, remove those UTR seq from first and last exons
        if strand == 'w':
            if utr5_info:
                print("5' UTR info:", utr5_info)
                alst = list(info_sorted[0]) ## First exon info converted to list
                utr5_end =  utr5_info[0][5] ## End of 5' UTR is the start of exaon1
                alst[4] = utr5_end+1 ## Exon start updated with 5' utr end, 1 added as cds starts next nt of end
                info_sorted[0] = tuple(alst) ## The first exon list is reconverted to tuple and updated to gene sorted info
                
            if utr3_info:
                print("3' UTR info:", utr3_info )
                blst = list(info_sorted[-1]) ## Last exon converted to list
                utr3_start = utr3_info[0][4] ## Start of 3' UTR
                blst[5] = utr3_start - 1 ## Last exon end updated with 3' UTR start, 1 substracted as cds ends 1nt before utr starts
                info_sorted[-1]=tuple(blst) ## The last exon list is reconverted to tuple and updated to gene sorted info
        
        elif strand =='c':
            if utr5_info:
                print("5' UTR info:", utr5_info)
                alst = list(info_sorted[0]) ## First exon info converted to list
                utr5_start =  utr5_info[0][4] ## Startof 5' UTR is the start of exon1
                alst[5] = utr5_start -1 ## Exon end updated with 5' utr start, 1 substracted before exon ends 1 nt before exon start
                info_sorted[0] = tuple(alst) ## The first exon list is reconverted to tuple and updated to gene sorted info
            
            if utr3_info:
                print("3' UTR info:", utr3_info)
                blst = list(info_sorted[-1]) ## Last exon converted to list
                utr3_end = utr3_info[0][5] ## End of 3' UTR
                blst[4] = utr3_end + 1 ## Last exon start updated with 3' UTR end, 1 added becasue last exon starts 1nt after utr end
                info_sorted[-1]=tuple(blst) ## The last exon list is reconverted to tuple and updated to gene sorted info
        else:
            print('What Strand is it man?')

        ## After UTRs are removed from first and last exon, now extract seqeunce
        for exon in info_sorted:
            print('Exon:', exon)
            start = exon[4]
            end = exon[5]
            length = int(end)-int(start)+1 ## 1 is added to adjust for substraction
            
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start),str(length),genomedb,chr_id,strand))####
            temp = cur.fetchall()
            #print(exon,temp)
            
            if strand == 'w':
                cds+=temp[0][0]
            elif strand == 'c':
                cds+=temp[0][0][::-1]
            else:
                print('What Strand is it man?')
            #print('This is cds:', cds)
        fh_out.write('>%s\n%s\n' % (gene,cds))
        
    fh_out.close()
    fh_in.close()
    
    return file_out

def genes2Trans(con,idFile,genomedb):
    
    print('Gene mode selected - UTR+exon+3UTR')
    
    cur = con.cursor()
    
    file_out = idFile.split('.')[0]+'.fa'
    fh_out = open(file_out,'w')
    
    fh_in = open(idFile,'r')
    if header == 'Y':
        fh_in.readline()

    for anent in fh_in:
        anent_splt = anent.split(sep)
        gene = anent_splt[int(geneID)-1]
        print ('\nExtracting sequence for ID: ',gene)
        
        cur = con.cursor()
        cur.execute ("SELECT min(abs(model)) FROM %s.gene_position WHERE gene = '%s' AND utr != 'T0'" % (genomedb,gene)) ## Select a model if several available for gene for sequence extraction, T0 entries of a gene have negative model so removed from here
        model = cur.fetchall()[0][0]
        print(model)
        cur.execute("SELECT * FROM %s.gene_position WHERE gene = '%s' AND model = %s AND utr = 'E'" % (genomedb,gene,'-' + str(model))) ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        utrE_info = cur.fetchall()
        cur.execute("SELECT * FROM %s.gene_position WHERE gene = '%s' AND utr = 'N' AND model = %s ORDER BY utr" % (genomedb,gene,model)) ## UTR N - Already contains 5', 3' and Extended UTR (as redundant exon entry must be present)
        gene_info = cur.fetchall() ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        print('Gene Info:',gene_info)
        strand = gene_info[0][7]
        chr_id = gene_info[0][6]
        #
        #if utrE_info:
        #    print('This is the extended UTR',utrE_info)
        #    gene_info2 =[] ## List without exons that match extended UTR
        #    exon_remove = [] ##List of conflicting exons to remove
        #    for i in utrE_info:
        #        exon_remove.append(i[2])
        #        for exon in gene_info:
        #            if exon[2] not in exon_remove:
        #                gene_info2.append(exon)
        #            else:
        #                pass
        #
        #else:
        #    gene_info2 = list(gene_info) ## copy original gene_info to gene_info2 to maintain same name below
        #    
        #print('This is gene info2:',gene_info2)


        info_sorted = sorted(gene_info, key=itemgetter(2)) ## On exon id
        print ('This is gene info sorted as per exon numbers:',info_sorted)
        trans = '' ## Empty Seq
        for exon in info_sorted:
            start = exon[4]
            end = exon[5]
            length = int(end)-int(start)+1 ## 1 is added to adjust for substraction
            
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start),str(length),genomedb,chr_id,strand))####
            temp = cur.fetchall()
            if strand == 'w':
                trans+=temp[0][0]
            else:
                trans+=temp[0][0][::-1]
                
            #print('This is cds:', cds)
        fh_out.write('>%s\n%s\n' % (gene,trans))
            
    fh_out.close()
    fh_in.close()
    
    return file_out

def gene2FASTA(con,idFile,genomeDB):
    
    print('Whole gene mode selected - 1kb + UTR + exon + Intron + 3UTR + 500 bases')
    
    cur = con.cursor()
    
    file_out = idFile.split('.')[0]+'.fa'
    fh_out = open(file_out,'w')
    
    fh_in = open(idFile,'r')
    if header == 'Y':
        fh_in.readline()

    for anent in fh_in:
        geneSeq ='' ## Gene seqeunce will be stored in this variable
        anent_splt = anent.split(sep)
        gene = anent_splt[int(geneID)-1]
        print ('\nExtracting sequence for ID: ',gene)
        
        ##Get lowest model number for gene
        cur.execute ("SELECT min(model_cnt) FROM %s.gene_master WHERE gene = '%s'" % (genomedb,gene)) ## Select a model if several available for gene for sequence extraction, T0 entries of a gene have negative model so removed from here
        model = cur.fetchall()[0][0]
        
        ## Get coords for gene
        cur.execute("SELECT * FROM %s.gene_master WHERE gene = '%s' AND model_cnt = %s" % (genomedb,gene,model)) ## ORDER BY will give N,5,3
        gene_info = cur.fetchall() ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        print('Gene Info:',gene_info)
        strand = gene_info[0][2]
        chr_id = gene_info[0][1]
        start = int(gene_info[0][3])
        end = int(gene_info[0][4])

        
        ## Get sequence of the gene with promoter (buffer as specified in the user settings)
        if strand == 'w':
            start2 = start-buff5
            end2 = end+buff3
            length = int(end2)-int(start2)+1  
            
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start2),str(length),genomedb,chr_id,strand))####
            temp = cur.fetchall()
            geneSeq+=temp[0][0]
        
        elif strand == 'c':
            start2 = start - buff3 ## start corresponds to 3' i.e. real end
            end2 = end + buff5 ## End corrsposnds to 5' i.e. real start
            length = int(end2)-int(start2)+1 
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start2),str(length),genomedb,chr_id,strand))####
            temp = cur.fetchall()
            geneSeq+=temp[0][0][::-1]
        
        else:
            print('What Strand is it man?')
        
        fh_out.write('>%s\n%s\n' % (gene,geneSeq))
        
    fh_out.close()
    fh_in.close()
        
    return file_out

def genePromo(con,idFile,genomeDB):
    print ('You choose to fetch the promoter of gene - 1.5 Kb + 5 UTR')
    
    cur = con.cursor()
    
    file_out = idFile.split('.')[0]+'.fa'
    fh_out = open(file_out,'w')
    
    fh_in = open(idFile,'r')
    if header == 'Y':
        fh_in.readline()

    for anent in fh_in:
        promSeq ='' ## Variablel to hold promoter sequence
        geneSeq ='' ## Gene seqeunce will be stored in this variable
        anent_splt = anent.split(sep)
        gene = anent_splt[int(geneID)-1]
        print ('\nExtracting sequence for ID: ',gene)
        
        ## Get lowest model number for gene
        cur.execute ("SELECT min(model) FROM %s.gene_position WHERE gene = '%s'" % (genomedb,gene)) ## Select a model if several available for gene for sequence extraction, T0 entries of a gene have negative model so removed from here
        model = cur.fetchall()[0][0]
        
        ## Get the first Exon/gene start as UTRs are not always annoatated
        cur.execute("SELECT * FROM %s.gene_position WHERE WHERE gene = '%s' AND model = %s AND utr = 'N'" % (genomedb,gene,model)) ## gene, model,exon_id,utr,strat,stop,chr_id,strand
        gene_info = cur.fetchall[0][0]
        print('Gene Info:',gene_info)
        strand = gene_info[0][7]
        chr_id = gene_info[0][6]
        
        ## Choose start co-ordinate fo first exon
        info_sorted = sorted(gene_info, key=itemgetter(2)) ## On exon id
        
        if strand == 'w':
            coding_start = info_sorted[0][4] ## First exon start co-ordinate
            start = coding_start-buff5
            length = buff5
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start),str(length),genomedb,chr_id,strand))####
            temp = cur.fetchall()
            promSeq+=temp[0][0]
            fh_out.write('>%s_%sprom\n%s\n' % (gene,buff5,cds))
            
        elif strand == 'c':
            coding_start = info_sorted[0][5] ## First exon end co-ordinate 
            start = coding_start ## Extract fom this point next few hundred nts as in buff5
            length = buff5
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start),str(length),genomedb,chr_id,strand))####
            temp = cur.fetchall()
            promSeq+=temp[0][0][::-1] ## So that orientation of promoter is changed from 0 to buff5 (1-2kb or whatever) to -buff5 to 0
            fh_out.write('>%s_%sprom\n%s\n' % (gene,buff5,cds))
        
        else:
            print('What Strand is it man?')

    fh_out.close()
    fh_in.close()
    
    return file_out

def main(genomedb,idFile):
    con = ConnectToDB(server)
    if seqType == '1':
        resFile = genes2CDS(con,idFile,genomedb)
    elif seqType == '2':
        resFile = genes2Trans(con,idFile,genomedb)
    elif seqType == '3':
        resFile = gene2FASTA(con,idFile,genomedb)
    elif seqType == '4':
        resFile = genePromo(con,idFile,genomedb)
    else:
        print('Please input the correct sequence type in the user option')
        print('Current available options are: 1 = CDS and 2 = Transcripts')
        
    print('\nThe results are in "%s"' % (resFile))
    
if __name__ == '__main__':
    main(genomedb,idFile)
    print ('\n**Script finished sucessfully**')
    
    sys.exit()
    
    

## Author Remarks
## Spent a day of my (June 9th,2014) life and PhD. writing this script - Once for fucking all

## v01
# Extract gene sequences - Use gene identifier to extract sequences of CDS, CDS with UTR and whole gene with promoters
