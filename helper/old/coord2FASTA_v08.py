#!/usr/local/bin/python3
#### Script to get sequences from genome DB - written by kakrana@udel.com

import sys
import mysql.connector as sql

######## SEETINGS ################
server              = 'raichu.dbi.udel.edu' 	## Server to use to connect to DB
genomedb            = 'STREPTOCHAETA_ISU1_genome' 	## Genome DB to use for extraction
coords              = str(sys.argv[1])

SpecStrand          = 'Y' 						## Y: If strand is available and sequence needs to be from specific strand OR  F: No strand information and you need forward strand B: Than you dont care about strand and want from both
strandpos           = '3' 						## Column for strand - excel format - should be 'w' and 'c'

header              = 'Y' 						## Does file have header
sep                 = '\t' 						## Specify seprator used in spreadsheed or input file
chrpos              = '2' 						## Column for chromosome - excel format
startpos            =  '10' 					## Column for start postion- excel format
endpos              = '11' 						## Column for end position - excel format

## Buffer Settings
flankType           = 0 						## 0:Static 1:Dynamic
buff                = 100							## Length of static buffer
buffRatio           = 3							## Ratio of input sequence length to use for dynamic flanking seq from 5' and 3' 

##################################
##coords = 'V2_Abun_Clustered_Union_Phase_21_score_30.csv'### Coord format should be comma separated


#################################
############ MODULES ############

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

def Coords2FASTA(con,coords,genomedb,flankType,buff):
    
    '''Get FASTA from our servers, based on start, end coords and strand'''


    cur = con.cursor()
    
    file_out = coords.split('.')[0]+'.fa'
    fh_out = open(file_out,'w')
    
    fh_in = open(coords,'r')
    if header == 'Y':
        fh_in.readline()
    
    for i in fh_in:
      
        coord = i.strip('\n').split(sep)
        print ("\nThis is the entry:\n",coord)
        chr_id = str(coord[int(chrpos)-1]).strip().replace("chr","")

    	## Flanking sequnces dynamic or static #########
        if flankType == 0: ## Static
            buff = buff
        elif flankType == 1: ## Dynamic
            length = int(end)-int(start)+1
            buff = length*buffRatio
            print ('Length of sequence: %s | Length of buffer: %s' % (length,buff))
        else:
            print("Please input correct flank option - [D] Dynamic or [S] Static")
        
        ### Coordinates for 'w' and 'c' strands #######
        if SpecStrand == 'Y':
            strand = str(coord[int(strandpos)-1]).strip()
            
            ## Start coulumn could have end co-ordinates in case of crick strand. Best is to determine the start and end
            ## before any operation

            alist = [] ## To hold start end co-ordinates and pick min. value as start and max. as end
            start = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
            end = str(coord[int(endpos)-1]).strip()
            alist.append(start)
            alist.append(end)

            astart = min(alist)
            aend = (max(alist))

            start2 = int(astart)-buff
            if start2 < 0: ## Added to fix negative number query to my SQL
                start2 = 1
            else:
                pass

            length = (int(aend)-int(astart))+(buff*2)+1 ### extra 1 to include last position, buffer of 100 bp from start and end
            end2   = start2+length ## Not used for query but good to have 

        elif SpecStrand == "B" or SpecStrand == "F":
            start = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
            start2 = int(start)-buff ### It is assumned that start column is the start i.e. min value 
            if start2 < 0: ## Added to fix negative number query to my SQL
                start2 = 1
            else:
                pass

            length = (int(end)-int(start))+(buff*2)+1 ### xtra 1 to include last position, buffer of 100 bp from start and end
            end2   = start2+length ## Not used for query but good to have 

        else:
            print("Please input correct option for specStrand")
            pass


        ## Check if buffer surpasses the chromosome length ###########
        cur.execute("select chr_id,length from %s.chromosome_master where chr_id= %s" % (genomedb,chr_id))
        temp = cur.fetchall()
        chr_len = temp[0][1]
        # print("Length of this chromosome:%s" % (chr_len))
        if int(end2) > int(chr_len):
            print("End cordinate %s surpasses chr length %s , defauting to max possible length" % (end2,chr_len))
            end2 = chr_len ## End of chromosome
        else:
            pass

        ### Fetch Seqeunces ########################################    
        cur = con.cursor()
        if SpecStrand == 'Y': ##Strand info available - Like if coordinates are for gene
            print ('This is chr_id:%s and Strand: %s and Start:%s and End:%s | Length of sequence: %s bp\n' % (chr_id,strand,str(start2),int(end2),length))
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start2),length,genomedb,chr_id,strand))
            temp = cur.fetchall()
            #print (temp)
            #print (temp[1][0][::-1])
            if strand == 'w': ## Distinguish between strands so that you can reverce 'crick' strand to give correct orientation
                fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0]))
            else:
                fh_out.write('>%s_c_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0][::-1])) ### Just reverse to give correct orientation


        elif SpecStrand == 'F': ### No strand information forward strand required
            print ('This is chr_id:%s and Start:%s and End:%s | Length of sequence: %s bp\n' % (chr_id,str(start2),int(end2),length))
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = 'w'" % (str(start2),length,genomedb,chr_id))####
            temp = cur.fetchall()
            #print (temp)### Both positive and negative strand reported
            fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0]))

 
        elif SpecStrand == 'B': ### No strand information both strands were required
            print ('This is chr_id:%s and Start:%s and End:%s | Length of sequence: %s bp\n' % (chr_id,str(start2),int(end2),length))
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s" % (str(start2),length,genomedb,chr_id))####
            temp = cur.fetchall()
            #print (temp)### Both positive and negative strand reported
            fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0]))
            fh_out.write('>%s_c_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[1][0][::-1])) ### Reversed to give correct orientation

        else:
            print("Choose correct option for 'SpecStrand' parameter")
            sys.exit()

    fh_out.close()
   
def main(genomedb,coords):
    con = ConnectToDB(server)
    Coords2FASTA(con,coords,genomedb,flankType,buff)

if __name__ == '__main__':
    main(genomedb,coords)
    sys.exit()


####V01 -> v02 11-Oct-13
###Inputs added for strand, seprator, and other column position sin input file
### If coords has strand info than choose 'SpecStrand' mode and specify 'strandpos' variable i.e position in excel or soreadsheet
### If coords info not available than 'SpecStrand' should be 'N' in this cse sequence from both strands reported back

## V02 -> v03
## Strand if not specified or used does not give an error as its moved to if stanrd specified section

## v03 -> v04
## Add dynamic flanking region assignment for each coord based on length of coordinates in input file

## v04 -> v05
## Added correct end with buffer to seq header
## Fixed if start or end surpass chromosome limits after adding buffer. In such cases set start to 0 and end to last of chromosme

## v05 -> v06
## Corrected extra addition of buff to length while fetching seqeunce

## v06->v07
## Added functionality to get from watson strand if strand is not known

## v07 -> v08 [Major fix]
## Crick strand is fetched correctly, script resolves start and end from given values, i.e. smaller value is start and large is end for the sake of mySQL
