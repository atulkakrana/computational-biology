#!/usr/local/bin/python3

### Script to get Annotation using coordinates

import sys
import mysql.connector as sql

#### Input ##########################
coordsFile = 'Final_PHAS_Loci_1e-07_21PHAS.csv'
phasLoci = 'Y' ### Is it a phased loci i.e two coords start and end

server = 'raichu.dbi.udel.edu'
db = 'RICE_MSU7_genome'

sep = '\t'
strand = 'N' ## Is strand info available as in Ping degradome validated phased loci
strandCol = 3 ### In 'w' and 'c' format

nameCol = 1 ## In excel format
chrCol = 3
startCol = 4
endCol = 5


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

def getAnnotation(con):
        
    cur = con.cursor()
    ### Test
    #cur.execute("SELECT gene,chr_id,strand,start,end,title,type FROM %s.gene_master where chr_id = 1 and strand ='w' and start <= 250000" % (db))
    ### Real
    cur.execute("SELECT gene,chr_id,strand,start,end,title,type FROM %s.gene_master" % (db))
    annoList = cur.fetchall() #### List to hold tuples of annotation information
    
    outfile = '%s.Anno' % (coordsFile)
    fh_out = open(outfile, 'w')
    fh_out.write('Name\tGene\tchr\tstrand\tTitle\n')
    
    fh_in = open(coordsFile,'r')
    fh_in.readline()
    
    for z in fh_in:
        ent = z.strip('\n').split(sep)
        name = ent[nameCol-1]
        chr_id = ent[chrCol-1]
        start = ent[startCol-1]
        end = ent[endCol-1]
        strand = ent[strandCol-1]
        if strand == 'Y':
            print ('\nPhased Loci:',name,chr_id,strand,start,end)
        else:
            print ('\nPhased Loci:',name,chr_id,start,end)
            
        anno = [] ## Temp list to store annotaton for each entry
        
        if strand == 'Y':
            for i in annoList:
                #print (chr_id,strand,start,end)
                #print (i[1],i[2],i[3],i[4],i[6])
                if int(chr_id) is (i[1]) and str(strand) is str(i[2]) and (i[3]) <= int(start) and int(end) <= (i[4]):
                    print ('Found',i)
                    anno.append((name,i[0],i[1],i[2],i[5],i[6])) ## Name of loci, genename, chr_id, strand, title, type
                    #fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,i[0],i[1],i[2],i[5],i[6]))
                    break
                
                elif int(chr_id) == i[1] and strand == i[2] and i[3] <= int(start) and int(start) <= i[4]:
                    print ('Found by start only',i)
                    anno.append((name,i[0],i[1],i[2],i[5],i[6]))
                    #fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,i[0],i[1],i[2],i[5],i[6]))
                    break
                    
                elif int(chr_id) == i[1] and strand == i[2] and i[3] <= int(end) and int(end) <= i[4]:
                     print ('Found by end only',i)
                     anno.append((name,i[0],i[1],i[2],i[5],i[6]))
                     #fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,i[0],i[1],i[2],i[5],i[6]))
                     break
                else:
                    #print ('Not Found')
                    pass
        
        else: ## No strand available
            
            for i in annoList:
                if int(chr_id) == i[1] and i[3] <= int(start) and int(end) <= i[4]:
                    print ('Found',i)
                    anno.append((name,i[0],i[1],i[2],i[5],i[6])) ## Name of loci, genename, chr_id, strand, title, type
                    #fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,i[0],i[1],i[2],i[5],i[6]))
                    break
                
                elif int(chr_id) == i[1] and i[3] <= int(start) and int(start) <= i[4]:
                    print ('Found by start only',i)
                    anno.append((name,i[0],i[1],i[2],i[5],i[6]))
                    #fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,i[0],i[1],i[2],i[5],i[6]))
                    break
                    
                elif int(chr_id) == i[1] and i[3] <= int(end) and int(end) <= i[4]:
                     print ('Found by end only',i)
                     anno.append((name,i[0],i[1],i[2],i[5],i[6]))
                     #fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,i[0],i[1],i[2],i[5],i[6]))
                     break
                else:
                    #print ('Not Found')
                    pass
        
            
        
        ### Write result for each entry
        if anno: ## If annotation found
            fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (name,i[0],i[1],i[2],i[5],i[6]))
        else: ## If annotation not found
            fh_out.write('%s\tinter\t%s\t%s\tinter\tinter\n' % (name,chr_id,start))
            
    fh_in.close()
    fh_out.close()
    
    return outfile
    

def main(server,db):
    con = ConnectToDB(server)
    outfile = getAnnotation(con)
    
if __name__ == '__main__':
    main(server,db)
    print ('***Script finished successfuly***\n')
    sys.exit()
    
    
#### LOG ####
### 21 Nov 2013 -> v01
## V01 -> V02
### Added annotation without strand
    

    
    
    