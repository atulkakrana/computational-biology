#!/usr/local/bin/python3

## Version 0.2 - 06/30/15
## Written by Atul Kakrana
## This strategy of reverse mapping is hardcoded for max 5 multiple hits - Change npaths value if required


''' This script takes co-ordinates and maps the seqeunces to new genome
to get new co-ordinates'''

import os,sys,subprocess,multiprocessing,datetime,time,shutil
import mysql.connector as sql

######## SEETINGS ################
server      = 'raichu.ddpsc.org'         ## Server to use to connect to DB
genomedb    = 'ASPARAGUS_UGA1_genome'       ## Genome DB to use for extraction
newGenome   = "/home/kakrana/99.genomes/asparagus/v2/modified/Asparagus.V2.0.genome.stable_modifed.fa"
coords      = str(sys.argv[1])
                                            ## New Genome to map to get new co-ordinates, file must be in . directory

gmapDB      = 0                             ## If 0: Do not make DB 1: Make the DB before mapping
bothStrands = 1                             ## Hits will be reported for both strands, like in case when fragment to be remapped has moved to other strand
SpecStrand  = 'Y'                           ## Y: If strand is available and sequence needs to be from specific strand | F: No strand information and you need forward strand | B: Than you dont care about strand and want from both - This mode is not used in this script
strandpos   = 58                            ## Column for strand - excel format - should be 'w' and 'c'

header      = 'Y'                           ## Does file have header
sep         = '\t'                          ## Specify seprator used in spreadsheed or input file
namepos     = 53
chrpos      = 55                             ## Column for chromosome - excel format
startpos    = 56                             ## Column for start postion- excel format
endpos      = 57                             ## Column for end position - excel format

## Buffer Settings
flankType   = 0                             ## 0:Static 1:Dynamic
buff        = 0                             ## Length of static buffer
buffRatio   = 3                             ## Ratio of input sequence length to use for dynamic flanking seq from 5' and 3' 
numProc     = 0
##################################

def ConnectToDB(server):

    '''Get Connection to the server'''
    
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
    
    '''Get (query) FASTA seqeunces from the co-ordinates'''

    cur = con.cursor()
    
    file_out = coords.split('.')[0]+'.fa'
    fh_out = open(file_out,'w')
    
    fh_in = open(coords,'r')
    if header == 'Y':
        fh_in.readline()
    
    nameList = [] ## Store names for uniqBED function
    for i in fh_in:
      
        coord = i.strip('\n').split(sep)
        print ("\nThis is the entry:\n",coord)
        start = str(coord[int(startpos)-1]).strip() ##Buffer of 100bp or as specified
        end = str(coord[int(endpos)-1]).strip()
        name = str(coord[int(namepos)-1]).strip()
        nameList.append(name)
        # print(name,start,end)

        ## Flanking sequnces dynamic or static
        if flankType == 0: ## Static
            buff = buff
        elif flankType == 1: ## Dynamic
            length = int(end)-int(start)+1
            buff = length*buffRatio
            print ('Length of sequence: %s | Length of buffer: %s' % (length,buff))
        else:
            print("Please input correct flank option - [D] Dynamic or [S] Static")
        
        start2 = int(start)-buff

        if start2 < 0: ## Added to fix negative number query to my SQL
            start2 = 1
        else:
            pass

        length = int(end)-int(start)+buff*2 ###Extra 1 to include last position, buffer of 100 bp from start and end
        chr_id = str(coord[int(chrpos)-1]).strip()
        end2 = start2+length
        
        cur.execute("select chr_id,length from %s.chromosome_master where chr_id= %s" % (genomedb,chr_id))
        temp = cur.fetchall()
        chr_len = temp[0][1]
        print(chr_len)
        if end2 > chr_len:
            end2 = chr_len ## End of chromosome
        else:
            pass
        
        cur = con.cursor()
        if SpecStrand == 'Y': ##Strand info available - Like if coordinates are for gene
            strand = str(coord[int(strandpos)-1]).strip().translate(str.maketrans("+-","wc"))
            print ('This is chr_id:%s and Strand: %s and Start:%s and End:%s | Length of sequence: %s bp\n' % (chr_id,strand,str(start2),int(end2)+buff,length))
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s AND strand = '%s'" % (str(start2),length,genomedb,chr_id,strand))
            temp = cur.fetchall()
            #print (temp)
            #print (temp[1][0][::-1])
            if strand == 'w': ## Distinguish between strands so that you can reverce 'crick' strand to give correct orientation
                # fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0]))
                fh_out.write('>%s\n%s\n' % (name,temp[0][0]))

            else:
                # fh_out.write('>%s_c_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0][::-1])) ### Just reverse to give correct orientation
                fh_out.write('>%s\n%s\n' % (name,temp[0][0][::-1]))


        elif SpecStrand == 'F': ### No strand information forward strand required
            print ('This is chr_id:%s and Start:%s and End:%s | Length of sequence: %s bp\n' % (chr_id,str(start2),int(end2)+buff,length))
            cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s" % (str(start2),length,genomedb,chr_id))####
            temp = cur.fetchall()
            #print (temp)### Both positive and negative strand reported
            # fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0]))
            fh_out.write('>%s\n%s\n' % (name,temp[0][0]))

 
        # elif SpecStrand == 'B': ### No strand information both strands were required
        #     print ('This is chr_id:%s and Start:%s and End:%s | Length of sequence: %s bp\n' % (chr_id,str(start2),int(end2)+buff,length))
        #     cur.execute("select substring(chromosome, %s, %s) from %s.chrom_sequence where chr_id = %s" % (str(start2),length,genomedb,chr_id))####
        #     temp = cur.fetchall()
        #     #print (temp)### Both positive and negative strand reported
        #     fh_out.write('>%s_w_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[0][0]))
        #     fh_out.write('>%s_c_%s:%s:%s\n%s\n' % (coord[0],chr_id,start2,end2,temp[1][0][::-1])) ### Reversed to give correct orientation

        else:
            print("Choose correct option for 'SpecStrand' parameter")
            sys.exit()

    fh_out.close()
    print("FASTA file %s generated with seqeunces to be mapped\n\n" % (file_out))

    return file_out

def mapper(fastaFile,newGenome):
    ''' mapped fasta file using gmap and returns 
    SAM file'''

    if gmapDB == 1:
        print("Creating a database of new genome - Relax")
        newGenomeDB = "%s.db" % (newGenome.split("/")[-1].rpartition('.')[0])
        db_folder = "%s/gmap_temp" % (os.getenv('HOME'))
        shutil.rmtree(db_folder,ignore_errors=True)
        os.mkdir(db_folder)

        print("FASTA Name: %s | DB folder: %s | DB Name: %s" % (newGenome,db_folder, newGenomeDB))
        retcode = subprocess.call(["gmap_build", "-D", db_folder, "-d", newGenomeDB, newGenome])
        
        if retcode == 0:
            print("The DB for new genome made sucessfully")
        else:
            print("There is some problem with the DB build - System will exit now")
            sys.exit()

    else:
        newGenomeDB = "%s.db" % (newGenome.split("/")[-1].rpartition('.')[0])
        db_folder = "%s/gmap_temp" % (os.getenv('HOME'))
        print("Will use gmap DB %s folder - Make sure it exists or run with gmapDb = 1" % (db_folder))

    print("\nMapping will start now - Relax")

    samName = "%s.sam" % fastaFile.rpartition('.')[0] 
    logName = "%s.log" % fastaFile.rpartition('.')[0]
    perc_id = 1.0
    npaths = 10 
    print("Working directory: %s" % (os.getcwd()))

    if bothStrands == 1:
        x = subprocess.Popen("gmap -D %s -d %s -f samse --mapboth --nosplicing --no-chimeras --ordered --min-identity %s -t %s -n %s %s > %s 2> %s" % (db_folder,newGenomeDB,str(perc_id),str(nproc),str(npaths),fastaFile,samName,logName),shell=True)
        x.communicate()
    elif bothStrands == 0:
        x = subprocess.Popen("gmap -D %s -d %s -f samse --nosplicing --no-chimeras --ordered --min-identity %s -t %s -n %s %s > %s 2> %s" % (db_folder,newGenomeDB,str(perc_id),str(nproc),str(npaths),fastaFile,samName,logName),shell=True)
        x.communicate()
    else:
        print("Please enter correct value of both strands - Script will exit")
        sys.exit()

    # retcode1 = subprocess.call(["gmap", "-D", db_folder, "-d", newGenomeDB, "-f", "samse", "--nosplicing", "--no-chimeras", "--ordered", "--min-identity", str(perc_id), "-t", str(nproc), "-n", str(npaths), fastaFile, ">", samName, "2>", logName])

    # if retcode1== 0:
    #     print("Working directory: %s" % (os.getcwd()))
    #     print("Mapping to new genome finished sucessfully\n")
    # else:
    #     print("There was some problem with mapping - System will exit now\n")
    #     sys.exit()

    return samName  

def sam2BED(samName):

    bamFile = "%s.bam" % samName.rpartition('.')[0]
    print("Conversion of SAM file %s to BAM file %s will start now - Relax\n" % (samName,bamFile))

    retcode = subprocess.call(["samtools","view", "-bS", samName, "-o", bamFile])
    
    if retcode == 0:
        print("Conversion of SAM to BAM finished sucessfully")
    else:
        print("There was some problem with SAM to BAM conversion - System will exit now")
        sys.exit()

    
    bedFile = "%s.bed" % bamFile.rpartition('.')[0]
    print("Conversion of BAM file %s to BED file %s will start now - Relax\n" % (bamFile,bedFile))

    bed = subprocess.Popen("bamToBed -i %s > %s" % (bamFile,bedFile),shell=True)

    # retcode1 = subprocess.Popen(["bamToBed", "-i", bamFile, "-o", bedFile])

    # print("Conversion of BAM to BED will start now - Relax")
    # if retcode1 == 0:
    #     print("Conversion of BAM to BED finished sucessfully")
    # else:
    #     print("There was some problem with BAM to BED conversion - System will exit now")
    #     sys.exit()

    print ("\nFor Final mapping file please see: %s" % (bedFile))
    time.sleep(3)

    return bedFile

def uniqBED(coords,bedFile):
    ''' This function will uniq the mappings if there are multiple hits'''

    print("Function: uniqBED")

    ## Outfile
    uniqFile = "%s.uniq.bed" % bedFile.rpartition('.')[0]
    fh_out = open(uniqFile,'w')

    ## Get query names
    fh_in = open(coords,'r')
    if header == 'Y':
        fh_in.readline()
    
    nameList = [] ## Store names for uniqBED function
    for i in fh_in:
        coord = i.strip('\n').split(sep)
        # print ("\nThis is the entry:\n",coord)
        name = str(coord[int(namepos)-1]).strip()
        nameList.append(name)
    print("List of %s names prepared to generate Uniq results" % (len(nameList)))
    print("snippet of the nameList:%s" % (nameList[1:5]))

    ## Get mappings from bedfile and test if this mapping assigned before ###
    assignedMaps = [] ## A list of sets with chr,start,end,strand for already assigned mappings
    bedList = bedParser(bedFile)

    acount = 0 ## Count of entries written to uniq bedFile
    for getname in nameList:
        for ent in bedList:
            # print(ent)
            if getname == ent[3]:
                print("Mapping found for %s:%s" % (getname,ent))
                achr,astart,aend,aname,aflag,strand = ent
                astrand = strand.translate(str.maketrans("+-","wc"))
                print(achr,astart,aend,astrand)
                akey = "%s-%s-%s-%s" % (achr,astart,aend,astrand)

                if akey not in assignedMaps:
                    # print("This is an unassigned mapping")
                    fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (achr,astart,aend,aname,aflag,astrand))
                    assignedMaps.append(akey)
                    acount+=1
                    break ## Mapping found for this entry no need to record further mappings if exits for this
                else:
                    print("This key %s already assigned to other entry" % (akey))
                    pass
            else:
                # print("Mapping is not for %s" % (getname))
                pass

    print("Entries in coord file:%s | Entries with uniq mappings assigned:%s" % (len(nameList),acount))

    fh_out.close()
    
    return None

def bedParser(bedFile):
    ''' parse the bedFile'''

    bedList = [] ## List to store results

    print("Reading BED file:%s" % (bedFile))
    fh_in   = open(bedFile,'r')
    bedRead = fh_in.readlines()
    fh_in.close()

    print("Parsing BED file:%s" % (bedFile))
    for i in bedRead:
        # print(i)
        achr,astart,aend,aname,aflag,astrand = i.strip("\n").split("\t")
        print(achr,astart,aend,aname,aflag,astrand)
        bedList.append((achr,astart,aend,aname,aflag,astrand))

    print("Entries parsed from bed file:%s" % len(bedList))
    
    return bedList   


def main(genomedb,coords):
    con = ConnectToDB(server)
    fastaFile = Coords2FASTA(con,coords,genomedb,flankType,buff)
    # fastaFile = "21-22NovelPriv.fa"
    samName = mapper(fastaFile,newGenome)
    bedFile = sam2BED(samName)
    # bedFile = "21-22UGA.bed"
    uniqBED(coords,bedFile)

if __name__ == '__main__':
    if numProc == 0:
        nproc = int(multiprocessing.cpu_count()*0.90)
    else:
        nproc = int(numProc)

    main(genomedb,coords)
    sys.exit()


## v01 -> v02
## Gmap now runs through Popen rather then call because stdout was missing writing to a file

## v02 -> v03 [major][stable]
## Fixed bug to generate new database and added bothStand mode
## Fixed an error with +/- strand conversion, it was being converted to W/C instead of w/c
## Added new functions that uniq the mappings for each entry
## Fixed bedFile was not being read because it takes time to get ready, and python moves forward detecting an empty file