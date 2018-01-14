
#!/usr/local/bin/python3

## Created by Parth Patel, DBI @ University of Delaware, Newark, Delaware 19717
## Completly changed by ATul Kakrana :)
## Requires HPAnnotate.v01.pl or above

import sys,os,subprocess,string,time

### USER VARIABLES
runMode       =  'M'              ## A: Automatic - Provide list of ids, FASTA file which has precursors, these will extracted aotomatically 
                                    ## and phasiRNAs file from phasiExtract. phasiRNAs for corresponding transcript will be extracted automatically
                                    ## M: Manual mode frovide FASTA file with singl eprecursor and phasiRNA file with corresponding phasiRNA, both should in FASTA format

idFile          = 'plot.list'       ## Required for mode 'A'
fastaFile       = '24phas.fa'       ## Required for mode 'A'
phasiFile       = 'ALL.txt.score_p1e-07_sRNA_24_out.clusterthres_0.95_phasi.csv'
mode            = 0                 ## 0: match direct reads only, 1: Match rev comp too, 3: Try to match direct, rev comp and comp
cleanup         = 1                 ## Delete temp files

#### FUNCTIONS ###########
def fastaReader(fastaFile):
    ''' Read FASTA and give a list of tuples with name and seqeuences'''

    fh_in       = open(fastaFile,'r')
    fasta       = fh_in.read()
    fasta_splt  = fasta.split('>')
    
    acount      = 0 ## count the number of entries
    empty_count = 0 ## To count empty 
    fastaList   = []
    
    for i in fasta_splt[1:]:
        ent     = i.split('\n')
        name    = ent[0].split()[0].strip()
        seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        alen    = len(seq)
        acount  +=  1

        if seq != '':
            fastaList.append((name,seq))
        else:
            empty_count+=1
            pass

    fh_in.close()

    print("Here is an examlpe of FASTAList:",fastaList[1:5])
    print('\nThe fastaList for "%s" with total entries %s has been prepared\n' % (fastaFile, acount))
    print('**There were %s entries found with empty sequences and were removed**\n' % (empty_count))
    
    return fastaList
    
def batchPlot(idFile,precursorList,phasiList):
    '''Takes a list of FASTA ids, fetched precursor from phased transcripts fasta file and
    then extracts phasiRNAs from phasiExtract output and feed to RNA plot'''
    fh_in = open(idFile,'r')
    idRead = fh_in.readlines()

    for i in idRead:
        query = i.strip("\n")
        
        ## Fetch precursor and phasiRNAs for this ID in temp files
        LONG_SEQ_FILE,LONG_SEQ          = fetchPrecursor(query,precursorList)
        SHORT_SEQ_FILE,SHORT_SEQ_FASTA  = fetchPhasi(query,phasiList)
        RNAFOLD(query,LONG_SEQ,SHORT_SEQ_FILE)
        # sys.exit()

    return None

def RNAFOLD(query,LONG_SEQ,SHORT_SEQ_FILE):
    '''Uses external scrip tto draw plot of provided precursor and phasiRNAs'''

    NAME = query.replace("|","_")
    print("Folding and plotting phasiRNAs for query:%s" % (NAME))

    retcode2 = subprocess.call(["perl","HPAnnotate.v01.pl",str(LONG_SEQ),str(SHORT_SEQ_FILE),str(NAME)]) # RUN PERL SCRIPT
    time.sleep(2)

    if retcode2 != 0:
        print("There was a problem in generating plots -Investigate !!")
    else:
        print("Folding and plotting done for:%s\n" % (NAME))
        pass

    return None

def tagExtracter(SHORT_SEQ_FILE):

    '''DEPRECATED-Extracts tags from FASTA file and provide input to script'''

    fh_in = open(SHORT_SEQ_FILE,'r')
    fasta = fh_in.read()
    fasta_splt = fasta.split('>')

    outfile = "%s.tags" % (SHORT_SEQ_FILE)
    fh_out = open(outfile,'w')

    for i in fasta_splt[1:]:
        # print("This is tag entry:",i)
        ent = i.split('\n')
        name = ent[0].split()[0].strip()
        seq = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        seq_rc = seq[::-1].translate(str.maketrans("TAGC","ATCG"))
        seq_c = seq.translate(str.maketrans("TAGC","ATCG"))
        if mode == 0:
            fh_out.write("%s\n" % (seq))
        elif mode == 1: 
            fh_out.write("%s\n%s\n" % (seq,seq_rc))
        elif mod == 2:
            fh_out.write("%s\n%s\n%s\n" % (seq,seq_rc,seq_c))

    fh_in.close()
    fh_out.close()

    return outfile

def fetchPrecursor(anid,precursorList):

    ''' Takes input name and fetched fasta entry for transcript'''
    
    print("Extracting precursor %s" % (anid))

    LONG_SEQ_FILE = '%s.long.temp' % (anid)
    fh_out = open(LONG_SEQ_FILE,'w')

    acount = 0 ## Count of matched
    resList = [] ## In case there are more then one precursors ony first precurspr will be used
    for i in precursorList:
        name,seq = i
        # print("Matching to ",name)

        if name.startswith(anid):
            print("Precursor %s match found:%s" % (anid,name))
            fh_out.write(">%s\n%s\n" % (name,seq.strip("\n")))
            resList.append((name,seq))
            acount+=1

    if len(resList) > 1:
        print("More then one precursor transcript with %s name found" % (anid))
        sys.exit()
    
    elif len(resList) == 0:
        print("No match for this transcript:%s found in precursors list" % (anid))
        sys.exit()
    
    if len(resList) == 1:
        pass

    fh_out.close()

    LONG_SEQ = resList[0][1]
    print(resList)
    # sys.exit()

    return LONG_SEQ_FILE,LONG_SEQ

def fetchPhasi(anid,phasiList):

    '''Takes input name and fetches phasiRNAs from phasiExtract results file'''
    print("Extracting phasiRNAs %s" % (anid))

    SHORT_SEQ_FASTA = '%s.phasi.temp' % (anid)
    SHORT_SEQ_FILE = '%s.tags.temp' % (anid)
    fh_out1 = open(SHORT_SEQ_FASTA,'w')
    fh_out2 = open(SHORT_SEQ_FILE,'w')

    acount = 0 ## Count of matched
    for i in phasiList:
        nameLine,seqLine = i
        name = nameLine.split(",")[0].strip("\n")
        seq = seqLine.split(",")[0].strip("\n")
        # print("Matching to ",name)

        if name.startswith(anid):
            print("Precursor %s match found:%s" % (anid,name))
            fh_out1.write(">%s\n%s\n" % (name,seq))
            acount +=1

            if mode == 0:
                fh_out2.write("%s\n" % (seq))
            elif mode == 1: 
                fh_out2.write("%s\n%s\n" % (seq,seq_rc))
            elif mod == 2:
                fh_out2.write("%s\n%s\n%s\n" % (seq,seq_rc,seq_c))
            else:
                print("Please input correct 'mode' input")
                sys.exit()

    if acount > 0:
        pass
    else:
        print("No match for this transcript:%s found in precursors list" % (anid))

    fh_out1.close()
    fh_out2.close()

    return SHORT_SEQ_FILE,SHORT_SEQ_FASTA

def singlePlot(LONG_SEQ_FILE,SHORT_SEQ_FILE):
    ''' Extract precursor and name, and tags from SHORT_SEQ file. Finally give these 
    to RNAFOLD to fold precursor and annotate phasiRNAs - This module is part of manual mode
    '''

    fh_in = open(LONG_SEQ_FILE,'r')
    longList = fh_in.read().split('>')
    query,seq,trash = longList[1].split("\n")
    print(query,seq)

    SHORT_SEQ_TAGS = tagExtracter(SHORT_SEQ_FILE)
    RNAFOLD(query,seq,SHORT_SEQ_TAGS)

    return None

def cleanTemp():
    garbage = [afile for afile in os.listdir('./') if afile.endswith (('.temp'))] ## Excluded-'chopped.trimmed.fastq' as used by RNA Runner
    
    for afile in garbage:
        if os.path.isfile(afile): ## Check to see its a file from bowtie and not tophat mapped folder - Untested
            print("Deleting %s" % (afile))
            os.remove(afile)
        else:
            print("Skiping cleanup, as its a directory %s" % (afile))

def main():

    if runMode == 'M':
        LONG_SEQ_FILE   = str(sys.argv[1])  #"Aspa_mapped_genomic_seq.fa" 
        SHORT_SEQ_FILE  = str(sys.argv[2])  #"SHORT.txt"
        singlePlot(LONG_SEQ_FILE,SHORT_SEQ_FILE)

    elif runMode == 'A':
        print("Caching precursors and phasiRNAs file")
        precursorList   = fastaReader(fastaFile)
        phasiList       = fastaReader(phasiFile)
        batchPlot(idFile,precursorList,phasiList)    
        ## Delete temp files
        if cleanup == 1:
            cleanTemp()

    else:
        print("Please input correct mode - Script will exit")
        sys.exit()

if __name__ == '__main__':
    main()
    print ("\nDone..\n")
    sys.exit()


### LOG
### v01 -> v02
### Added a manual mode in which user can give precursor file and phasiRNA tags as an command line input rather then a list of ids as in automatic mode





### README ###############################
##########################################

### MODE: AUTOMATIC
### Provide a fasta file of transcripts/precursors from where seqeunces will be extracted for ids in idFIle
### Provide phasiRNAs as extracted from phasiExtract from which phasiRNAs correspoding to precursors will be extracted
### Script will generate plots for all ids in idList w/o any intervenation

### MODE: MANUAL
### Provide a fasta file of one single precursor as first command line argument
### Provide a fasta file of phasiRNAs as second command line argument
### script will generate plot for single precursor
### This mode is required for inffered IRs for which precursor is made by combining two transcripts

        
    
    
    
