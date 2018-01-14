#!/usr/local/bin/python3
#### Script to get sequences from genome DB - written by kakrana@udel.com

######## IMPORT ################
import sys,os,subprocess,multiprocessing,datetime,time
from multiprocessing import Process, Queue, Pool


######## SETTINGS ###############

fastaFile = sys.argv[1]
# fastaFile = 'Final_24PHAS_Loci_ALL_dynamicbuff3_v4.fa'
# fastaFile = 'test.fa'
gap = 12            ## Default:12 (from emboss website)
match = 3           ## Default:3
mismatch = -4       ## Default:-4
threshold = 300     ## Default:50
parallel = 1 		## Default: 1 (Yes)


############ FUNCTIONS ###########

def readFASTA(afile):
    '''Read Fasta file and generate a list of sequences for 
    parallel processing '''

    print("Generating list of fasta seqeunces\n")

    fh_in = open(afile,'r')
    fasta = fh_in.read()
    fasta_splt = fasta.split('>')
    
    fastaList = [] ## Store entries to be processed sequentially
    acount = 0 ## count the number of entries
    empty_count = 0
    for i in fasta_splt[1:]:
        ent = i.split('\n')
        name = ent[0].split()[0].strip()
        seq = ''.join(ent[1:])##Sequence in multiple lines
        if seq:
            fastaList.append((name,seq))
            # fh_out.write('>%s\n%s\n' % (name,seq))
            acount+=1
        else:
            empty_count+=1
            pass          
            acount+=1

    fh_in.close()
    print("Total entries in File: %s | Total empty entries in file: %s" % (acount,empty_count))
    print("Total entries in fastList: %s\n" % (str(len(fastaList))) )

    return fastaList

def batcheinverted(fastaList):
    ''' Running einverted on a FASTA file with number of sequences 
    directly through einverted gives wrong results.So, this script runs einverted one by one on each sequence in FASTA file. 
    '''

    for anent in fastaList:
        name = anent[0]
        seq = anent[1]
        
        tempInput = "tempSeq.fa"
        fh_out = open(tempInput,'w')
        fh_out.write('>%s\n%s\n' % (name,seq))
        fh_out.close()

        outseq = "%s.inv.fasta" % name
        outinv = "%s.inv" % name

        retcode = subprocess.call(["einverted", "-sequence", tempInput, "-gap", str(gap), "-threshold", str(threshold), "-match",str(match),"-mismatch", str(mismatch), "-outfile",outinv, "-outseq",outseq ])

        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('\n****einverted for %s complete****\n' % (name) )
        else:
            print('Something wrong happened while running einverted for sequence: %s - - Debug for reason' % (name))
            sys.exit()

        ## Cleanup entry specifc FASTA file
        if os.path.exists(tempInput):
            os.remove(tempInput)

    pass

def einverted(atup):

    name = atup[0]
    seq = atup[1]
    
    tempInput = "%s_tempSeq.fa" % name
    fh_out = open(tempInput,'w')
    fh_out.write('>%s\n%s\n' % (name,seq))
    fh_out.close()

    outseq = "%s.inv.fasta" % name
    outinv = "%s.inv" % name

    retcode = subprocess.call(["einverted", "-sequence", tempInput, "-gap", str(gap), "-threshold", str(threshold), "-match",str(match),"-mismatch", str(mismatch), "-outfile",outinv, "-outseq",outseq ])

    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        print('\n****einverted for %s complete****\n' % (name) )
    else:
        print('Something wrong happened while running einverted for sequence: %s - - Debug for reason' % (name))
        sys.exit()

    ## Cleanup entry specifc FASTA file
    if os.path.exists(tempInput):
        os.remove(tempInput)

    pass

def parseRes():

    '''Reads the result file for each sequence and generated a summary file 
    with entries for sequences (and files) that have results passing provided threshold'''

    print ("Parsing results\n")
    ## File to record results
    resOut= "RES_%s_%s.csv" % (fastaFile,datetime.datetime.now().strftime("%m_%d_%H_%M"))
    fh_out = open(resOut,'w')
    fh_out.write("EntryName,Score,Matches,Mismatches\n")

    ## Get list of files
    invRes = [file for file in os.listdir('./') if file.endswith ('.inv')]

    for afile in invRes:
        try:
            if os.stat(afile).st_size > 0:
                print ("Sequence %s - Results !!!!" % (afile.split(".")[0]))
                fh_in = open(afile,'r')
                fileRead = fh_in.readlines()
                resBlock = fileRead[1]
                print(resBlock)

                resBlock_splt = resBlock.split(":")
                score = resBlock_splt[1].split()[1]

                # print(resBlock_splt[1].split())
                matches,gaps = resBlock_splt[2].split(",")

                fh_out.write("%s,%s,%s,%s" % (afile.split(".")[0],score,matches,gaps))

            
            else:
                # print ("Sequence %s - No results" % (afile.split(".")[0]))
                pass
        except OSError:
            print ("No result file for sequence %s found - Please check" % (afile))
            print("System will exit")
            sys.exit()

    fh_out.close()

    return resOut

def combineClean():

    print("Combining result file\n")

    fastaFiles = [file for file in os.listdir('./') if file.endswith ('.inv.fasta')]
    invFiles = [file for file in os.listdir('./') if file.endswith ('.inv')]
	
    ## Combine inverted repeats alingments results
    invComb = './All.inv'
    inv_out = open(invComb ,'w')
    print("Combining 'inverted' results file\n")
    for x in invFiles:
        print (x)
        invfile = open('./%s' % (x), 'r')
        #targfile.readline()
        data = invfile.read()
        invfile.close()
        inv_out.write(data)
    inv_out.close()

    ## Combine inverted repeats sequences FASTA file
    fastaComb = './All.inv.fa'
    fasta_out = open(invComb ,'w')
    print("Combining 'inverted' fasta file\n")
    for x in fastaFiles:
        print (x)
        fastafile = open('./%s' % (x), 'r')
        #targfile.readline()
        data = fastafile.read()
        fastaile.close()
        fasta_out.write(data)
    fasta_out.close()

    ## Delete individual inverted files and fasta files
    print("COmbined files ready - Deleting individual files\n")
    garbage = fastaFiles + invFiles
    for file in garbage:
        print("Deleting %s" % (file))
        os.remove(file)
        
    return invComb,fastaComb

def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    nprocPP = round((accel/int(nspread))+1) #
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

def main():
    fastaList = readFASTA(fastaFile)

    if parallel == 0:
        print("Running einverted repeat in non-parallel mode\n")
        batcheinverted(fastaList)
    elif parallel == 1:
        print("Running einverted repeat in parallel mode\n")
        PP(einverted,fastaList)
    else:
        print("Please input correct value for 'parallel' variable' - System will exit now")
        sys.exit()
    
    resOut = parseRes()
    invComb,fastaComb = combineClean()

if __name__ == '__main__':
    nspread = 1
    if parallel == 1:
        accel = int(multiprocessing.cpu_count()*0.80)
    else:
        accel = int(accel)
    main()
    sys.exit()


## v01
## Script to run EMBOSS utility 'einverted' in batch mode. Using the 'einverted' directly on set of sequences gives wrong results

## V01->v02
## Added parallel processing for large number of queries or scaffolds
## Added a function to combine result files and delete individual files