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
threshold = 500     ## Default:50
maxRepLen = 3000    ## Default: 2000
parallel = 1 		## Default: 1 (Yes)

### Steps #####
fastaList_step = 1
einverted_step = 1
resParse_step  = 1
combClean_step = 1

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
        name = anent[0].replace(":","_")
        seq = anent[1]
        
        tempInput = "tempSeq.fa"
        fh_out = open(tempInput,'w')
        fh_out.write('>%s\n%s\n' % (name,seq))
        fh_out.close()

        outseq = "%s.fa.temp" % name
        outinv = "%s.inv.temp" % name

        retcode = subprocess.call(["einverted", "-sequence", tempInput, "-gap", str(gap), "-threshold", str(threshold), "-match",str(match),"-mismatch", str(mismatch), "-maxrepeat",str(maxRepLen), "-outfile",outinv, "-outseq",outseq ])

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

    name = atup[0].replace(":","_")
    seq = atup[1].strip()
    
    print("Seqeunce being analyzed for inverted repeats: %s\n" % (name))
    tempInput = ("%s_tempSeq.fa" % (name))
    fh_out = open(tempInput,'w')
    fh_out.write('>%s\n%s' % (name,seq))
    fh_out.close()

    outseq = "%s.fa.temp" % name
    outinv = "%s.inv.temp" % name

    retcode = subprocess.call(["einverted", "-sequence", tempInput, "-gap", str(gap), "-threshold", str(threshold), "-match",str(match),"-mismatch", str(mismatch),"-maxrepeat",str(maxRepLen), "-outfile",outinv, "-outseq",outseq ])

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
    invRes = [file for file in os.listdir('./') if file.endswith ('.inv.temp')]

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

                fh_out.write("%s,%s,%s,%s" % (afile.split(".")[0],score,matches,gaps)) ## name of file,score,match,gaps

            
            else:
                # print ("Sequence %s - No results" % (afile.split(".")[0]))
                pass
        except OSError:
            print ("No result file for sequence %s found - Please check" % (afile))
            print("System will exit")
            sys.exit()

    fh_out.close()

    return resOut

def parseRes2():

    print ("Parsing results\n")
    ## File to record results
    resOut= "RES_%s_%s.csv" % (fastaFile,datetime.datetime.now().strftime("%m_%d_%H_%M"))
    fh_out = open(resOut,'w')
    fh_out.write("EntryName,Score,Matches,Perc,Gaps,AlignLen,5'Start,5'End,3'start,3'end,Loop\n")

    ## Get list of files
    invRes = [file for file in os.listdir('./') if file.endswith ('.inv.temp')]
    for afile in invRes:
        print("\nSeqeunce being parsed:%s" % (afile))
        try:
            if os.stat(afile).st_size > 0:
                print ("Sequence %s - Results !!!!" % (afile.split(".")[0]))
                fh_in = open(afile,'r')
                invs = fh_in.read().split("\n\n")
                # print("Empty line splitted:",invs)
                
                for i in invs:
                    # print('\nInverted:',i.strip('\n'))
                    invLines = i.strip('\n').split('\n')
                    # print(invLines)
                    
                    resBlock_splt = invLines[0].split(":")
                    # print(resBlock_splt)
                    score = resBlock_splt[1].split()[1]
                    # print(resBlock_splt[1].split())
                    matchesInfo,gapsInfo = resBlock_splt[2].split(",")
                    # print(matchesInfo,gapsInfo)
                    gaps = gapsInfo.split()[0]
                    
                    # print(matchesInfo.strip().split())

                    matches = matchesInfo.strip().split()[0]
                    matched,total = matches.split("/")
                    alignLen = int(total)*2
                    perc = round(int(matched)/int(total),2)

                    # matches,garbage1,perc,garbage2 = matchesInfo.strip().split()
                    # alignLen = int(matches.split("/")[1])*2
                    # print(score,matches,gaps)

                    arm5 = invLines[1].strip()
                    start5,seq5,end5 = arm5.split(" ")

                    arm3 = invLines[3].strip()
                    end3,seq3,start3 = arm3.split(" ")

                    loop = int(start3)-int(end5)

                    fh_out.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (afile.split(".")[0],score,matches,str(perc),gaps,str(alignLen),start5,end5,start3,end3,loop)) ## name of file,score,match,gaps
            else:
                # print ("Sequence %s - No results" % (afile.split(".")[0]))
                pass

        except OSError:
            print ("No result file for sequence %s found - Please check" % (afile))
            print("System will exit")
            sys.exit()

    return resOut

def combineClean():

    print("Combining result file\n")

    fastaFiles = [file for file in os.listdir('./') if file.endswith ('.fa.temp')]
    invFiles = [file for file in os.listdir('./') if file.endswith ('.inv.temp')]
	
    ## Combine inverted repeats alingments results
    invComb = './Inverted_Alings.inv'
    inv_out = open(invComb ,'w')
    print("\nCombining 'inverted' results file")
    for x in invFiles:
        print (x)
        invfile = open('./%s' % (x), 'r')
        #targfile.readline()
        data = invfile.read()
        invfile.close()
        inv_out.write(data)
    inv_out.close()

    ## Combine inverted repeats sequences FASTA file
    fastaComb = './Inverted_Seqs.fa'
    fasta_out = open(fastaComb ,'w')
    print("\nCombining 'inverted' fasta file")
    for x in fastaFiles:
        print (x)
        fastafile = open('./%s' % (x), 'r')
        #targfile.readline()
        data2 = fastafile.read()
        fastafile.close()
        fasta_out.write(data2)
    fasta_out.close()

    # Delete individual inverted files and fasta files
    print("Combined files '%s and %s' ready - Deleting individual files\n" % (invComb,fastaComb))
    garbage = fastaFiles + invFiles
    for file in garbage:
        print("Deleting %s" % (file))
        os.remove(file)

    # Delete one temp seq file
        
    return invComb,fastaComb

def PP(module,alist):
    print('***********Parallel instance of %s is being executed*********' % (module))
    
    start = time.time()
    nprocPP = round((accel/int(nspread))+1) #
    print('\nnprocPP:%s\n' % (nprocPP))
    npool = Pool(int(nprocPP))
    npool.map(module, alist)

def main():

    if fastaList_step == 1:
        fastaList = readFASTA(fastaFile)
    else:
        print("No FASTA file is read as the required step is turned-off\n")

    if einverted_step == 1:
        if parallel == 0:
            print("Running einverted repeat in non-parallel mode\n")
            batcheinverted(fastaList)
        elif parallel == 1:
            print("Running einverted repeat in parallel mode\n")
            PP(einverted,fastaList)
        else:
            print("Please input correct value for 'parallel' variable' - System will exit now\n")
            sys.exit()
    else:
        print("No inverted repeat analysis is performed as the step is turned off\n")
    
    if resParse_step == 1:
        resOut = parseRes2()
    else:
        print("Results are not combined as the step is turned off\n")

    if combClean_step == 1:
        invComb,fastaComb = combineClean()
    else:
        print("Results will not be combined from individual seqeunces and no cleanup will be performed, as the step is off\n")
        print("This might generate thousands of file leaving your system unstable if you tried to open the results folder\n")

if __name__ == '__main__':
    nspread = 1
    if parallel == 1:
        accel = int(multiprocessing.cpu_count()*0.80)
    else:
        pass
    main()
    sys.exit()


## v01
## Script to run EMBOSS utility 'einverted' in batch mode. Using the 'einverted' directly on set of sequences gives wrong results

## V01->v02
## Added parallel processing for large number of queries or scaffolds
## Added a function to combine result files and delete individual files

##v02 -> v03
## Added parser to parse multiple alingments in one sequence

##v02 -> 04
## Fixed the error that "File cannot be created" in parallel mode
