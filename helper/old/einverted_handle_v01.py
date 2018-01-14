#!/usr/local/bin/python3
#### Script to get sequences from genome DB - written by kakrana@udel.com

######## IMPORT ################
import sys,os,subprocess,multiprocessing,datetime


######## SETTINGS ###############
fastaFile = 'Final_24PHAS_Loci_ALL_v4_dynamicbuff.fa'
# fastaFile = 'test.fa'
gap = 12            ## Default:12 (from emboss website)
match = 3           ## Default:3
mismatch = -4       ## Default:-4
threshold = 500     ## Default:50


############ FUNCTIONS ###########

def readFASTA(afile):

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

def parseRes():

    '''Reads the result file for each sequence and generated a summary file 
    with entries for sequences (and files) that have results passing provided threshold'''

    ## File to record results
    resOut= "RES_%s_%s" % (fastaFile,datetime.datetime.now().strftime("%m_%d_%H_%M"))
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



def main():
    fastaList = readFASTA(fastaFile)
    batcheinverted(fastaList)
    resOut = parseRes()
    pass

if __name__ == '__main__':
    main()
    sys.exit()



## Script to run EMBOSS utility 'einverted' in batch mode. Using the 'einverted' directly on set of sequences gives wrong results
##