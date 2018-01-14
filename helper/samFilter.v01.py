#!/usr/local/bin/python3
#### Script to filter SAM files - written by kakrana@udel.com

######## Imports ############
import os,sys

######### Settings ############

samFile     = 'map_out.sam'
includeFlag = [256,]        ## 0 i.e. no flag is included by default
excludeFlag = []            ## Used when you want to exclude specific flags
minQual     = 255           ## Default should be low cutoff 30 and high cutoff 255
unMap       = 1             ## 0: Include them in SAM FILE 1: Seprate out in Other SAM file


####### Function #############

def filterSAM(samFile,includeFlag,excludeFlag,minQual):

    fileOut = "%s.filter.sam" % samFile
    fh_out  = open(fileOut,'w')

    fh_in =open(samFile,'r')
    samRead = fh_in.readlines()

    if unMap == 1:
        fileOut2 = "%s.unmap.sam" % samFile
        fh_out2 = open(fileOut2,'w')

    acount = 0      ## All entries
    bcount = 0      ## passed entries
    xcount = 0      ## Entries with exclude flags
    for i in samRead:
        x = i.strip("\n")
        # print ("sam entry:%s " % x)


        if x.startswith('@'):
            # print("Header: %s" % (x))
            fh_out.write("%s\n" % x)
        else:
            print("Mapped entry: %s" % (x))
            acount += 1 
            
            ent     = x.split("\t")
            flag    = int(ent[1])
            qual    = int(ent [4])

            if includeFlag:
                
                ## Matches our flags
                if (flag in includeFlag) and (flag not in excludeFlag):
                    if int(qual) >= minQual:
                        fh_out.write("%s\n" % x)
                        bcount +=1
                ## No flag set - default mapped reads
                elif flag == 0 and (flag not in excludeFlag):
                    if int(qual) >= minQual:
                        fh_out.write("%s\n" % x)
                        bcount +=1

                else:
                    pass


            else: ## includeFlag is empty
                if (flag not in excludeFlag):
                    if qual >= minQual:
                        fh_out.write("%s\n" % x)
                        bcount +=1

            if unMap == 1: ## Filter out unmapped flags in separet SAM file
                if flag == 4:
                        fh_out2.write("%s\n" % x)
                        xcount +=1



    print ("\nTotal mapped entries: %s | Passed entries: %s | Unmapped: %s (if unMap enabled) )" % (acount,bcount,xcount))
    print ("Results saved in %s file\n" % (fileOut))
    fh_in.close()
    fh_out.close()

    if unMap == 1:
        fh_out2.close()

    return fileOut


def main():
    resFile = filterSAM(samFile,includeFlag,excludeFlag,minQual)


if __name__ == '__main__':
    main()
    sys.exit()

