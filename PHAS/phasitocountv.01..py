#!/usr/local/bin/python3

### Script takes phasiExtract results, extracts abundance from header and write in unique tags in tag count format
### script might need modify for new phasiExtract folrmat that includes p-value and hits

import sys,os

# infile = "Hl_21Phasi.fa"
infile  = str(sys.argv[1])
abunPos = 8 ## In python format could be a positive or negative number | default = -1


def phasitocount(infile):
    '''
    Read file and extract tag and abundance 
    '''

    #### I/O ###############################################
    #######################################################
    tagcountF   = "%s.txt" % infile.rpartition(".")[0]
    fh_out      = open(tagcountF, 'w')

    fh_in = open(infile,'r')
    fasta = fh_in.read()
    fasta_splt = fasta.split('>')

    ##### Read fasta file and convert to tag count ########
    #######################################################
    acount      = 0 ## count the number of entries
    empty_count = 0
    alist       = []
    for i in fasta_splt[1:]:

        ent     = i.split('\n')
        name    = ent[0].split()[0].strip()
        abun    = name.split("_")[abunPos]
        seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        alen    = len(seq)
        alist.append((seq,abun))
    print("This is snippet of tag list",alist[1:5])

    ### Sort list on abundnace
    blist = sorted(alist,key=lambda x: int(x[1]),reverse=True) ## Sorted on bitscore
    print("This is snippet of abun sorted tag list",blist[1:5])

    ### Now remove redundant ###############################
    #######################################################
    uniqL = [] ### Set to store unique tags
    for i in blist:
        atag,aabun = i
        if atag not in uniqL:
            uniqL.append(atag)
            fh_out.write("%s\t%s\n" % (atag,aabun))
        else:
            # print("Tag Already recorded")
            pass

    fh_in.close()
    fh_out.close()


    return tagcountF

def main():
    tagcountF = phasitocount(infile)


if __name__ == '__main__':
    main()


#####