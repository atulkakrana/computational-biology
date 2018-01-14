#!/usr/local/bin/python3

import sys

## Script converts comma and tab seprated file to FASTA format

##Signature
print("\nScript by /-\ ^|^ |_| |_ for Blake's Lab\n")

##Config
file_in     = str(sys.argv[1])
#file_in = 'Mt3.5_non_TIR_jixian_protein.csv'

mode        = str(sys.argv[2]) ## 0: for csv2FASTA 1: for FASTA2csv

## for csv2FASTA
name_row    = 1 ## Excel format
seq_row     = 2 ## Excel format
head        = 'Y' ## If header first line of csv file will be skipped

## for FASTA2csv and csv2FASTA both
sep = '\t'

##Change separater type in line22

def csv2FASTA(file_in,sep,head):
    
    ## Open file
    fh_in= open(file_in,'r')
    if head == 'Y':
        fh_in.readline()
    ## Write FASTA format to another file
    file_out    = "%s.fa" % (file_in.rpartition(".")[0])
    fh_out      = open(file_out,'w')
    
    ## Extract name and sequence
    for row in fh_in:
        print ("Ent:",row)
        entry_strpd =row.strip('\n')
        entry       =entry_strpd.split(sep)
        name        =entry[name_row-1]
        seq         =entry[seq_row-1]
        seq         =seq.upper()
    #    print(name, seq)
    #    print(">%s\n%s\n" % (name,seq))
        fh_out.write(">%s\n%s\n" % (name,seq))
    
    fh_in.close
    fh_out.close
    
    return file_out

def FASTA2csv(file_in,sep):
    
    ## Open file
    fh_in       = open(file_in,'r')
    ## Write FASTA format to another file
    file_out    = "%s.mod.csv" % (file_in.rpartition(".")[0])
    fh_out      = open(file_out,'w')
    
    file_data   = fh_in.read()
    
    ## Extract name and sequence
    for entry in file_data.split('>')[1:]: ## First entry is empty
    #    print (row)
        entry_splt  = entry.split('\n')
        name        = entry_splt[0] ## '>' not included
        seq         = ''.join(x.strip() for x in entry_splt[1:]) ## Sequence in multiple lines
        #seq=seq.upper()
    #    print(name, seq)
    #    print(">%s\n%s\n" % (name,seq))
        fh_out.write("%s%s%s\n" % (name,sep,seq))
    
    fh_in.close
    fh_out.close
    
    return file_out

def main(file_in):
    if mode == '0':
        csv2FASTA(file_in,sep,head)
        print('%s converted to FASTA' % (file_in))
    elif mode == '1':
        FASTA2csv(file_in,sep)
        print('%s converted to CSV' % (file_in))
    else:
        print('Please enter correct mode')
        print('You entered: %s' % (mode))
        sys.exit


if __name__ == "__main__":
    main(file_in)
    sys.exit()


## V2 -> v3
## Added FASTA 2 csv
## Added header removal

## v3.1 -> v3.2
## CLeaned up code
## improved output name
## input rows are in excel format

    


