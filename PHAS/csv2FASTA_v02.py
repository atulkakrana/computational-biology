#!/usr/local/bin/python3

import sys


##Signature
print("\nScript by /-\ ^|^ |_| |_ for Blake's Lab\n")


##Config
file_in = str(sys.argv[1])
#file_in = 'Mt3.5_non_TIR_jixian_protein.csv'
name_row = 0
seq_row = 1
##Chnage separater type in line22


##Open file
fh_in= open(file_in,'r')

##Write FASTA format to another file
out_file = file_in.split('.')[0]+'.fa'
fh_out=open(out_file,'w')

##Extract name and sequence
for row in fh_in:
#    print (row)
    entry_strpd=row.strip('\n')
    entry=entry_strpd.split(',')
    name=entry[name_row]
    seq=entry[seq_row]
    seq=seq.upper()
#    print(name, seq)
#    print(">%s\n%s\n" % (name,seq))
    fh_out.write(">%s\n%s\n" % (name,seq))
    
fh_in.close
fh_out.close

