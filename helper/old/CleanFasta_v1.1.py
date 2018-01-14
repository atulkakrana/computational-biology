#!/usr/local/bin/python3

##This module cleans the FASTA file header
##Trunchates the header till first whitespace and pipe, you can change the nature below

##Usage: python3 clean_fasta_header_v1.py ENTER_FASTA_FILE


import re
import sys

##config 



def FASTAClean(filename):
    fh_in=open(filename, 'r')
    #write file
    out_file = ('%s_Clean.fa' % (filename.split('.')[0]))
    fh_out =open(out_file, 'w')
    
    print ('\nCleaning "%s" FASTA file\n' % (filename))
    
    fasta = fh_in.read()
    fasta_splt = fasta.split('>')
    acount = 0 ## count the number of entries
    empty_count = 0
    for i in fasta_splt[1:]:
        ent = i.split('\n')
        name = ent[0].split()[0].strip()
        seq = ''.join(ent[1:])##Sequence in multiple lines
        if seq:
            fh_out.write('>%s\n%s\n' % (name,seq))
            acount+=1
        else:
            empty_count+=1
            pass          
        acount+=1
    
    fh_in.close()
    fh_out.close()  

    print('\nThe fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file, acount))
    print('**There were %s entries found with empty sequences and were removed**' % (empty_count))
    return out_file

def main():
    #inp_file_name = 'medicago.cdna'
    inp_file_name = str(sys.argv[1])
    #inp_file_name = input('Please enter the file name:')
    FASTAClean(inp_file_name)
    
    
if __name__ == '__main__':

    main()
    sys.exit()


'''
Created on Apr 23, 2012

@author: atul
'''
##1 -> 1.1
##Added function to remove gene entries with empty sequences
##Changed the output file name by removing '.fa' in between of name

