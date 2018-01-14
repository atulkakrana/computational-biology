##This module cleans the FASTA file header
##Trunchates the header till first whitespace and pipe, you can change the nature below

##Usage: python3 clean_fasta_header_v1.py ENTER_FASTA_FILE


import re
import sys

##config 

#inp_file_name = 'medicago.cdna'
inp_file_name = str(sys.argv[1])
#inp_file_name = input('Please enter the file name:')


#read file
fh_in=open(inp_file_name, 'r')
#write file
out_file = ('%s_new_head.fa' % (inp_file_name))
fh_out =open(out_file, 'w')

acount = 0 ## count the number of entries
for i in fh_in:
    if re.match('>', i):
        header = i.split()##change here if separater is not white space
        new_head = header[0].split('|')[0]## Change here for number of fields you want strating from 0
        fh_out.write('%s\n' % new_head)
        acount+=1
#        print(i)
#        print(new_head)
    else:
        fh_out.write('%s' % i)
        
fh_in.close()
fh_out.close()

print('\nThe fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file, acount))
    







'''
Created on Apr 23, 2012

@author: atul
'''

