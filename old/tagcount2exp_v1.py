##This script is to convert tag_count files back to original expression fasta file by adding tag number of times in tag count file
##Usage : python3 tagcount2exp_v1.py <input>


import re
import sys

#inp_file_name = 'medicago.cdna'
inp_file_name = str(sys.argv[1])

#read file
fh_in=open(inp_file_name, 'r')
#write file
out_file = ('%s_expression.fa' % (inp_file_name))
fh_out =open(out_file, 'w')

##Write to file
print('Writing expression file for %s tagcount file' % (inp_file_name))
print('\n~~PLEASE BE PATIENT~~~')
#fh_out = open('./PARE/%s_PARE_tags.fa' % (lib[0]), 'w')##Naming file with lib_ids name
tag_num = 1
for ent in fh_in:##All the entries of the library
    #if len(ent[0]) == 20:
    ent = ent.split('\t')
    tag_count = int(ent[1])
    for count in range(tag_count):##Number of times the tag_count file
        fh_out.write('>%s\n%s\n' % (tag_num, ent[0]))
        tag_num += 1
  
  
print('\nThe expression file from %s is ready to use' % (inp_file_name))
fh_in.close()
fh_out.close()

