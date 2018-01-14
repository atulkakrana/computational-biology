#!/usr/local/bin/python3

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
####Please feed the target finder file either combined or for individual miRNAs through command line
####Usage: ./TarFindParse_v2.py NAME_OF_FILE
###Result: A file with TF_parsed.csv at end is the resulting file
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# - - - - - - F U N C T I O N S- - - - - - - - - - - - - - - 

import csv
import sys
import re
import os

# - - - - - U S E R    V A R I A B L E S - - - - - - - -
##Config
#cleave_file = 'comb_res_target_finder'
cleave_file = str(sys.argv[1])####To take arguments from command line

# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -


##First clean the results file - option added on 7-Feb-2013

TF_res_clean = './TF_results_all_cleaned'
fh_out2 = open(TF_res_clean, 'w') ###Intermidiate file to write clean results

fh_in2=open(cleave_file, 'r')##Can give errors due to file opening and closing
No_res_re = re.compile('No results')
for i in fh_in2:
    if re.search(No_res_re, i): ##checks for file saying : No results for ath420
        pass
    else:
        fh_out2.write(i)
fh_in2.close()
fh_out2.close()


##Open a files to be written in CSV format
out_file = ('%s_TF_parsed.csv' % (cleave_file))
fh_output=open( out_file, 'w') # File to compare entries with validated ones
#csv_out=csv.writer(fh_output, delimiter=',')

##Write header to outputfile
fh_output.write('miRNA name,target_gene,target_location,mirna_seq,tar_seq,targetfinder_score\n')

#fh_output2=open('/home/atul/Analysis/unvalidated_scoring_input', 'w')# File to score directly at this step with comparing with validated ones
#csv_out2=csv.writer(fh_output2, delimiter=',')

#Open the input file for operations($$$$$$$$$$$$)
#Please enter the NAME OF CLEAVELAND FILE to be parsed here
fh_in=open(TF_res_clean, 'r')
csv_in=csv.reader(fh_in, delimiter='\t')
csv_table=[]

##Populate the list
for row in csv_in:
	csv_table.append(row)
	
#print(csv_table[0:12])
	
#print(csv_table[0:12])##+6 to get next block i.e 6 lines makes an entry
entries=int(len(csv_table)/6)
print('Total number of entries in input file:',entries,'\n')


#Extract Entry, miRNA, Target
result_entries=0
a_count=0###incremented to move the next block
b_count=6
for i in range(entries):
	ent_loc=csv_table[a_count:b_count]##A miR and target entry complete
#	print(ent_loc)
	info_block=ent_loc[0][0].split(',')##block containing names and other info
#	print(info_block)
	mirna=info_block[0].split('=')[1]
	target=info_block[1].split('=')[1]
	score=info_block[2].split('=')[1]
	range=info_block[3].split('=')[1]
	
	tar_seq=ent_loc[2][0].split()[2]
	mir_seq=ent_loc[4][0].split()[2]
		
	
#	print(mirna,target,score,range, mir_seq, tar_seq)
##	break

	fh_output.write('%s,%s,%s,%s,%s,%s\n' % (mirna,target,range,mir_seq,tar_seq,score))
#	csv_out.writerow([mirna]+[target]+[range]+[mir_seq]+[tar_seq]+[score])
	a_count+=6
	b_count+=6
	result_entries+=1
	
print('Number of entries in Target finder parsed file:', result_entries)
print('The output file with parsed results generated is "%s"' % (out_file))
	
# Close output file: parsed_out
os.remove('./TF_results_all_cleaned')##Delete intermidiate file
fh_output.close()
fh_in.close()





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
