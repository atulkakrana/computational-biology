# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

import sys




# - - - - - U S E R    V A R I A B L E S - - - - - - - -
##Config

cleave_file = str(sys.argv[1])
#cleave_file = 'gma_mirBASE_all.fa_soy_full_genome_REAL_targets_7'


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -



# - - - - - - F U N C T I O N S- - - - - - - - - - - - - - - 

import csv

# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -

##Count the number of entries
#Please enter the NAME OF CLEAVELAND FILE to be parsed here
fh_in=open(cleave_file, 'r')
csv_input=fh_in.read()
entries=int((csv_input.count('>'))/2)##Because that symbol '>' is repeated twice
print ("The total number of miRNA entries in cleaveland file:",entries,'\n')
fh_in.close()


##Open a files to be written in CSV format
fh_output=open('cleaveland_parsed', 'w') # File to compare entries with validated ones
csv_out=csv.writer(fh_output, delimiter=',')

#fh_output2=open('/home/atul/Analysis/unvalidated_scoring_input', 'w')# File to score directly at this step with comparing with validated ones
#csv_out2=csv.writer(fh_output2, delimiter=',')

#Open the input file for operations($$$$$$$$$$$$)
#Please enter the NAME OF CLEAVELAND FILE to be parsed here
fh_in=open(cleave_file, 'r')
csv_in=csv.reader(fh_in, delimiter='\t')
csv_table=[]

##Populate the list
for row in csv_in:
	csv_table.append(row)


#Extract Entry, miRNA, Target
result_entries=0
a_count=0
for i in range(entries):
#	print('An entry:', csv_table[0:5])
	ent_loc=csv_table[0+a_count][0]
	miDNA=csv_table[0+a_count][1]
	cleave_site=csv_table[0+a_count][3]
	tar_loc=csv_table[2+a_count][0]
	mi_loc=csv_table[4+a_count][0]
	
	
	block=csv_table[0+a_count][2]
#	print(block)
	
	block_splt=block.split('|')
#	print(block_splt[-1])

	sub_block=block_splt[-1].split(':')
	chr_val=sub_block[0]
	chromo_up=chr_val.strip(' ')
	chromo=chromo_up.lower()## to change 'Chr' to 'chr' as in cleaveland file
	
#	print(chromo)
	

#	print(ent_loc, miDNA, mi_loc, tar_loc, chromo, cleave_site)
#	break
	csv_out.writerow([miDNA]+[chromo]+[cleave_site]+[ent_loc]+[tar_loc]+[mi_loc])
	a_count+=6
	result_entries+=1
	
print('Number of entries in cleaveland parsed file:', result_entries)
	
# Close output file: parsed_out
fh_output.close()
fh_in.close()


###Module to compare two files

	
	







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
