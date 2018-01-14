# - - - - - H E A D E R - - - - - - - - - - - - - - - - -



# - - - - - U S E R    V A R I A B L E S - - - - - - - -
##Config
cleave_file = 'comb_res_target_finder'


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -



# - - - - - - F U N C T I O N S- - - - - - - - - - - - - - - 

import csv

# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -

##Open a files to be written in CSV format
fh_output=open('targetfinder_results_parsed.csv', 'w') # File to compare entries with validated ones
#csv_out=csv.writer(fh_output, delimiter=',')

##Write header to outputfile
fh_output.write('miRNA name,target_gene,target_location,mirna_seq,tar_seq,targetfinder_score\n')

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

	fh_output.write('>%s,%s,%s,%s,%s,%s\n' % (mirna,target,range,mir_seq,tar_seq,score))
#	csv_out.writerow([mirna]+[target]+[range]+[mir_seq]+[tar_seq]+[score])
	a_count+=6
	b_count+=6
	result_entries+=1
	
print('Number of entries in cleaveland parsed file:', result_entries)
print('The output file with parsed results generated is "targetfinder_results_parsed.csv"' )
	
# Close output file: parsed_out
fh_output.close()
fh_in.close()


###Module to compare two files

	
	







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - S U B R O U T I N E S - - - - - - - - - - - - - -



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
