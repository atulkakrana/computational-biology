##Note the parser need s to be re-configured if the PARE validation output have different coulumn structure/numbers

## This script parses the results from PARE validation pipeline - jixian
##Before running this script you may want to filter the PARE pipeline results using PARE_validation_optimization.py
## Immediate old version - parsing_validated_soya.py


import csv 


#Open the validated targets file using csv reader
fh_in=open('PARE_validated_miRNAs', 'r')
csv_validated=csv.reader(fh_in, delimiter='\t')

# Open the file to write parsed results
fh_out=open('PARE_valid_parsed', 'w')
csv_validated_parsed=csv.writer(fh_out, delimiter=',')

#Create an empty list so as to read the file row by row and append to this list- loop run directly on csv_val now
#csv_val_table=[]

##Populate the list
#for row in csv_val:
#    csv_val_table.append(row)
    
validated_entry_count=0
#for i in csv_val_table:
for i in csv_validated:
    miDNA=i[1]
    chromo_up=i[4]
    chromo=chromo_up.lower()## to change 'Chr' to 'chr' as in cleaveland file
    cleave=i[10]
#    print(miDNA, chromo, cleave)
    csv_validated_parsed.writerow([miDNA]+[chromo]+[cleave])
    validated_entry_count+=1

print('The total number of validated entries were:', validated_entry_count)

fh_in.close()
fh_out.close()


'''
Created on Apr 14, 2012

@author: atul
'''
