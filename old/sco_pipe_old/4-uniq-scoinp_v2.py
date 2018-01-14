import csv

##Modules
## 1. Remove redundant entries-with all the attributes same i.e DNA, chr#, cleavage site, entry name, tar, miRNA
## 2. Write different scoring inputs for final scoring steps i.e sco_inp_21, sco_inp_22 etc 

##Reading parsed output having redundant entries
fh_in=open('validated_out', 'r') 
parsed_in=csv.reader(fh_in, delimiter=',')
#for i in parsed_in:
#    print(i)

## Output file
fh_output=open('scoring_input', 'w') 
csv_out=csv.writer(fh_output, delimiter=',')

##Method 1 and 2 - check old version named 'remove_redundant.py'
##Method:3- To find unique entries from the 'Matched i.e validated entries from cleaveland output'

#outlist=[]## To get all the unique values as a list, replaced by taking output to a file using csv.writer in loop below

added_keys=set()## A set to store first 3 elements from input file: miRNA-DNA, chr# and cleavage site and than use it to compare further entries in file

parsed_out_count=0## TO keep count of unique entries
for row in parsed_in:
    lookup=tuple(row[:3])
    if lookup not in added_keys:##That means such a tuple has not been recorded yet and is unique
        csv_out.writerow(row[3:6])
        parsed_out_count+=1
#        outlist.append(row)
        added_keys.add(lookup)## Once a new entry is found it is recorded so as to compare and neglect further entries
    else:
        pass
        
print('The number of unique entries found and will be used for scoring:', parsed_out_count )

#print(outlist)
#print('The total number of entries after removing duplicates:', len(outlist))

fh_in.close()
fh_output.close()

###############  Module 2- Segregate miRNAs on the basis of size before actual scoring - eliminates the need to run count_mir_len

##The output files corresponding to the length of miRNA
fh1_out=open('scoring_inp_21', 'w')
csv1=csv.writer(fh1_out, delimiter=',')
fh2_out=open('scoring_inp_22', 'w')
csv2=csv.writer(fh2_out, delimiter=',')
fh3_out=open('scoring_inp_23', 'w')
csv3=csv.writer(fh3_out, delimiter=',')
fh4_out=open('scoring_inp_24', 'w')
csv4=csv.writer(fh4_out, delimiter=',')


##The input file to be read and filtered according to their length
fh_in1=open('scoring_input', 'r')
csv_in1=csv.reader(fh_in1, delimiter=',')

## Read the file line by line and do operation
for i in csv_in1:
    a=i[2].count('-')##Count the bulges and gap    
#    print('*No. of Gaps and bulges:', a)
#    print('miR:', i[2], ' | miR nt:', len(i[2]))
    original_len=int(len(i[2])-a)
#    print('Correct length:', original_len)
    if original_len == 21:
        csv1.writerow(i[0:3])
    elif original_len == 22:
        csv2.writerow(i[0:3])
    elif original_len == 23:
#        print(i[2])
        csv3.writerow(i[0:3])
    elif original_len == 24:
        csv4.writerow(i[0:3])
    else:
        print('An miRNA of length :', original_len, 'was found')
        print(i[2])
        pass
    
print('\nPlease use the output file of required miRNA length for next step of scoring')



    




'''
Created on Mar 1, 2012

@author: atul
'''


'''
Created on Apr 14, 2012

@author: atul
'''
