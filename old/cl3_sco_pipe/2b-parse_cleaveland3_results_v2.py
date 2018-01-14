##Note the parser needs to be re-configured if the PARE validation output have different coulumn structure/numbers

## This script parses the results Cleaveland 3 pipeline for further analysis



import csv 

##Config: Enter the filename you got from final cleaveland3 analysis
cl3_result='Gmax_cleaveland3_final_res' 

#Open the Cleaveland 3 final results file using csv reader
fh_in=open(cl3_result, 'r')
fh_in.readline()##Header line wasted
csv_results=csv.reader(fh_in, delimiter='\t')

# Open the file to write parsed results
fh_out=open('cl3_results_parsed.csv', 'w')
#Write Header
fh_out.write('miRNA name,target_gene,target_location,targetfinder_score,p-value\n')
#cl3_results_parsed=csv.writer(fh_out, delimiter=',')

 
entry_count=0

for i in csv_results:
#    print(i)    
    mirna=i[0]
    target=i[1]
    range=i[3]
    score=i[2]
    p_val=(i[6])
    fh_out.write('>%s,%s,%s,%s,%s\n' % (mirna,target,range,score,p_val))
#    cl3_results_parsed.writerow([mirna]+[target]+[range]+[score]+[p_val])
    entry_count+=1
#    print(mirna,target,range,score,p_val)
    

print('The total number of result entries were:', entry_count)
print('The output file "cl3_results_parsed.csv" is generated for further analysis')

fh_in.close()
fh_out.close()


'''
Created on Apr 14, 2012

@author: atul
'''
