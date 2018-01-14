##Make sure TF comined file had an empty line at end


import mysql.connector as sql
import re
import sys
import os
import subprocess
import csv 


##This script is for NBLRR analysis


def TargetFinderParse(inp_file):

    ##Open a files to be written in CSV format
    TF_parsed_file = 'TF_res_parsed.csv'
    fh_output=open(TF_parsed_file, 'w') # File to compare entries with validated ones   
    ##Write header to outputfile
    fh_output.write('miRNA name,target_gene,target_location,mirna_seq,tar_seq,targetfinder_score\n')
        
    ##Open concatenated target finder results file for reading    
    fh_in=open(inp_file, 'r')##Can give errors due to file opening and closing
    csv_in=csv.reader(fh_in, delimiter='\t')
    csv_table=[]

    ##Populate the list
    for row in csv_in:
        csv_table.append(row)
    #print (csv_table)
    
    #print(csv_table[0:12])        
    #print(csv_table[0:12])##+6 to get next block i.e 6 lines makes an entry
    
    entries=int(len(csv_table)/6)
    print('Total number of entries in input file:',entries,'\n')


    #Extract Entry, miRNA, Target
    scr_re = re.compile('score=\d{1,2}')
    ran_re = re.compile('range=\d*-\d*') 
    result_entries=0
    a_count=0###incremented to move the next block
    b_count=6
    for i in range(entries):
        score =''
        ent_loc=csv_table[a_count:b_count]##A miR and target entry complete
        info_block=ent_loc[0][0].split('|')##block containing names and other info
#        print('%s\n' % (info_block))
        
        mirna=info_block[0].split(',')[0].split('=')[1]
        target=info_block[0].split(',')[1].split('=')[1].strip(' ')
#        print(info_block[3])
        
        scr = re.search(scr_re, str(ent_loc[0]))## make string of list entry before regex search
        score = scr.group().split('=')[1]##scr.group()returns the matched value
        
        ran = re.search(ran_re, str(ent_loc[0]))## make string of list entry before regex search
        arange = ran.group().split('=')[1]  
        
        tar_seq=ent_loc[2][0].split()[2]
        mir_seq=ent_loc[4][0].split()[2]
        
        print(mirna, target, score, arange, mir_seq, tar_seq)


        fh_output.write('>%s,%s,%s,%s,%s,%s\n' % (mirna,target,score,arange,mir_seq,tar_seq))
        a_count+=6
        b_count+=6
        result_entries+=1
    #    
    print('Number of entries in %s file: %s' % (inp_file,result_entries))
    print('The output file with parsed results generated is "targetfinder_results_parsed.csv"' )
        
    # Close output file: parsed_out
    fh_output.close()
    fh_in.close()
    
    return TF_parsed_file

##Main

inp_file = '/data2/homes/kakrana/bugs/1279/1.TF_ANALYSIS/at_target_finder_results/test'
TF_parse_file = TargetFinderParse(inp_file)



'''
Created on Jul 31, 2012

@author: atul
'''
