#!/usr/bin/python3


import os
import sys
import operator
import itertools
import time
import csv
import sys
import argparse
from operator import itemgetter



#################################################### Settings #################################################################
master ='./results/IL_RESULTS_GLOB_P60_Batch_wo_P30C2_oldWBv1.tsv'
query = './test'
delimiter = '\t' ##Delimiter in master
delimiter2 = '\t' ##Delimter in query
position = 2 ## Column numer in masther that has key
position2 = 9 ## Column number in query that has key

#out_file = './results/test_ranked.csv'

match = 'N'##Use matching 

TS = 'Y'
cols=[20,8]##list that hold info about first column to start sorting and number of samples

global rank_col
# rank_col = int(2)##Column to start adding ranks from

FC = 'Y' ##has fold change columns being added before
fc_cols = [4,8]###List that hold info about FC columns, first column to start sorting and second for number of samples

################################################### Functions ###################################################################
#################################################################################################################################

def Extractor(master,query):
    
    ##Read files once for all
    if match == 'Y':
    
        fh_in_query = open(query, 'r')
        fh_in_master = open(master, 'r')
        query_read = fh_in_query.readlines()
        master_read = fh_in_master.readlines()
        
        out_match_file = ('%s_match_%s' % (query,master))
        fh_match_out = open(out_match_file, 'w')
        
        out_mis_file = ('%s_mis_%s' % (query,master))
        fh_mis_out = open(out_mis_file, 'w')
        
        
        
        query_set = set()## To uniq query set
        query_count = 0
        for i in query_read:
            i = i.strip('\n')
            #print(i.split('\t')[0].upper())
            #key = i.split('\t')[0].upper()
            #print('Column from query to match : %s' % (position2))
            key = i.split('\t')[int(position2)-1].upper()###Key is the value on which matching will be done, upper case, index converted to python format
            #print(key)
            query_set.add(key)
            query_count+=1
        print (query_set)    
        print ('**\\nNumber of query processed: %s || Number of non-redumdant queries: %s**\n\n' % (query_count,len(query_set)))
            
        
        #miss_set =set()
        #miss_set =  query_set ##A copy of query set is made which will be used to remove the keys which were found and finally those not found will be left
        
        match_count = 0###To keep count of matched entries
        mis_count = 0 ###TO keep count of mismatches
        key_count = 1###It shows on screen that how many of total keys being searched
        for key in query_set:
            print('Key which is being searched %s out of total keys %s' % (key_count, int(len(query_set))))
            key_count += 1
            
            ent_count = 0 ## entries in master file
            for entry in master_read:
                ent_count+=1
                #print ('%s entry from matched file being analyzed' % (ent_count))
                entry_strp = entry.strip('\n')
                entry_splt = entry_strp.split(delimiter)##index converted to python format
                #print (entry_splt[1].upper())
                if key==entry_splt[int(position)-1].upper():#Both key and query should be at same case
                    fh_match_out.write(entry)
                    match_count += 1
                    
        ## To extract that do not match    
        for entry in master_read:
            entry_strp = entry.strip('\n')
            entry_splt = entry_strp.split(delimiter)##index converted to python format
            if entry_splt[int(position)-1] not in query_set:
                    fh_mis_out.write(entry)
                    mis_count += 1
                    
                    
        fh_match_out.close()
        fh_mis_out.close()
        fh_in_query.close()
        fh_in_master.close()
        
        print('Total number of keys in query: %s | Total number of entries in master file:  %s' % (len(query_set), ent_count))
        print('Total number of matching entries found: %s | Total number of mis-matching entries found:  %s' % (match_count, mis_count))
        
        res_file = out_file
        
    else:
        res_file = master## If no matching need to be done than master file is 
        pass
    
    return res_file

def GenRanks(res_file,cols):
    fh_in = open(res_file,'r')
           
    
    #head = fh_in.readline().strip('\n').split('\t')##Waste Header
    #print (head)
    res = [line.strip('\n').split('\t') for line in fh_in]##Line converted to list and read
    #print (res[0])
    
    print ('The number of results in file: %s\n' % (len(res)-1))
    
    res_list = list(res) ##List to hold results while preocessing
    t_col = int(cols[0])-1 ##First column to start sorting, number refers to column in excel sheet
    fc_col = int(fc_cols[0])-1##First column to start sorting for FC, updated with every new column added while T-stat ranikng
    rank_col = 3 ##PYTHON FORMAT - Start adding ranks at this column@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ##for i in range (2):        
    
    if TS == 'Y':       

        t_rank_start = rank_col###As refrence to rember from where TStat ranks start
        
        for i in range (cols[1]):##for number of columns
    
            head = res_list.pop(0) ## Remove header before sorting and store here
            #print ('Header:\n',head)
            #print ('\nFirst line of result:\n',res[0])
           
            if head[t_col].split('.')[0] == 't': ### Identifies t-stats as column, sort on high to low t-stats only
                print("\n***Sorting on Column: %s'%s'***\n" % (t_col,head[t_col]))
               
                #print ('\nFirst line of result:\n',res[0])
                res_list.sort(key= lambda x: float(x[t_col]),reverse = True)##Sort from high to low, header is excluded
                #print ('Sorted results:\n',res[:2])
               
                ##Generate ranks
                rank = 1
                for ent in res_list:
                    ent.insert(rank_col,rank)##Insert at column you started
                    #print (ent)
                    #res_list.append(ent)##Adding to final list
                    rank+=1
               
               
                ##Add header to sorted list
                title = 'Rank.'+ head[t_col]##New Rank column header
                head.insert(rank_col,title)##Insert new column to header
                #print (new_head,'\n', head)
                res_list.insert(0,head)### Sorted results with header pasted back
                #print ('Header after processing:\n',res_list[0])
               
                ##Increment columns
                t_col+=2## One column is moved due to addition of ranks and column next to last needs to be considered for ranking
                fc_col +=1##updates new column for FC ranking as a new rank column is added
                rank_col+=1
               
                #If ranking complete on one feature i.e t-stat or FC
                if rank_col-t_rank_start == cols[1]:##ranks columns became equal to columns to rank
                   
                    head = res_list.pop(0)##Remove header before min ranks
                    for ent in res_list:
                        
                        block = ent[t_rank_start:rank_col]
                        # print (block)
                        min_rank = min(block)
                        ent.insert(rank_col,min_rank)
                   
                    head.insert(rank_col,'Min.Rank.TStats')
                    res_list.insert(0,head)### Sorted results with header pasted back
                
                    ##Increment columns
                    fc_col +=1##updates new column for FC ranking as a new rank column is added, Dont know why I have to keep it to 2 instead of 1 but it works this way only
                    rank_col+=1 ##Rank column updated if ranking on FC needs to be performed than start at this col
                    print('Columns passed on for Rank %s and FoldChange %s' % (rank_col,fc_col))##Both can be equal because FC sorting done and rank inserted on same column
    
            else:
                print("\n***Error: Please check column used for T-stat sorting: '%s'***\n" % (head[t_col]))
                #print('\n**Please check the column %s specifed for sorting - Bye!**\n' % (col))
                sys.exit()     
            
    ##Rank on basis of P-val and Fold change    
    if FC == 'Y':
        
        # for i in range(1):  
        fc_rank_start = rank_col ##Column from where FC ranks were added, used to calcultae min rank
        
        for i in range (fc_cols[1]):## For all the columns
            head = res_list.pop(0) ## Remove header before sorting and store here
            #print ('Header:\n',head)
            #print ('\nFirst line of result:\n',res[0])
            
            if head[fc_col].split('.')[0] == 'FC': ### Identifies t-stats as column, sort on high to low t-stats only                
                pval_title = ('p.value.%s' % head[fc_col].split('.')[-1])## Matching p_value column
                #print('pval_title:', pval_title)
                ##for i in head:
                ##   print(i)
                ##   if i.strip() == pval_title:
                ##       print ('Yes')

                                    
                #print ("Here is header:\n",head)
                pval_col = head.index(pval_title)##Get column number to sort
                # print('\nThis is the p value title %s and Column number %s\n' % (pval_title, pval_col))

                print("\n***Sorting on Columns: '%s' and '%s'***" % (head[pval_col],head[fc_col]))
                res_list.sort(key= lambda x: (float(x[int(pval_col)]),-float(x[fc_col])))##Sort from high to low, header is excluded
                # print ('Sorted results:',res[-5:-1],'\n')
                
                ##Generate ranks
                rank = 1
                print ('***Generating P-val+FC ranks for time point: %s***\n' % (head[fc_col].split('.')[-1]))
                for ent in res_list:
                    ent.insert(rank_col,rank)##Insert at column you started or left by T-stat ranking
                    # print (ent)
                    # res_list.append(ent)##Adding to final list
                    rank+=1
                    
                 ##Add header to sorted list
                title = 'Rank.'+ head[fc_col]##New Rank column header
                head.insert(rank_col,title)##Insert new column to header
                #print (new_head,'\n', head)
                res_list.insert(0,head)### Sorted results with header pasted back
                # print ('Header after processing:\n',res_list[0])
                
                ##Increment columns
                # print ('Incrementing')
                fc_col +=2##updates new column for FC ranking as a new rank column is added
                rank_col+=1

                # ##If ranking complete on one feature i.e t-stat or FC
                if rank_col-fc_rank_start == fc_cols[1]: ## Calculate min rank when columns equal to numer of samples have been added
                    print ('Calculating min rank')
                    head = res_list.pop(0)##Remove header before min ranks
                    for ent in res_list:
                        block = ent[fc_rank_start:rank_col]
                        min_rank = min(block)
                        ent.insert(rank_col,min_rank)

                    head.insert(rank_col,'Min.Rank.FC')
                    res_list.insert(0,head)### Sorted results with header pasted back
            
            else:
                print("\n***Error: Please check column used for P-val and FC based sorting: '%s'***\n" % (head[fc_col]))
                sys.exit() 
    
    #print (res_list[1])
            
    return res_list

def WriteRes(res_list,res_file):
    
    out_file = ('.%s_Rank.tsv' % (res_file.split('.')[-2]))
    print ('\nThe result file is: .%s_Rank.tsv\n' % (res_file.split('.')[-2]))
    res_out = open (out_file, 'w')
    writer = csv.writer(res_out, delimiter='\t', quoting=csv.QUOTE_ALL)
    for i in res_list:
        #print(i)
        writer.writerow(i)
    
    res_out.close()
    
    return res_out

def main():  
  
    res_file = Extractor(master,query)
    
    res_list = GenRanks(res_file,cols)
    WriteRes(res_list,res_file)
    

if __name__ == '__main__':
    
    start = time.time()###time start
    main()
    end = time.time()
    print ('The script run time is %s' % (round(end-start,2)))
    print('The run has completed sucessfully.....CHEERS! - Exiting........\n')
    sys.exit()

