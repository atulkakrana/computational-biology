#!/usr/local/bin/python3

###This module is to extract entries from one file using a list of values in other file
##Please do not use commandline valued for delimiter instead mention in script
##As output you will get two files one with entries matched query and other entries that do not match to query
##Please keep in mind that this did not replace ;Extractheaders' script which is used to extract FASTA sequences or entries matching a header

 
import argparse

##- - - - - - - -- - - - - - Arguments - - - - - - - - - - - - - -##
parser = argparse.ArgumentParser()

##Command line options
parser.add_argument("-m", "--master", dest = "master", help="master file")
parser.add_argument("-q", "--query", dest = "query", help="queries to be extracted")
parser.add_argument("-d", "--delimiter", dest = "delimiter",default = '\t', help="delimiter in master")
parser.add_argument("-p", "--position", dest = "position",default = '1', help="position/column of value in master")
parser.add_argument("-d2", "--delimiter2", dest = "delimiter2", default = '\t', help="delimiter in query")
parser.add_argument("-p2", "--position2", dest = "position2", default = '1', help="position/column of value in query")

args = parser.parse_args()

print( "master {} query {} delimiter {} delimiter2 {} position {} position2 {} ".format(
        args.master,
        args.query,
        args.delimiter,
        args.delimiter2,
        args.position,
        args.position2,     
        ))

def Extractor(master,query):
    
    ##Read files once for all
    
    fh_in_query = open(query, 'r')
    fh_in_master = open(master, 'r')
    query_read = fh_in_query.readlines()
    master_read = fh_in_master.readlines()
    
    out_match_file = ('%s_match_%s' % (query,master))
    fh_match_out = open(out_match_file, 'w')
    
    out_mis_file = ('%s_mis_%s' % (query,master))
    fh_mis_out = open(out_mis_file, 'w')
    
    
    
    query_set = set()## To uniq query set
    for i in query_read:
        i = i.strip('\n')
        #print(i)
        key = i.split('args.delimiter2')[int(args.position2)]###Key is the value on which matching will be done
        #print(key)
        query_set.add(key)
        
    
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
            entry_splt = entry_strp.split('args.delimiter')
            #print (entry_splt)
            if key==entry_splt[int(args.position)]:
                fh_match_out.write(entry)
                match_count += 1
                
    ## To extract that do not match    
    for entry in master_read:
        entry_strp = entry.strip('\n')
        entry_splt = entry_strp.split(',')
        if entry_splt[int(args.position)] not in query_set:
                fh_mis_out.write(entry)
                mis_count += 1
                
                
    fh_match_out.close()
    fh_mis_out.close()
    fh_in_query.close()
    fh_in_master.close()
    
    print('Total number of keys in query: %s | Total number of entries in master file:  %s' % (len(query_set), ent_count))
    print('Total number of matching entries found: %s | Total number of matching entries found:  %s' % (match_count, mis_count))
    
    return query_set, out_match_file
    
    
    

##Main


query_set,out_file = Extractor(args.master, args.query)             
                
        
    
    
        



'''
Created on Sep 8, 2012

@author: atul
'''
