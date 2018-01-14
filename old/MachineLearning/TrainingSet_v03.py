#!/usr/local/bin/python3
### Input valid and TF pased file and set match to 'Y' and Uniq to 'N' - Run script - From pos_pattern_out_ext filter the results and name as pos_pattern_out_ext_filter.csv
### Now Uniq the entries from this new file and existing neg_pattern_out.csv by setting match to 'N' and Uniq to 'Y' and rerun script

####MODIFIED TO UNIQ VALIDATION RESULTS ON BASIS OF MIRSEQ BINDING SITE AND CLEVAGE SITE

import re
import sys
import os
import subprocess
import csv 
import glob

######### INPUT #########
Match = 'N' ### Will give you positive and negative files
TF_file = 'TF_res_parsed.csv'
valid_file = 'sco_inp_ext_geno' ## Sco_inp_ext_geno has reverse mapped coordinates, ass abundance info from sco_inp_geno file and than run this

Uniq = 'Y'#### Input positive file (filtered) and negative file - to uniq interaction in both on basis of miRNA and Target seq - This is followed by scoring script
fls = ['pos_pattern_out_ext_filter.csv','neg_pattern_out.csv']
sep = '\t'

SizeSeparate == 'N'## Break results in different parts

###Get positive and negative files with target information by matching valid results gainst target finder all predictions
def MatchCLandTF(valid_file,TF_file):

    ##This module matches CL files (one from each PARE library) and TF file (single file) to generate files (one for each PARE library) with sequence information for scoring  
    ##Because there are multiple files from CL output, one for each PARE library from the database
    ##Prepare output files and name according to library 
    
    positive_file = ('./pos_pattern_out.csv')## This file should be combined at the end of module to unify results from all library
    fh_out_pos=open(positive_file, 'w') 

    negative_file = ('./neg_pattern_out.csv')
    fh_out_neg=open(negative_file, 'w')
        
    ### Matching starts#################################################################################
    validated_tuples = []###List to hole validated entries
    info_tuples = [] ###List to hold additional information for positive set
    ##Make Dictionary of mir seq and tar seq from valid file - To compare with TF file and dive it into Positive and negative
    #print ('\nThe file %s is being matched with TargetFinder Results')
    with open(valid_file) as fh1:
        #fh1.readline()### Kill the header
        csv_reader1 = csv.reader(fh1)
        for row in csv_reader1:
            #print(row)
            #print(row[7],row[8])
            validated_tuples.append(tuple((row[7],row[8]))) #### miRNA seq, Target seq
            #info_tuples.append(tuple((row[0],row[1],row[7],row[8]))) ### Used to map important information later
    print('The total number of entries i.e tuples in %s are %s' % (valid_file, len(validated_tuples)))

    ##Read TF file and compare with the column 1-3 in 'Validates_tuples' list 
    pos_list = []####The function of list to hold matching target finder entries which will be used later to join TF and CL entries
    neg_list = []### Not validated 
    with open(TF_file) as fh2:
        fh2.readline()
        csv_reader2 = csv.reader(fh2)
        match_count=0
        for row in csv_reader2:
            #print (row)
            #print(row[3],row[4])
            if tuple((row[3],row[4])) in validated_tuples: ### Matched with miRNA, Target seq
                pos_list.append(row[0:5])
                fh_out_pos.write('%s,%s,%s,%s,%s,%s\n' % (row[0],row[1],row[2],row[3],row[4],row[5]))
                match_count+=1
                
            else:
                #print(row[1:5])
                fh_out_neg.write('%s,%s,%s,%s,%s,%s\n' % (row[0],row[1],row[2],row[3],row[4],row[5]))
                neg_list.append(row[0:5])

    print('Total number of entries joined with target finder final results:%s\n' % (match_count))
#        print('Please use the output files %s and supplemental file %s  for further analysis' % (matched_file,matched_supp ))

    fh_out_pos.close()
    fh_out_neg.close()
    
    
    
    ####Add p-value, abundances from geno file back to positive file#####
    positive_file_ext = ('./pos_pattern_out_ext.csv')## This file version will have p-value and abundances at end
    fh_out_pos=open(positive_file_ext, 'w')
    pos_list_ext = []### Includes extended information of postive file
    
    with open(positive_file) as fh3: ###Open original poistive file and make key
        csv_reader3 = csv.reader(fh3)
        for row in csv_reader3:
            info_tuples.append(tuple((row[0][1:],row[1],row[3],row[4]))) ### Key of mirname (without >), genename,mir seq, geneseq

    with open(valid_file) as fh4: ### Use key to match with valid file/extend geno
        csv_reader4 = csv.reader(fh4)
        for row in csv_reader4:
            print ((row[0],row[1],row[7],row[8]))
            if tuple((row[0],row[1],row[7],row[8])) in info_tuples: ### Matched with miRname, Gene name, miRNA, Target seq
                pos_list_ext.append((row[0],row[1],row[4],row[7],row[8],row[9],row[10],row[11],row[12]))
                fh_out_pos.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (row[0],row[1],row[4],row[7],row[8],row[9],row[10],row[11],row[12]))
                match_count+=1
            else:
                print('Not found')

    return positive_file,negative_file,positive_file_ext,pos_list,neg_list,pos_list_ext


## Get the uniq miRNA-target pattern
def UniqPattern(fls):
    
    #master_list = [pos_list_ext,neg_list]
    #namelist = ['pos','neg']
    #name_num = 0 
    #for alist in master_list:
    for afile in fls:   
        print ('\n**File being analyzed: %s**' % (afile))
        fh_in = open(afile,'r')
        
        #outfile = '%s_uniq_pattern' % (namelist[name_num])
        outfile = '%s_uniq_pattern' % (afile.split('.')[0])
        fh_out=open(outfile, 'w') 
        
        ###
        added_keys=set()## A set to store first 3 elements from input file: miRNA seq, chr# and cleavage site and than use it to compare further entries in file
        parsed_out_count=0## TO keep count of unique entries
        #for ent in alist:
        for line in fh_in:
            ent =line.strip('\n').split(sep)            
            #print(ent[3],ent[4])
            #genename = ent[1].split('.')[0]##To avoid different variations of same gene to be counted as uniq
            #print (genename)
            lookup=tuple((ent[2],ent[3],ent[5]))##binding_dite, cleavage site and miR_seq
            if lookup not in added_keys:##That means such a tuple has not been recorded yet and is unique
                fh_out.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[1],ent[2],ent[3],ent[4],ent[5],ent[6],ent[7],ent[8],ent[9],ent[10],ent[11]))###format same as traget finder 
                #csv_out.writerow(ent[0:5])##MiRNA name, miRNA seq and target seq
                parsed_out_count+=1
        #        outlist.append(row)
                added_keys.add(lookup)## Once a new entry is found it is recorded so as to compare and neglect further entries
            else:
                pass
                
        print('The number of unique entries found and will be used for scoring: %s\n'% (parsed_out_count))
        
        #print(outlist)
        #print('The total number of entries after removing duplicates:', len(outlist))
        fh_in.close()
        fh_out.close()
        
      
        ###############  Module 3- Segregate miRNAs on the basis of size before actual scoring - eliminates the need to run count_mir_len
        if SizeSeparate == 'Y':
            ##The output files corresponding to the length of miRNA
            sco_21 = '%s_sco_pat_21' % (afile.split('.')[0])
            fh1_out=open(sco_21, 'w')
            sco_22 = '%s_sco_pat_22' % (afile.split('.')[0])
            fh2_out=open(sco_22, 'w')
            sco_23 = '%s_sco_pat_23' % (afile.split('.')[0])
            fh3_out=open(sco_23, 'w')
            sco_24 = '%s_sco_pat_24' % (afile.split('.')[0])
            fh4_out=open(sco_24, 'w')
    
            ##The input file to be read and filtered according to their length
            fh_in1=open(outfile, 'r')
                  
            ## Read the file line by line and do operation
            for ent in fh_in1:
                #print (ent.strip('\n'))
                i = ent.strip('\n').split(',')
                #print(i,'\n',i[3])
                a=i[3].count('-')##Count the bulges and gap in miRNA   
            #    print('*No. of Gaps and bulges:', a)
            #    print('miR:', i[2], ' | miR nt:', len(i[2]))
                original_len=int(len(i[3])-a)###Length of miRNA
            #    print('Correct length:', original_len)
                if original_len == 21:
                    fh1_out.write('%s,%s,%s\n' % (i[0],i[3],i[4]))
                    #csv1.writerow(i[1]+i[4:5])
                elif original_len == 22:
                    fh2_out.write('%s,%s,%s\n' % (i[0],i[3],i[4]))
                    #csv2.writerow(i[1]+i[4:5])
                elif original_len == 23:
                    fh3_out.write('%s,%s,%s\n' % (i[0],i[3],i[4]))
            #        print(i[2])
                    #csv3.writerow(i[1]+i[4:5])
                elif original_len == 24:
                    fh4_out.write('%s,%s,%s\n' % (i[0],i[3],i[4]))
                    #csv4.writerow(i[1]+i[4:5])
                else:
                    print('An miRNA of length :', original_len, 'was found')
                    print(i[1],i[3])
                    pass            
            print('\nPlease use the output file of required miRNA length for next step of scoring')
            
            
            ###Make a combines file of 21 and 22 mers
            #os.cat()
            fh_in1.close()
            fh1_out.close()
            fh2_out.close()
            fh3_out.close()
            fh4_out.close()
            
            
        #return outfile 

def main():
    
    if Match == 'Y':
        positive_file,negative_file,positive_file_ext,pos_list,neg_list,pos_list_ext = MatchCLandTF(valid_file,TF_file)
    if Uniq == 'Y':
        UniqPattern(fls)
    

if __name__ == '__main__':
    main()
    print ('Cheers ..!_!.. Script finished succesfully')
    sys.exit()
    
    

###V01 -> v02
###Sco_inp_ext_geno file have co-ordinates reveerse mapped and therefore cannot be used to compare with target finder
## Using Sco_inp_ext_geno as input is important because it has p-value and window information and can be directly filtered
##unlike sco_input_extend. It is therefore to use Sco_inp_ext_geno as input maching keys were changed to just miR and Tar seq.
###The idea is if one interaction is validated than it does not matter that waht was the binding site ot was PARE mapping was from different site
####pattern identification only cares about a  validated interaction thats all.

####Broken into two steps
##Match - Match provides two files those interaction which were validated and those which were predicted but not validated
##uniq - With positive and negative files - One might further need to filter both to increase confidence - Therefore Uniq is a separate process - Feed filtered files from 'match' in input parameter
    
    




