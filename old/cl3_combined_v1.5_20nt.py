##This script makes use of Cleaveland3 pipeline to predict PARE validated targets. This version of script is written for biologists as some significant changes for them are made like below:
##Originally final results were filtered on basis of miR seq+targetname+cleavage site but for biologist miRseq was replaced by miR name
##You may want to run CL2 2nd script which creates an excel sheet of PARE reads in small and large windows because there could be non-specificity in targets.
##MiRNAs from same family can have same sequence and thus targets so you have to filter by using CL2 Jixian 2nd script
##This script was written by Atul kakrana: atulkakrana@gmail.com

##


import mysql.connector as sql
import re
import sys
import os
import subprocess
import csv 
import glob


##Config - This part needs to be improved for user friendliness
#1. Input cDNA/transciptome file name
trascriptome_file = 'mt3.5_transcripts.fa_new_head.fa'##The FASTA header will be cleaned automatically

##2. miRNA File
miRNA_file = 'mt_pnl_mir.fa'##The header will be cleaned automatically

##3. Server name in 'connectToDB' module
##Please change the name of server where your DB exists

##4. PARE database names
##Please change the PARE db name at two queries in 'TagAbundanceFile' module

##5. Cleaveland folder
##Change CleaveLand address if required in 'CleavelandAnalysis' module

##6. Index filename
##Please change index filename in 'mapdd2trans' module. Index should be made on cDNA file to be used in this study and must have clean header, use 'clean_fasta_header_v1' for that.

##7. Bowtie folder
##Change folder address if required in 'mapdd2trans' module

##8. You may want to change length of PARE tags in 'CleavelandAnalysis'


#########AA###################################TT###############################################UU############################################################LL##################


##Module to clean headers
def CleanHeader(filename):
    #read file
    fh_in=open(filename, 'r')
    #write file
    out_file = ('%s_new_head.fa' % (filename))
    fh_out =open(out_file, 'w')
    
    print ('Processing "%s" file to clean FASTA headers\n' % (filename))
    
    acount = 0 ## count the number of entries
    for i in fh_in:
        if re.match('>', i):
            header = i.split()##change here if separater is not white space
            new_head = header[0].split('|')[0]## Change here for number of fields you want strating from 0
            fh_out.write('%s\n' % new_head)
            acount+=1
    #        print(i)
    #        print(new_head)
        else:
            fh_out.write('%s' % i)
        
    fh_in.close()
    fh_out.close()
    return out_file

    print('The fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file, acount))

##Module to connect to DB--------INPUT SERVER NAME HERE-------------------
def ConnectToDB():
    print ('\nTrying to connect to mySQL server')
    # Try to connect to the database
    try:
        con=sql.connect(host='raichu.dbi.udel.edu', user='kakrana', passwd='livetheday')
        print ('\nConnected to Database\n')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con

####Module to get PARE data-----INPUT PARE DB NAME HERE----DONE TWICE----------------
def TagAbundanceFile(con):
    ##Get the number of libraries and make a list of their lib_ids:
    # Set cur to be the cursor so we can execute a query
    os.mkdir('./PARE')
    cur = con.cursor()
    cur.execute('select distinct(lib_id) from LEGUME_priv_PARE_v4.run_master')##@@@@@@@@@@@@@@@@@CHANGE PARE DB NAME HERE@@@@@@@@@@@@@
    libs = cur.fetchall()
#    print (libs)
    print ('\nTotal number of PARE libraries found: %s\n' % (len(libs)))
    for lib in libs:##For all the libraries
#        print (lib[0])
        print ('Caching tag and count information from server for lib %s' % (lib[0]) )
        cur.execute('select tag, raw_value from LEGUME_priv_PARE_v4.run_master where lib_id = %s', (lib[0],))##@@@@@@@@@@@@@@@@@CHANGE PARE DB NAME HERE@@@@@@@@@@@@@
        lib_info = cur.fetchall()
#        print(lib_info[0])
        
        ##Write to file
        print('Writing PARE Fasta file for lib: %s\n' % (lib [0]))
        fh_out = open('./PARE/%s_PARE_tags.fa' % (lib[0]), 'w')##Naming file with lib_ids name
        tag_num = 1
        for ent in lib_info:##All the entries of the library
            if len(ent[0]) == 20:              
                for count in range(ent[1]):##Number of times the tag_count file
                    fh_out.write('>%s\n%s\n' % (tag_num, ent[0]))
                    tag_num += 1
            else:
                #print ('Length is not 20nt')
                pass
            
    return libs

####Module to map degradome to transcriptome/geome----INPUT INDEX NAME HERE---------------
def mapdd2trans(libs):##Cleaveland 3 pipeline 
    os.mkdir('./dd_map')
    mismatch = str(1)
    report_seq = str(2)
    for lib in libs:
        print ('\nThe library %s is being mapped to transcriptome index file' % (lib[0]))
        dd_file = ('./PARE/%s_PARE_tags.fa' % (lib[0]))
        map_out = ('./dd_map/%s_map' % (lib[0]))
#        fh_in = openopen('%s_PARE_tags.fa' % (lib[0]), 'r')
        pipe = subprocess.Popen(["bowtie", "-v", mismatch,"--best", "--strata", "-k", report_seq, "-f", "/data2/homes/kakrana/tools/bowtie-0.12.7/indexes/mt3.5_transcripts_index", dd_file, map_out])##@@@@@@CHANGE INDEX FILENAME HERE@@@@BOWTIE FOLDER CHANGE

    return libs##Just like that

def CleavelandAnalysis(libs,filename):##Both summarize density and Run final cleaveland script | filename is clean Header cDNA_file
#    os.cd('./dd_map')
    print ('\nStarting Cleaveland3 based analysis\n')
    os.mkdir('./dd_density')
    os.mkdir('./cl_results')
    seq_len = str(20)##For mapping PARE data on transcriptome
    pval = str(0.025)##For final cleveland command
   
    for lib in libs:
        map_file = ('./dd_map/%s_map' % (lib[0]))
        density_file = ('./dd_density/%s_density' % (lib[0]))
        fh_out = open(density_file, 'w')
        run_info_file = ('./dd_density/%s_density_run_info' % (lib[0]))
        print ('Summarizing degradome map for PARE lib:%s' % (lib[0]))
##Summarize Degradome
        pipe = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/CleaveLand3_map2dd.pl", "-d", map_file, "-t", filename, "-s", seq_len, "-f", "bowtie"], stdout=subprocess.PIPE, universal_newlines=True)##@@@@@@@@CHANGE CLEAVELAND FOLDER HERE
        results = pipe.stdout.read()
        fh_out.write('%s' % (results))


        ## Cleaveland analysis final run
        cl_results = ('./cl_results/%s_cl_results' % (lib[0]))
        fh_out2=open(cl_results, 'w')
        cl_results_info = ('./cl_results/%s_cl_results_info' % (lib[0]))
        print ('Running final cleaveland analysis step for PARE lib:%s\n' % (lib[0]))
        pipe2 = subprocess.Popen(["/data2/homes/kakrana/tools/3-Cleaveland/CleaveLand3_analysis.pl", "-d", density_file, "-t", "./target_finder_results/", "-p", pval], stdout=subprocess.PIPE, universal_newlines=True)##@@@@@@@@CHANGE CLEAVELAND FOLDER HERE
        results2 = pipe2.stdout.read()
        fh_out2.write('%s' % (results2))

def CleavelandParse(libs):
       
    for lib in libs:
        cl_results = ('./cl_results/%s_cl_results' % (lib[0]))##filename you got from final cleaveland3 analysis
        #cl_results = ('./cl_results/%s_cl_results' % (lib))##Test line
        print('\nCleaveland result file for lib %s is being parsed' % (cl_results))
        #Open the Cleaveland 3 final results file using csv reader
        fh_in=open(cl_results, 'r')
        fh_in.readline()##Header line wasted
        csv_results=csv.reader(fh_in, delimiter='\t')   
    
        # Open the file to write parsed results
        cl_parsed = ('./cl_results/%s_cl_results_parsed.csv' % (lib[0]))
        #cl_parsed = ('./cl_results/%s_cl_results_parsed_test.csv' % (lib))##Test line
        fh_out=open(cl_parsed, 'w')
        #Write Header
        fh_out.write('miRNA name,target_gene,target_location,Cleavage_site,targetfinder_score\n')
    
     
        entry_count=0        
        for i in csv_results:
        #    print(i)    
            mirna=i[0]
            target=i[1]
            arange=i[3]
            cleavage = i[4]
            score=i[2]
            p_val=(i[6])
            fh_out.write('>%s,%s,%s,%s,%s\n' % (mirna,target,arange,cleavage,score))
            entry_count+=1           
        
        print('The total number of result entries were:', entry_count)
        print('The output file "%s" is generated for further analysis\n' % (cl_parsed))
    
        fh_in.close()
        fh_out.close() 
        
    return(libs)##Just like that, of no use


def TargetFinderParse(TFfoldername):
    ##This module combines all the result files from target finder and than parses the single result file to be used for scoring  
    os.mkdir('./scoring')
    ##Cat all files in target finder folder
    
    files_to_cat = TFfoldername+'/*' 
    fls = glob.glob(r'./target_finder_results/*')###NEED to include files to cat in this line instead of explicitly mentioning folder
    TF_res_comb = ('./scoring/TF_results_all')
    TF_out = open(TF_res_comb, 'w') 
    for x in fls:
        TFfile = open(x, 'r')
        data = TFfile.read()
        TFfile.close()
        TF_out.write(data)
    TF_out.close()
    
#    TF_out.write('%s' % (results))
#    TF_out.close()    
    
    ##Prepare output files
    
    ##Open a files to be written in CSV format
    TF_parsed_file = './scoring/TF_res_parsed.csv'
    fh_output=open(TF_parsed_file, 'w') # File to compare entries with validated ones   
    ##Write header to outputfile
    fh_output.write('miRNA name,target_gene,target_location,mirna_seq,tar_seq,targetfinder_score\n')
    
    ##Open concatanated file for reading    
    fh_in=open(TF_res_comb, 'r')##Can give errors due to file opening and closing
    csv_in=csv.reader(fh_in, delimiter='\t')
    csv_table=[]

    ##Populate the list
    for row in csv_in:
        csv_table.append(row)
    
    #print(csv_table[0:12])        
    #print(csv_table[0:12])##+6 to get next block i.e 6 lines makes an entry
    
    entries=int(len(csv_table)/6)
    print('\nTotal number of entries in combined TargetFinder results file:',entries)


    #Extract Entry, miRNA, Target
    result_entries=0
    a_count=0###incremented to move the next block
    b_count=6
    for i in range(entries):
        ent_loc=csv_table[a_count:b_count]##A miR and target entry complete
    #    print(ent_loc)
        info_block=ent_loc[0][0].split(',')##block containing names and other info
    #    print(info_block)
        mirna=info_block[0].split('=')[1]
        target=info_block[1].split('=')[1]
        score=info_block[2].split('=')[1]
        loc=info_block[3].split('=')[1]
        
        tar_seq=ent_loc[2][0].split()[2]
        mir_seq=ent_loc[4][0].split()[2]
        
    #    print(mirna,target,score,range, mir_seq, tar_seq)
    ##    break

        fh_output.write('>%s,%s,%s,%s,%s,%s\n' % (mirna,target,loc,mir_seq,tar_seq,score))
        a_count+=6
        b_count+=6
        result_entries+=1
        
    print('Number of entries in %s file: %s\n' % (TF_parsed_file, result_entries))
        
    # Close output file: parsed_out
    fh_output.close()
    fh_in.close()
    
    return TF_parsed_file

def MatchCLandTF(libs,TF_parsed_file):

    ##This module matches CL files (one from each PARE library) and TF file (single file) to generate files (one for each PARE library) with sequence information for scoring  
    
    ##Because there are multiple files from CL output, one for each PARE library from the database
    for lib in libs:
        
        ##Prepare output files and name according to library 
        matched_file = ('./scoring/%s_validated_out' % (lib[0]))## This file should be combined at the end of module to unify results from all library
        #matched_file = ('./scoring/%s_validated_out_test' % (lib))## test line
        fh_output=open(matched_file, 'w') 
        #csv_out=csv.writer(fh_output, delimiter=',')
        
        matched_supp = ('./scoring/%s_validated_out_supp' % (lib[0]))
        #matched_supp = ('./scoring/%s_validated_out_supp_test' % (lib))##test line
        fh_out_supplemental=open(matched_supp, 'w')
        csv_out_sup=csv.writer(fh_out_supplemental, delimiter=',')
            
        
        ### Matching starts
        validated_tuples = []
        
        ##Read cleaveland parsed file and write column 1-3 in validated tiples for matching with target finder file
        
        cl_parsed = ('./cl_results/%s_cl_results_parsed.csv' % (lib[0]))
        #cl_parsed = ('./cl_results/%s_cl_results_parsed_test.csv' % (lib))##Test line
        print ('\nThe file %s is being joined with TargetFinder Results' % (cl_parsed))
        with open(cl_parsed) as fh1:
            fh1.readline()
            csv_reader = csv.reader(fh1)
            for row in csv_reader:
                validated_tuples.append(tuple(row[0:3]))
        print('The total number of entries i.e tuples in %s are %s' % (cl_parsed, len(validated_tuples)))

        
        ##Read TF file and compare with the column 1-3 in 'Validates_tuples' list 
        TF_parsed_file = './scoring/TF_res_parsed.csv'
        TF_match_list = []####The function of list to hold matching target finder entries which will be used later to join TF and CL entries
        with open(TF_parsed_file) as fh2:
            fh2.readline()
            csv_reader = csv.reader(fh2)
            match_count=0
            for row in csv_reader:
                if tuple(row[0:3]) in validated_tuples:
                    #csv_out.writerow(row[0:6])###+ add cleavage information
                    TF_match_list.append(row[0:5])                    
                    #csv_out_sup.writerow(row[0:6])            
                    match_count+=1
                else:
                    pass
        print('Total number of entries joined with target finder final results:%s\n' % (match_count))
#        print('Please use the output files %s and supplemental file %s  for further analysis' % (matched_file,matched_supp ))

        with open(cl_parsed) as fh1:###Here we join TF and CL entries to obtain miR and tar sequence
            fh1.readline()
            csv_reader = csv.reader(fh1)
            for ent in csv_reader:
                key = tuple(ent[0:3])
                #print(ent[0],ent[1],ent[2],ent[3])
                for i in TF_match_list:
                    #print(ent[0],ent[1],ent[2],ent[3])
                    #print (i)
                    if tuple(i[0:3]) == key:
                        fh_output.write('%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],ent[1],ent[2],ent[3],ent[4],i[3],i[4]))
                    else:
                        pass
                    
                    
        
    #Now combine the validated out results from all the libraries so that we can collectively remove redundant files in the next step
    
    validated_comb = './scoring/validated_all_lib_comb'
    validated_out = open (validated_comb ,'w')
       
    fls = glob.glob(r'./scoring/*validated_out')
    print ('Combining all the matched aka validated files from CL3 results of different PARE libraries to single file before removing redundant entries\n')
    for x in fls:
         CLfile = open(x, 'r')
         data = CLfile.read()
         CLfile.close()
         validated_out.write(data)
    validated_out.close()
    
#    validated_out.write('%s' % (results))
#    validated_out.close()
#    
    return validated_comb

def UniqScoringInp(validated_comb):
    
    ## 1. Remove redundant entries-with all the attributes same i.e DNA, chr#, cleavage site, entry name, tar, miRNA
    ## 2. Write different scoring inputs for final scoring steps i.e sco_inp_21, sco_inp_22 etc 
    
    ##Reading parsed output having redundant entries
    fh_in=open('./scoring/validated_all_lib_comb', 'r') 
    fh_in.readline()##Waste header line##There is no header
    parsed_in=csv.reader(fh_in, delimiter=',')
    
    ## Output file
    fh_output=open('./scoring/scoring_input_extend', 'w') 
    csv_out=csv.writer(fh_output, delimiter=',')
    
    fh_output2=open('./scoring/scoring_input', 'w') 
    #csv_out2=csv.writer(fh_output, delimiter=',')
    ##Method:3- To find unique entries from the 'Matched i.e validated entries from cleaveland output'
    
    added_keys=set()## A set to store first 3 elements from input file: miRNA-DNA, chr# and cleavage site and than use it to compare further entries in file
    
    parsed_out_count=0## To keep count of unique entries
    print('\nRemoving redundant entries aka more than one instance of an miRNA (PARE Validated)')
    for row in parsed_in:
    #    print(row)
        genename = row[1].split('.')[0]##To avoid different variations of same gene to be counted as uniq
        lookup=tuple((row[0],genename,row[2]))##miR name + Target Gene+position of cleavage on gene, earlier miR sequence was used instead of miR name but that removes miR belonging to same family as they have same sequence
        if lookup not in added_keys:##That means such a tuple has not been recorded yet and is unique
            csv_out.writerow(row)###Complete information
            fh_output2.write('%s\t%s\t%s\t%s\t%s\n' % (row[0],row[5],row[1],row[3],row[4]))##Information relevent to connect with Jixans pipeline script 2
            #csv_out.writerow([row[0]]+[row[3]]+[row[4]])##MiRNA name, miRNA seq and target seq
            parsed_out_count+=1
    #        outlist.append(row)
            added_keys.add(lookup)## Once a new entry is found it is recorded so as to compare and neglect further entries
        else:
            pass
            
    print('The number of unique entries found and will be used for scoring:', parsed_out_count )
    print('\n->>Two results file generated at this step "scoring_input" and "scoring_input_extend"<<-')
    
    #print(outlist)
    #print('The total number of entries after removing duplicates:', len(outlist))
    
    fh_in.close()
    fh_output.close()
    
    ###############  sub-module: Segregate miRNAs on the basis of size before actual scoring - eliminates the need to run count_mir_len
    
    
    print ('\nClassifying results on basis of miR size | Size is calculated after gap adjustment')
    print('Dont Panic! the results file generated in last step are retained\n')
    
    
    ##The output files corresponding to the length of miRNA
    fh1_out=open('./scoring/scoring_inp_21', 'w')
    csv1=csv.writer(fh1_out, delimiter=',')
    fh2_out=open('./scoring/scoring_inp_22', 'w')
    csv2=csv.writer(fh2_out, delimiter=',')
    fh3_out=open('./scoring/scoring_inp_23', 'w')
    csv3=csv.writer(fh3_out, delimiter=',')
    fh4_out=open('./scoring/scoring_inp_24', 'w')
    csv4=csv.writer(fh4_out, delimiter=',')
    
    
    ##The input file to be read and filtered according to their length
    fh_in1=open('./scoring/scoring_input_extend', 'r')
    csv_in1=csv.reader(fh_in1, delimiter=',')
    
    ## Read the file line by line and do operation
    for i in csv_in1:
        a=i[5].count('-')##Count the bulges and gap    
        original_len=int(len(i[5])-a)
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
            print('An miRNA %s: %s of length %s was found, it normal' % (i[0][1:],i[5],original_len))
            #print(i[0])
            pass
        
    print('\nPlease use the output file of required miRNA length from "scoring" folder for your analysis\n')
    
    ##Close all the files of this module
    fh1_out.close()
    fh2_out.close()
    fh3_out.close()
    fh4_out.close()
    fh_in1.close()


##- - - -- - - - - - - MAIN - - - - - - - - - - - - -
##Clean FASTA headers
trans_file_clean = CleanHeader(trascriptome_file)
miRNA_file_clean = CleanHeader(miRNA_file)

##Get PARE data
con = ConnectToDB()
libs=TagAbundanceFile(con)

##Cleaveland analysis, index should be made before as not a part of script
mapdd2trans(libs)
CleavelandAnalysis(libs, trans_file_clean)

##Processing
CleavelandParse(libs)
TF_parsed_file = TargetFinderParse('./target_finder_results')
validated_comb = MatchCLandTF(libs,TF_parsed_file)
UniqScoringInp(validated_comb)



'''
Created on Jun 2, 2012

@author: atul
'''

