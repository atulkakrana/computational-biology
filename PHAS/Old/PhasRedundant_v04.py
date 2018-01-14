#!/usr/local/bin/python3
###Goes through 'valid.cluster.XX' files in a folder and remove redundant PHAS loci

import os
import glob
import sys
import difflib
import time

###Chech PHAS loci for redundancy

res_folder = './DEG_VALID'
p_val = '1e-07'
genome = 'RICE_MSU7'##Name for PHAS cluster naming

def RemoveRedundant(res_folder, p_val):
    
    fh_out = open('PHASRedundant.log','w')
    ##Read files
    fls = glob.glob(r'./%s/*YES.by.PARE.txt' % (res_folder))
    print (fls,'\n')
    print ('Total files to analyze: %s' % (len(fls)))
    
    main_dict = {} ## declare empty dictionary
    anum = 1 ## To name PHAS loci
    for afile in fls: ###
        print ('**Analyzing file: %s' % (afile))
        tmp_dict = {}## Dictionary to store values for one file - recycled after every file
        neg_list = []## List to store keys that needs to be removed
        fh_in = open(afile, 'r')
        alib = afile.split('/')[-1].split('.')[0]###
        
        ###First Instance - Fill up the dictionary with entries from filrst file
        if not main_dict.values(): ## First file will populate dictionary
            lines = [i for i in fh_in if i[:-1]] ## Remove empty lines one liner
            for ent in lines:
                #print (ent)
                ent_splt = ent.strip('\n').split('\t')
                #print (ent_splt)
                #break
                #print (ent)
                if float(ent_splt[1]) <= float(p_val): ###Less than specified cut-off
                    key = '%s-%s-%s' % (ent_splt[2],ent_splt[3],ent_splt[4])###Chrid, strat and end makes key
                    start = ent_splt[3]
                    end = int(ent_splt[4])+1 ### 1 added because when opening a range using start and end, end number is not included in range
                    value  = (list(range(int(ent_splt[2]+str(start)),int(ent_splt[2]+str(end)))),ent_splt,alib) ### The range has chromosme number in begining to make range unique to chromosme
                    #print ('Key:%s and Val%s' % (key,value))
                    ###Check for redundancy, if above cutoff add to dict and write ent  to file
                    tmp_dict[key] = value
                    
                else:
                    print ('***All entries with specified cutoff analyzed for - %s' % (afile))
                    break
            
            main_dict.update(tmp_dict) ## Update the main dict
            
            
        ##Second and further Instance - Match the entries of new file with dictioary and add new one#############
        else: ##Dictionary has key and value populated from first file and now remove redundancy
            lines = [i for i in fh_in if i[:-1]] ## Remove empty lines one liner
            for ent in lines:
                ratiodict= {}###dict to hold ratios of comaprarision of this entry with all in dictionary
                ent_splt = ent.strip('\n').split('\t')
    
                ###Compare with dict entries and get ratio
                if float(ent_splt[1]) <= float(p_val): ###Less than specified cut-off
                    start = ent_splt[3]
                    end = int(ent_splt[4])+1 ### 1 added because when opening a range using start and end, end number is not included in range
                    value  = (list(range(int(ent_splt[2]+str(start)),int(ent_splt[2]+str(end)))),ent_splt,alib)####Ent_splt[2] is chr_id
                    
                    for i in main_dict.values():##Compare with all dictionary values
                        #print (i, main_dict[i])
                        sm=difflib.SequenceMatcher(None,i[0],value[0])
                        ratiodict[str(i)]=round(sm.ratio(),2)###Make a dict of main dict entries and their comparision ratio with current new entry
                        #ratiolist.append(round(sm.ratio(),2))
                else:
                    print ('***All entries with specified cutoff analyzed for - %s' % (afile))
                    break
                
                #Decide if entry is different enough to be added
                mainkey = max(ratiodict,key=ratiodict.get)###Key from main_dict with max comparable ratio for current entry
                #print(mainkey)
                mainkey_decode= mainkey.split(']')[0].split('[')[1].split(',')
                mainkey_remade = '%s-%s-%s' % (mainkey_decode[0].strip()[0],mainkey_decode[0].strip()[1:],mainkey_decode[-1].strip()[1:])
                
                #print(mainkey.split('[')[2].split(','))
                #mainkey_decode = mainkey.split('[')[2].split(',')
                #mainkey_remade = '%s-%s-%s' % (mainkey_decode[1][2],mainkey_decode[1][3],mainkey_decode[1][4])
                maxratio = ratiodict[mainkey]### Max ratio
                #print(mainkey_remade,maxratio)### If maxratio is zero same entry will appear again here
                print (ent,maxratio)
                if maxratio <= 0.35: ###Treat as new loci
                    #print ('Adding new key')
                    key = '%s-%s-%s' % (ent_splt[2],ent_splt[3],ent_splt[4])
                    tmp_dict[key]=value ## Life = one file
    
                    
                elif maxratio > 0.35: #### Choose the longest loci
                    #print('Selecting longest loci')
                    if len(mainkey_decode) < len(value[0]):### New Loci is longer
                        #print(value[0],mainkey_decode)
                        neg_list.append(mainkey_remade)
                        tmp_dict[key]=value
                    
                    else: ## The loci in dictionary is longer
                        pass
                else:
                    print('Redundant')
                    pass
                
            main_dict.update(tmp_dict) ### Update the main dict
            ########################Test####################
            #print ('Dictionary')
            #for i in main_dict.keys():
            #    print (i)
            #print ('Dictionary print finish')
            ################################################
            
            for akey in neg_list:
                #print (akey)
                try:
                    del main_dict[akey]
                    #print (akey, 'Key found in main dict and is being removed')
                except KeyError:
                    #print (akey, 'Key not found')
                    pass
                
            pass
        print ('**Number of current Phased loci: %s**\n' % (len(main_dict)))
        fh_out.write('**Number of current Phased loci: %s**\n' % (len(main_dict)))
    
    print ('Number of final phased loci: %s' % (len(main_dict)))
    fh_out.write('Number of final phased loci: %s' % (len(main_dict)))
    fh_out.close()
    return main_dict

def writer(main_dict):
        
    ##Prepare output file
    outfile1 = 'Final_PHAS_Loci.csv'
    fh_out1 = open(outfile1,'w')
    fh_out1.write('Name\tp-val\tChr\tStart\tEnd\tStrand\n')
    outfile2 = 'PHASLociID'
    fh_out2 = open(outfile2,'w')### For our Genome viewer, No header required
    
    anum = 1
    for value in main_dict.values():
        #print(value)
        #print('Phas-%s\t%s\t%s\t%s\t%s\t%s\n' % (anum,value[1][1],value[1][2],value[1][3],value[1][4],value[1][6]))
        fh_out1.write('Phas-%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (anum,value[1][1],value[1][2],value[1][3],value[1][4],value[1][6],value[2]))
        fh_out2.write('%s.%s.%s.%s\n' % (genome,value[1][2],value[1][3],value[1][4]))
        anum +=1
        
    fh_out1.close()
    fh_out2.close()
    return outfile1,outfile2

def main():
    main_dict = RemoveRedundant(res_folder,p_val)
    outfile1,outfile2 = writer(main_dict)
    

if __name__ == '__main__':
    start = time.time()
    main()
    print ('It took', round(time.time()-start,2), 'sec')
    print ('\nCheers script finished sucessfully\n')
    sys.exit()



##v01 -> v03
##Added functionality to compare between similar loci and select the longest one
##v03 -> v04
## Corrected bug with 'key' - Range was being made one integer less as last integer is never counted and with every run (file) there is change in key (one integer less) it is
###Therefore key from first run not matched woth earlier and all loci being retained

        
        
    
