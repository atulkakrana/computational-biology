#!/usr/local/bin/python3

import os
import glob
import sys
import difflib

###Chech PHAS loci for redundancy

res_folder = './files'
p_val = '1e-07'

def RemoveRedundant(res_folder, p_val):
    
    ##Prepare output file
    outfile = 'Final_PHAS_Loci.csv'
    fh_out = open(outfile,'w')
    fh_out.write('Name\tp-val\tChr\tStart\tEnd\tStrand\n')
    
    ##Read files
    fls = glob.glob(r'./files/*YES.by.PARE.txt')
    print (fls)
    
    main_dict = {} ## declare empty dictionary
    anum = 1 ## To name PHAS loci
    for afile in fls: ###
        
        tmp_dict = {}## Dictionary to store values for one file - recycled after every file
        
        fh_in = open(afile, 'r')
        
        
        ###First Instance - Fill up the dictionary with entries from filrst file
        if not main_dict.values(): ## First file will populate dictionary
            lines = [i for i in fh_in if i[:-1]] ## Remove empty lines one liner
            for ent in lines:
                #print (ent)
                ent_splt = ent.strip('\n').split('\t')
                #print (ent_splt)
                #break
                if float(ent_splt[1]) <= float(p_val): ###Less than specified cut-off
                    key = '%s-%s-%s' % (ent_splt[2],ent_splt[3],ent_splt[4])
                    value  = list(range(int(ent_splt[3]),int(ent_splt[4])))
                    #print ('Key:%s and Val%s' % (key,value))
                    ###Check for redundancy, if above cutoff add to dict and write ent  to file
                    tmp_dict[key] = value
                    print('Phas%s\t%s\t%s\t%s\t%s\t%s\n' % (anum,ent_splt[1],ent_splt[2],ent_splt[3],ent_splt[4],ent_splt[6]))
                    fh_out.write('Phas%s\t%s\t%s\t%s\t%s\t%s\n' % (anum,ent_splt[1],ent_splt[2],ent_splt[3],ent_splt[4],ent_splt[6]))
                    anum += 1
                    
                else:
                    print ('All entries with specified cutoff analyzed')
                    break
            
            main_dict.update(tmp_dict) ## Update the main dict
            
            
        ##Second Instance - Match the entries of new file with dictioary and add new one
        else: ##Dictionary has key and value populated from first file and now remove redundancy
            lines = [i for i in fh_in if i[:-1]] ## Remove empty lines one liner
            for ent in lines:
                ent_splt = ent.strip('\n').split('\t')
                #print (ent_splt)
                
                if float(ent_splt[1]) <= float(p_val): ###Less than specified cut-off
                    key = '%s-%s-%s' % (ent_splt[2],ent_splt[3],ent_splt[4])
                    value  = list(range(int(ent_splt[3]),int(ent_splt[4])))
                    for i in main_dict.values():
                        #print ('Dictionary entry being matched:%s \n' % (i))
                        sm=difflib.SequenceMatcher(None,i,value)
                        if sm.ratio() >= 0.6:
                            print ('Adding new key')
                            tmp_dict[key]=value ## Life = one file
                            print('Phas%s\t%s\t%s\t%s\t%s\t%s\n' % (anum,ent_splt[1],ent_splt[2],ent_splt[3],ent_splt[4],ent_splt[6]))
                            fh_out.write('Phas%s\t%s\t%s\t%s\t%s\t%s\n' % (anum,ent_splt[1],ent_splt[2],ent_splt[3],ent_splt[4],ent_splt[6]))
                            anum += 1
                        
                        else:
                            print('Redundant')
                            pass
                else:
                    print ('All entries with specified cutoff analyzed')
                    break
            
            main_dict.update(tmp_dict) ### Update the main dict
    



def main():
    RemoveRedundant(res_folder,p_val)
    

if __name__ == '__main__':
    main()
    sys.exit()

        
        
    
