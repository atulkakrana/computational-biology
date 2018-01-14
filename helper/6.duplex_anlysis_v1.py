##script written by Atul


def uniq_mir(fh_in):
##Take miRNA names from the scoring input file, as they will be repeated so just take uniq names
##The rational is that these miRNA are validated by PARE data and we need to see what is the length of their duplex
##So, we just need the name of miRNA.
    mir_set = set()
    for i in fh_in:
        mir = i.split(',')[0][1:].lower()
#        print(mir)
        mir_set.add(mir)
#    print(len(mir_set))
    return mir_set

#This module finds the mir_minor length of validated miRNAS from the mir-miR* list and compares it with miR length, 
#if length is not equal than its counted as asymmetry 
def find_dup_len(mir_set):
    
    ##Read and compare with miRNAs validated from PARE data 
    acount = 0 ##Counts the number of miRNAs found in mir-Mir* list
    asym_count = 0 ##Counts the number of miRNA duplexes that show asymmetry i.e mir != mir*
    miR_min_NA_count = 0 ## To keep count of miRNA duplexes where miR* was not found in input mir-mir* list
    for ent in fh_in2: 
        ent = ent.strip('\n')       
        ent = ent.split(',')
        mir_name = ent[0].lower()
        if mir_name in mir_set:
            mir = ent[1]
#            print(mir_name)
            mir_len = len(mir)
            mir_min = ent[2]## miR minor
            mir_min_len = len(mir_min)
#            print(mir,len(mir), mir_min, len(mir_min))
#            print(mir_name)
            acount +=1
            if mir_min == 'NO PREDICTION':
#                print ('mir* not found')
                miR_min_NA_count +=1
            elif mir_min_len != mir_len:
                asym_count +=1
    perc_sym = 100*asym_count/(acount-miR_min_NA_count)
    print('Total entries:%s | Asymmetric entries:%s | Percentage Asymmetry: %s' % (acount,asym_count,perc_sym))
    fh_out.write('%s\t%s\t%s\t%s\n' % (acount,asym_count,miR_min_NA_count,perc_sym))
    return asym_count/acount*100
            
            
            
            
        
    

##For Common entries to count nucleotides in miR


##MAIN

##Make a list of Xmers you want to analyze
mer_list = [21,22]

fh_out = open('duplex_info', 'w')
##Header 
fh_out.write('Entries\tAsym_Entries\tMissing_mir*_ent\tPerc_asym\n')

for i in mer_list:
    ##Input miR-miR*.csv and scoring_22mer
    ##Scoring Input file
    file_in = ('scoring_inp_%s' % (i))
    fh_in = open(file_in , 'r')

    ##File with miR-mir* information
    fh_in2 = open('/home/atul/At_miR_miR*.csv', 'r')
    
    mir_set = uniq_mir(fh_in)
    asym_perc = find_dup_len(mir_set)
    


fh_out.close()
fh_in.close()
fh_in2.close()


'''
Created on Jul 11, 2012

@author: atul
'''
