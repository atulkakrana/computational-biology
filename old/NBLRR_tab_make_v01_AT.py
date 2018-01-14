##This script is to make table of genes targetted by several miR(NB_LRR targetting) by using parsed results from NBLRRv0.1 (below)

#>mtr-miR1507,AT3G14470.1,7,847-868,GAUCUACUACAUACCUUGCUCC,CUAGAUGAUCUGUGGAAUGAAA
#>mtr-miR1507,AT5G48620.1,7,1744-1765,GAUCUACUACAUACCUUGCUCC,CUCGAUGAUGUAUGGAAAAAAG

##Fix- 0 is not written for no target instead - is written

##Config

inp_file = 'mt_all_targetsparsed_v2.csv'
out_file = 'mt_all_NBLRR_TF_tab'


def TableMaker(file):
    
    ##Open file
    fh_in = open(file, 'r')
    fh_in.readline()##Waste the header line
    
    inp_file = fh_in.readlines()## Complete file is read because it will used may times and now we do not have to open and close again
    ##Create a set of unique targets
    tar_set = set()## A set to store unique targets 
    mir_set = set() ## A set to store unique miRs  
    for ent in inp_file:
        print(ent)
        tar = ent.split(',')[1][0:-2]###Gene models were not considered
        print (tar)
        mir = ent.split(',')[0][1:]
#        print(tar)
        tar_set.add(tar)
        mir_set.add(mir)        
    print('The total number of NBLRR targets are: %s\n' % (len(tar_set)))
    print('The order of miR in table is : %s' % (mir_set))
    
    
    ##Counters
    
    records = ([['-' for col in range(int(len(mir_set))+1)] for row in range(len(tar_set))])## This will keep records of all targets and the score of miRs targetting them
    tar_count = 0## This will keep count of uniq targets and also track the specified position on where to store scores for multiple miRs in records
    
    
    ent_count = 0## This will just keep count of entries processed, keep in mind that it is the entries multiplied by number of targets analyzed
    match_count = 0 ###must be equal to entries in inp_file
    
    
    for tar in tar_set:##For every uniq 'target' in file for which a single entry will be made
#        print(tar)

        for ent in inp_file:## Go through every entry and fetch entries matching target taken from set
#            print (ent.split(',')[1][0:-2])
            if tar == ent.split(',')[1][0:-2]:##If entry has target than find miR and recoed its score at specific position
                match_count += 1## Increase that match count
                ent = ent.split(',')                
                mir = ent[0][1:]
                score = ent[2]
                records[tar_count][0] = tar
                
                mir_pos = 1## This will add scores to 'target tuple' at mirna specific position, ,mir_pos =0 is for miRNA Name
                for miR in mir_set:
#                    print(miR)## Checked the pattern it follows always the same
                    if miR == mir:##If miRNA from uniq mirset is same as mir in the entry than record the score in 'records' at the same position as in mir_set
                        if records[tar_count][mir_pos] == '-': ## If the position is still not filled
                            records[tar_count][mir_pos] = score
                        elif  records[tar_count][mir_pos] != '-' and records[tar_count][mir_pos] > score:##If position is filled with score and same miR targets same 'target' than choose which one lower and fill
                            records[tar_count][mir_pos] = score
                        else:
                            pass       
                        mir_pos += 1## The position in target tuple is increased after every mirna analyzed
                    else:
                        mir_pos +=1 ## The position in target tuple is increased after every miRNA analyzed even if doent matched with mirna from entry
            
            ent_count += 1    
        tar_count +=1        
#                print(ent)
    print ('Targets processed = %s | Entries processed : %s | Matches found : %s' % (tar_count, ent_count, match_count))
            
    
    return records,mir_set


def WriteOutput(records,mir_set):
    alist = list(mir_set)
    fh_out = open (out_file, 'w')
    fh_out.write('Target\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (alist[0],alist[1],alist[2],alist[3],alist[4],alist[5],alist[6]))
    #fh_out.write('Target\tmiR2118c\tmiR5213\tmiR5213*\tmiR2118a\tmiR2118b\tmiR1507*\tmiR1507\n')
    for tuple in records:
        print(tuple)
        fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (tuple[0],tuple[1],tuple[2],tuple[3],tuple[4],tuple[5],tuple[6],tuple[7]))
    
    fh_out.close()
        
    

###Main
records,mir_set = TableMaker(inp_file)
WriteOutput(records,mir_set)





'''
Created on Aug 8, 2012

@author: atul
'''
