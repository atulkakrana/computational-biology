#!/usr/bin/python

# - - - - - H E A D E R - - - - - - - - - - - - - - - - - -
 
    #1.Parsing
    #2.Scoring
    #    1. Create a nested list to record feature(mat,wob, gap,bulge, mis) at every position of sequence
    #    2. Scoring loops
    #    3. Changing scoring orientation in regarding to miRNA and calculate percentages
    #3.Histogram


# - - - - - F U N C T I O N S - - - - - - - - - - - - - - - 
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
#import matplotlib.table as table

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - U S E R    V A R I A B L E S  - - - - - - - - -
##Config

sco_inp = str(sys.argv[1])
#sco_inp = 'scoring_inp_21_jixian_2'##Scoring input file


## Other Filenames - no need to change

res_file = ('%s_result.tsv' % (sco_inp))##Results file
sup_file = ('%s_res_sup.tsv' % (sco_inp))##Supplemental results file
comp_file = ('%s_res_comp.tsv' % (sco_inp))##Complete results file
plot_file = ('%s_figure.png' % (sco_inp))##Plot results file

## Position variables for loop only
mat=0##match count
mis=0##mismatch count
gap=0##gap count
bul=0##bulge count
wob=0##wobble count
#read=0 ##reads the miRNa and target position for matching from start till end, not in the orientation in which final miRNA position will be reported
gpos=0 ## to keep track of positions from 1-22 relative to miRNA i.e -1miRNA=poscount1
poslists= list()## To make a nested list for number of base pairs and at every bp values for Mismatch, Wob, Bul, Match, gap is maintained

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - G L O B A L  V A R I A B L E S  - - - - - - - -
## Initialize feature count
gentrycount=0## Entry being read from the parsed file or you can also say line being read
poscount=0 ## to keep track of positions from 1-22 relative to miRNA i.e -1miRNA=poscount1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


################################ MODULE 1 ##################################
##1 : Module to count the number of entries
##Change the INPUT SCORING FILENAME here:@@@@@@@@@@
fh=open(sco_inp, 'r')
entries=fh.read()
#print ('The entries are:',entries)
entcount=int(entries.count('>'))
print('Total entries in file are:',entcount)
fh.close()

############################### MODULE2 ####################################
##2 : Module to find the longest miRNA length in file and use that length to define a scope of poslists matrix

alist=list()## an empty list of miRNA just for this module
##Change the INPUT SCORING FILENAME here:@@@@@@@@@@
fh=open(sco_inp, 'r')
## INDIVIDUAL LOOP STARTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1

for ent in range (entcount):
    anentry=fh.readline()
#    print('This is single entry:',anentry)
    anentry_strp=anentry.strip()
#    print('Strip entry:', anentry_strp)
    anentry_splt=anentry_strp.split(',')
#    print('An entry split', anentry_splt)
    
    miR=anentry_splt[2]
#    print(miR)
#    print(len(miR))
   
#    miRlen=len(spltd[anent])
    alist.append(len(miR))   
    
##END OF LOOP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~1

#print(alist)
maxlen=max(alist)##longest miRNA in file, this base pair length will be used to create poslists matrix for each bp 
print('Max miRNA length:',maxlen)

################################ MODULE 3 #####################################
##3: Module to create a poslist matrix to be filled later equal to length of longest miRNA

##START OF INDIVIDUAL LOOP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2
##For making scoring list of length equal to maximum characters encountered in longest entry

for x in range(maxlen):      
    poslists.append([0,0,0,0,0]) ##Make Position variables for each bp  of miRNA
    
## END OF LOOP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~2
    
print('The size of poslists when created:',len(poslists))
fh.close() ## If not closed here than interferes with next LOOP1

################################# MODULE 4 ####################################
##4: Module to extract single entry from parsed file and feed miRNA and target to scoring loop below. 
###A lot of information like length of miRNA and target and their sequences can be output from this loop

##Write supplemental info from LOOP1 to a file, miRNA, length of miRNA, Target, length of target, file opened in Append mode
print('\n``````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````')
print('WARNING: Re-running the script on same file appends supplemental results in the "results_supplemental" file')
print('``````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````\n')

fh_csv=open(sup_file, 'w')
sup_result=csv.writer(fh_csv, delimiter='\t')
##Header for Supplemental File
sup_result.writerow(['EntryCount', entcount,'','',''])
sup_result.writerow(['Entry','miRNA',' miRLen', 'Target',' TargetLen'])

#Open the parsed file to read
##Change the INPUT SCORING FILENAME here:@@@@@@@@@@
fh=open(sco_inp, 'r')
                   
##START OF MAIN LOOP~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~3
##THE 'SINGLE ENTRY READING FROM PARSED FILE'  LOOP STARTS FROM HERE

for line in fh:
    line_strp = line.strip('\n')
    ent = line_strp.split(',')
#    check_tar = ent[1]
    pre_tar=list(ent[1])
    pre_tar.reverse()
    tar=''.join(pre_tar)
    tarlen = len(tar)
    
#    check_mi =ent[2]
    pre_mi = list(ent[2])
    pre_mi.reverse()
    mi=''.join(pre_mi)
    milen = len(mi)
#    print(tar,mi)
    
    ##maintaining the number of entries processed
    gentrycount+=1  

    ##Writing the required values to supplemental file
    sup_result.writerow([gentrycount]+[mi]+[milen]+[tar]+[tarlen])

    ##check for length mismatch
    if tarlen!=milen:
        print(' miRNA and target alignment length does not matches', gentrycount)
        break## NEEDS FIX-SHOULD COME OUT OF SCRIPT
    else:
        pos=0## Initialized zero before scoring loop and incremented only at end of Scoring loop so that doesn't get affected by skipping of frame by GAP in scoring loop
        pass
    
################################# MODULE  5 ########################################
##5: Module to score the single entry fed from GIANT LOOP 3
   
##LOOP4-INNER STARTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~4
## Under GIANT LOOP3 for reading and matching through all miRNA positions for single entry/alignment
   
    for bp in range(len(mi)):
        if mi[bp]=='-':
            poslists[pos][4]+=1##GAP is counted at this point as frame is shifted to next and scoring is done on same position as after skipping it becomes to score next bp      
            ##To score gap as in old version change above to poslists[pos-1][4]+=1
            continue
    #    elif mi_strp[bp]=='':# [read=1] equals to len(mi)
    #        break
        else:
    #        print('The base pairs read after skipping', bp)
            pass
                           
        
        ##Complementing the target for matching with miRNA in next if loop
        tar_c=tar[bp].translate(str.maketrans("AUGC","UACG"))
#        print(tar,tar_c)
        
        ##SCORING LOOP6 STARTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                       
        if mi[bp]==tar_c:
            mat+=1
    #        print ('Match',read)
            poslists[pos][3]+=1##See above for keys to value location in poslists matrix
            
        ## If a mismatch is found than it is checked for reason in following order:: wobble (G:U) pairing,  gap, bulge, base mismatch     
        elif mi[bp]!=tar_c:
            
            ##NESTED LOOP7 STARTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            ##Wobble testing
            if  mi[bp]=='G' and tar[bp]=='U':
                wob+=1
                poslists[pos][1]+=1##See above for keys to value location in poslists matrix
    #            print('Wobble:',read)
            elif mi[bp]=='U' and tar[bp]=='G':
                wob+=1
                poslists[pos][1]+=1
    #            print('Wobble:', read) 
        
            ## The BULGE  in miRNA is defined by the gap'-' in target as target lacks matching causing a bulge to appear in miRNA 
            elif tar[bp]=='-':
                bul+=1
                poslists[pos][2]+=1##See above for keys to value location in poslists matrix
    #            print('Bulge:',read)
               
            ## Score MISMATCH
            else:
                mis+=1
                poslists[pos][0]+=1##See above for keys to value location in poslists matrix
    #            print('Mismatch:',read)      
            ##Nested LOOP 7 ENDS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
            
        else:
            pass
        ## SCORING LOOP 6 ENDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        pos+=1## Incremented at end so not effected by skipping, 0 at initial and after scoring only incremented
        
#        print(read)##Was used just to trouble shoot the loop leaking problem, i.e just first bp was being read and filled in all matrix positions

print("Length of poslits after scoring module:", len(poslists), "| Must be equal to 'Max mi RNA Length' :", maxlen)

#Closing file handles opened before LOOP1 in MODULE 4
fh.close()
fh_csv.close()


################################### MODULE 6 #####################################
##6. Module to writing the Poslists in required orientation

#Open file to write in csv format
outfile  = open(res_file, 'w')
csv_writer = csv.writer(outfile, delimiter='\t')

## Module to calculate percentage of appearance of Match, Wobble, gap, Bulge and BP mismatch by going through each bp in POSLISTS
poslist_reduced_len=0
for location in range(len(poslists)):
    a=poslists[location]## a takes the  first base pair values to convert into percentage
##Taking the values of mismatch, wobble, bulge, match and gap into 'pos' variables to calculate sum
    pos1=float(a[0])    
    pos2=float(a[1])    
    pos3=float(a[2])    
    pos4=float(a[3])   
    pos5=float(a[4])    
    summed=pos1+pos2+pos3+pos4+pos5
    #    print(summed)
## Loop to delete extra rows that were generated due to extra length of miRNA by GAPS
    if summed==0: ##if there is no feature scored than its the extra position and need not to be included in result 
#        print('This row has no element')
        pass
    else:
        perc=[x/entcount*100 for x in a]
        csv_writer.writerow(perc)
        poslist_reduced_len+=1
print('The final length of poslist is reduced to the actual length of miRNA:', poslist_reduced_len)
outfile.close()


#################################### MODULE 7 ####################################
##7. Module to append the results and supplemental info (captured in loop 1)
    
###Append supplemental info to csv_writer
### The variable 'tf' stands for 'total file'
tf = open(comp_file,'a')
src1 = open(sup_file,'r')
src2 = open(res_file, 'r')
a=src1.read()
b=src2.read()
#Write in order results and than supplemental files
tf.write(b)
tf.write('\n')
tf.write(a)
#Close the files for this module
tf.close()
src1.close()
src2.close()
   
    ######################################## MODULE 8 ###############################

##7:  Module to draw stacked histogram plots   

## Make empty list for every feature
##These are five features need to stacked over each other for different position from 1-22 on X axis
mis_list=list()
wob_list=list()
bul_list=list()
mat_list=list()
gap_list=list()

##Read file Results file@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
fh2_csv=open(res_file, 'r')
csv_reader=csv.reader(fh2_csv, delimiter='\t')

N=0 ## actual number of miRNA base pairs i.e length of miR will be counted by this variable, every row in input file represents features from each bp and thus rows are counted
for row in csv_reader:
    mismatch=(row[0])
    mis_list.append(float(mismatch))
    
    wobble=(row[1])
    wob_list.append(float(wobble))
    
    bulge=(row[2])
    bul_list.append(float(bulge))
    
    match=(row[3])
    mat_list.append(float(match))
    
    gapped=(row[4])
    gap_list.append(float(gapped))
    
    ##Count the number of positions (rows in input file for plotting graph)
    N+=1
    
ind=np.arange(N)
width = 0.30

##empty lists to be used for bottom function
bul_list_bottom = []
mat_list_bottom = []
gap_list_bottom = []

for i in range(N):
    bul_list_bottom.append(mis_list[i]+wob_list[i])
    mat_list_bottom.append(mis_list[i]+wob_list[i]+bul_list[i])
    gap_list_bottom.append(mis_list[i]+wob_list[i]+bul_list[i]+mat_list[i])

##plotting variables
p1=plt.bar(ind, mis_list, width, color = 'r',)
p2=plt.bar(ind, wob_list, width, color = 'm', bottom=mis_list)
p3=plt.bar(ind, bul_list, width, color = 'b', bottom=bul_list_bottom)
p4=plt.bar(ind, mat_list, width, color = 'w', bottom=mat_list_bottom)
p5=plt.bar(ind, gap_list, width, color = 'y', bottom=gap_list_bottom)


plt.ylabel('Percentage (total miRNAs:%s)' % entcount, fontproperties=font_manager.FontProperties(size=10))
plt.xlabel('miRNA position', fontproperties=font_manager.FontProperties(size=10))
plt.title('miRNA-target feature scoring')
plt.xticks(np.arange(22), np.arange(1,23) )
plt.yticks(np.arange(0,121,5))
plt.legend((p1[0], p2[0], p3[0], p4[0],p5[0]), ('MISMATCH','WOBBLE','BULGE','MATCH','GAP'), loc=1, prop=font_manager.FontProperties(size=7))


#addtable(   )
#table(rowLoc='left', colLabels=None, colColours=None, colLoc='center', loc='bottom', bbox=None)

##give name of OUTPUT GRAPH@@@@@@@@@@@@@@@@@@@@@@@@@@@
plt.savefig(plot_file, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)

#plt.show()     

# Close the file opened to draw histograms in this module
fh2_csv.close()
print('\nThree files generated are: "results", "results_supplemental" and "results_figure"')
print('\n````````````````````````````````````````````````````````End of script````````````````````````````````````````````````````````````')
            
            
###Remove intermediate files
# os.remove('Path/To/File.ext')


'''
Created on Apr 13, 2012

@author: setu
'''
