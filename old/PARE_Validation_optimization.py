#!/usr/bin/python
##Author: Atul Kakrana | kakrana@udel.edu

##This script is to filter the results from Cleaveland 2 PARE validation pipeline. After the Jixians pipeline  you get 
##the  'Output_tasi_target_LINE_soy_jiaxian_validated_known_results_soy_full_genome_REAL_targets_7' kind of file and to 
##filter the miRNA targets on the basis PARE libraries using three parameters 'a', 'b' and 'c' (see config section below)

#Usage: change the file_name in config section for the file from Jixians 'miR_target_PARE_Phase_V3_Cleaveland-3_Atul.pl' script

# - - - - - - - - - - - - -  H E A D E R - - - - - - - - - - - - - - - - - -

##Config
file_name = ''
a=6 ## PARE score cutoff (Less the better), default = 5 (stringent)
b=5 ##Small window value cutoff(More the better), default = 5 stringent
c=0.50 ## Ws/Wl ratio cutoff(More the stringent), default = 0.5 (stringent)
d=10##Total number of libraries in Pare validated file i.e (ws+wl)/2 columns in input file

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - -U S E R    V A R I A B L E S  - - - - - - - - - - - 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - F U N C T I O N S - - - - - - - - - - - - - -  
##Signature
print("\nScript by /-\ ^|^ |_| |_ for Blake's Lab\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - -M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - -

fh_in= open(file_name, 'r')
fh_out = open('PARE_validated_miRNAs', 'w')

firstline = fh_in.readline()## Skip first line
fh_out.write('%s\t%s\t%s\n' % (firstline.strip('\n'),'Sum Ws','Sum Wl'))## Print in output file

ws_sum_list=[]## To keep track of sum and find max value in end for printing
wl_sum_list=[]


line_count = 0
filtered_count = 0
for line in fh_in:
#    print(line)
    line_strp=line.strip('\n')
    ent = line_strp.split('\t')
    mir = ent[0]
    score = float(ent[5])
    lib_values = ent[14:]
#    print(lib_values)
    
    ## Filter 1 : PARE score <= 5
    if score <= a: ## 'a' from config
#        print ('Score is less than 5:', score)
        acnt = 0## count to index the alternative ws and wl values from lib_values
        ws_sum = 0## Sum of values in small window
        wl_sum = 0## Sum of values in large window
        ratio = 0## Ratio of the ws/Wl        
        for i in range (d): # Scoring the library values
            ws_sum += int(lib_values[acnt])# Sum of small window
            wl_sum += int(lib_values[acnt+1])# Sum of large window
            acnt+=2
#        print('Ws sum:', ws_sum, '| Wl Sum:', wl_sum, '\n')  
  
        ## Filter 2: Ws_sum >= 5
        if ws_sum >= b:## 'b' from config
#            print('ws_sum is more than 5:', ws_sum)
            ratio = float(ws_sum/wl_sum)
        
        ## Filter 3: Ws/Wl i.e ratio    
        if ratio >= c:## 'c' from config
#            print ('Filter 3- Ws Sum:', ws_sum,'Wl Sum:', wl_sum, 'Ratio:',ratio,'\n')
#            fh_out.write(line)
            fh_out.write('%s\t%d\t%d\n' % (line_strp,ws_sum,wl_sum))
            ws_sum_list.append(ws_sum)
            wl_sum_list.append(wl_sum)
            filtered_count += 1
            
            
                    
            
            
    line_count += 1


print('\nTotal number of entries analyzed:', line_count)
print('The max sum from Small window:', max(ws_sum_list))
print('The max sum from Large window:', max(wl_sum_list))
print('Total number of entries passed filter:', filtered_count,'\n')
print("Check the output in 'PARE_validated_miRNAs'")  
            
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
