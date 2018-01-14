#!/usr/bin/python

##Signature
print("\n\t\tScript by /-\ ^|^ |_| |_  for Blake\'s Lab'\n")

##This function calculates the numbers of 22 mers in a 500nt cluster for allthe chromosomes. You may want to change database settings in function 1 and length of miRNA you wish to study in Function 3 in query at line 58. The output files will be generated in the folder of script.

##Keep in Mind Atul
##1. This script can work on both tag position summary and tag position table, interchange the table names
##2. Change the minimum number of 22mer in cluster in in 'outputwrite' function. Currently, it is set to 10



import MySQLdb as mdb
import sys
import re
#import os


#Function 1: Connects to db
def connectToDB():
    try:
#        conn=mdb.connect(host='pikachu.dbi.udel.edu', user='kakrana', passwd='livetheday', db='reza')
        conn=mdb.connect(host='127.0.0.1', user='kakrana', passwd='livetheday', db='AT_sbsML_sRNA_v4')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0],e.args[1])
        sys.exit(1)

    return conn

##Function 2: This function calculates the total number of chromosomes and length of each chromosome
def dbrange(conn):
    bininfo = []##Empty list to store number of bins and last bin length
    cur = conn.cursor()
    cur.execute("SELECT MAX(chr_id) FROM tag_position_summary;")##Calculates the maximum number of chromosomes i.e total number of chromosomes
    maxchr= cur.fetchall()
#    print maxchr[0][0]

    for chromo in range (maxchr[0][0]):##For all the chromosomes
#        print chr_num          
        cur.execute("SELECT max(position) FROM tag_position_summary WHERE chr_id = %s;" % (chromo+1))##Get the length of each chromosome | chromo+1 because first loop as entry 0 so add one to make it chrmosome 1
        maxlen  = cur.fetchall()
#        print maxlen[0][0]
        bins=int(maxlen[0][0]//500)##List index were used as the maxlen is a tuple
        lastbin=int(maxlen[0][0]%500)##List index were used as the maxlen is a tuple
        
        bininfo.append((chromo+1,bins,lastbin))##chr num, total number of clusters, last cluster remaining length
#    print(bininfo)
    return bininfo
#    return bins
#    return lastbin



## Function 3: This functions counts the abundance of 22mers in bin
def abundance22merbin(conn, bininfo):
    abundance_list=[]
    for ent in bininfo:## i.e number of chromosomes as each chromosome has a tuple in bininfo 
        print 'Clustering and recording 22mers in Chromosome %s\n' % (ent[0])     
#        print ent
        pos=1##Variable to hold position on chromosome through the loop, set to 1 when next chromosome is analyzed
        cluster = 1
        for i in range(ent[1]):##ent[1] is number of bins## change range(ent[1]) after debugging for real analysis
#            print ent[0]##aka chromosome number
            cur = conn.cursor()
            cur.execute("select count(*) from tag_position_summary where chr_id= %d and length(tag) = 22 and (position between %s and %s);" % (ent[0], pos, pos+499))
            mircount = cur.fetchall()
#            print mircount[0][0]
#            abundance_list.append((ent[0],pos,pos+499,mircount))
            abundance_list.append(('Chr'+str(ent[0])+'_clust'+str(cluster), int(mircount[0][0])))##Chr_num+bin number as unique identifier
            pos+=500
            cluster+=1
        
        ###After the exact multiple bins were analyzed remainder sequence at end will be analyzed by this query
        finalbin_num = ent[1]+1
        finalbin_start = (ent[1]*500)+1##After the clusters of length 500 are analyzed,last cluster start position is calculated
        finalbin_end= (finalbin_start+ent[2])-1
        cur2 = conn.cursor()
        cur2.execute("select count(*) from tag_position_summary where chr_id= %s and (position between %s and %s);" % (ent[0], finalbin_start, finalbin_end))
        mircount_last = cur2.fetchall()
        abundance_list.append(('Chr'+str(ent[0])+'_clust'+str(finalbin_num), int(mircount_last[0][0])))##Chr_num_bin number as unique identifier

#        print abundance_list
    return abundance_list           
                        
##This function sorts the results chromosome wise from the 'abundance list' and records them into separate files 
 
def outputwrite(abundance_list):
    chr_num_re=re.compile('\d{1,2}')
    max_chr_num=re.findall(chr_num_re,abundance_list[-1][0])##Last chromosome number extracted from abundance list
#    print abundance_list[-1]
#    print max_chr_num[0]###Index used to get integer
    for i in range (1,int(max_chr_num[0])+1):##For number of chromosomes
        print 'The results are being sorted and written for Chromosome %s\n' % (i)
        fh_out=open('22mers_in_bin_chr%s' % (i), 'w')##Open one file for each chromosome
        fh_out.write('Identifier\tCluster\tFrequency\n')
        clust_num=1
        for ent in abundance_list:## Go through each entry of abundance list
#            print ent[0]        
            if re.search('^Chr%d_' % (i), ent[0]):##If entry has Chromosome number than write to file
#                print 'yes'
                if ent[1]> 10: ##CHANGE MINIMUM 22mer value HERE
                    fh_out.write('%s\t%s\t%s\n' % (ent[0],clust_num, ent[1]))
                    clust_num+=1


def main():
    conn = connectToDB()
    bininfo = dbrange(conn)
    abundance_list=abundance22merbin(conn,bininfo)
    results=outputwrite(abundance_list)
    

    
    
if __name__ == "__main__":
    main()











'''
Created on May 10, 2012

@author: atul
'''
