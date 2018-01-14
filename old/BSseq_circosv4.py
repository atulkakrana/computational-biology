##V3-Improvements
#1. Calculates cluster and records on the basis of libraries and not on Chromosomes
#2. Also records clusters with 0 value as multiple samples should show equal number of clusters for each sample

##v4- Improvements
#1. User inputs the library_ids, if not than all the libraries are used
#2. Variables like db, user, password and binsize are now given from command line

##Script written by Atul:atulkakrana@gmail.com

import mysql.connector as sql
import sys
import re
import argparse



##- - - - - - - - - - - - - GloBal Variables/Config- - - - - - - - - - - ##

lib_ids = [1495,1496]

##- - - - - - - -- - - - - - Arguments - - - - - - - - - - - - - -##
parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20000
parser.add_argument("-host", "--hostname", dest = "hostname", default = "tarkan.dbi.udel.edu", help="Server name")
parser.add_argument("-db", "--database", dest = "db", default = "MG_priv_BSseq", help="Database name")
parser.add_argument("-u", "--username",dest ="username", help="User name")
parser.add_argument("-p", "--password",dest = "password", help="Password")
parser.add_argument("-size", "--binsize",dest = "binsize", help="Size", type=int)

args = parser.parse_args()

print( "Hostname {} db {} User {} Password {} size {} ".format(
        args.hostname,
        args.db,
        args.username,
        args.password,
        args.binsize
        ))


##- - - - - - - - - - - - - Modules - - - - - - - - - - - - - - - ##

def ConnectToDB():
    print ('Trying to connect to mySQL server')
    # Try to connect to the database
    try:
        con=sql.connect(host=args.hostname, user= args.username, passwd= args.password)
        print ('\nConnected to Database\n')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con

def DBRange(con):
    bininfo = []
    cur = con.cursor()
    cur.execute("SELECT MAX(chr_id) FROM %s.C_table;" %(args.db))##Calculates the maximum number of chromosomes i.e total number of chromosomes
    maxchr= cur.fetchall()
    
    for chromo in range (maxchr[0][0]):##For all the chromosomes
        cur = con.cursor()
        cur.execute("SELECT max(position) FROM %s.C_table WHERE chr_id = %s;" % (args.db,chromo+1))##Get the length of each chromosome | chromo+1 because first loop as entry 0 so add one to make it chrmosome 1
        maxlen  = cur.fetchall()

      
        bins=int(maxlen[0][0]//args.binsize)##List index were used as the maxlen is a tuple
        lastbin=int(maxlen[0][0]%args.binsize)##List index were used as the maxlen is a tuple
            
        bininfo.append((chromo+1,bins,lastbin))##chr num, total number of clusters, last cluster remaining length
        
#    print(bininfo)
    return bininfo##List holdings number of bins and size of last/remainder bin

def GetCTable(con, bininfo, lib_ids):
    ##Check for library input, if empty use all libraries as default
    if not lib_ids:
        print('No libraries were specified, using all the libraries as default')    
        #Get the lib_ids
        cur=con.cursor()
        cur.execute("SELECT distinct(lib_id) FROM %s.C_table;" % (args.db))
        lib_ids=cur.fetchall()
    
    C_list=[]##Big list holding information for every cluster as a tuple: lib_id, chrNo, cluster_start, cluster_end, value
    libraries =[]
    for lib_id in lib_ids:
#        lib_id = lib[0] ##Uncomment incase you are automatically getting lib_ids
        libraries.append(lib_id)
        print('The library of id %s is being analyzed' %(lib_id))
        
        for ent in bininfo:## i.e number of chromosomes as each chromosome has a tuple in bininfo 
            print ('\nCaching the values of "cmethylated" and "totalc" in Chromosome %s\n' % (ent[0]))     
    
            pos=1##Variable to hold position on chromosome through the loop, set to 1 when next chromosome is analyzed
            cluster = 1
        
            for i in range(ent[1]):##ent[1] is number of bins i.e for every bin## use : for i in range(ent[1]) : after testing is complete
                cur = con.cursor()
                cur.execute("select sum(cmethylated),sum(totalc),count(*) from %s.C_table where chr_id = %s and lib_id=%s and position between %s and %s;" % (args.db,ent[0],lib_id,pos,pos+(args.binsize-1)))
                bindata = cur.fetchall()
                #print(bindata[0])
                if bindata[0][2] == 0:
                    #print(lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),0 )
                    C_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),0 ))##Chr_num+bin number as unique identifier
                    pos+=args.binsize
                    cluster+=1
    #            if not bindata:
    #            if bindata is None:
    #            if bindata.rowcount==0:
                    #print ('NO rows found')
                    
    
                else:
                    #print(lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),float( (int(bindata[0][0])*100)/int(bindata[0][1]) ))
                    C_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),float( (int(bindata[0][0])*100)/int(bindata[0][1]) ) ))##Chr_num+bin number as unique identifier
                    pos+=args.binsize
                    cluster+=1
    #                pass
            
            ##Final cluster of irregular length    
            finalbin_num = ent[1]+1
            finalbin_start = (ent[1]*500)+1##After the clusters of length 500 are analyzed,last cluster start position is calculated
            finalbin_end= (finalbin_start+ent[2])-1
    #        print(finalbin_start, finalbin_end)
            cur2 = con.cursor()
            cur2.execute("select sum(cmethylated),sum(totalc),count(*) from %s.C_table where chr_id = %s and lib_id=%s and position between %s and %s;" % (args.db, ent[0],lib_id,finalbin_start, finalbin_end))
            bindata_last = cur2.fetchall()
            #print(bindata_last)
            if bindata_last[0][2] is 0:
                #print (lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,0)
                C_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,0))##Chr_num_bin number as unique identifier
                #print('No rows found')
                
            else:
                #print (lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,float( (int(bindata_last[0][0])*100)/int(bindata_last[0][1]) ) )
                C_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,float( (int(bindata_last[0][0])*100)/int(bindata_last[0][1]) ) ))##Chr_num_bin number as unique identifier
            
    return C_list,libraries

def outputwrite(C_list,libraries):
#    chr_num_re=re.compile('\d{1,2}')
#    max_chr_num=re.findall(chr_num_re,C_list[-1][0])##Last chromosome number extracted from abundance list
    for lib in libraries:
        
        
#    print C_list[-1]
#    print max_chr_num[0]###Index used to get integer
#    for i in range (1,int(max_chr_num[0])+1):##For number of chromosomes
        print ('The results are being sorted and written for library %s\n' % (lib))
        fh_out=open('C_data_in_bin_lib%s' % (lib), 'w')##Open one file for each chromosome
#        fh_out.write('Chromosome\tClust Start\tClust End\tConverted_value\n')
        clust_num=1
        for ent in C_list:## Go through each entry of abundance list
#            print ent[0]        
            if ent[0] == lib:##If entry has Chromosome number than write to file
#                print 'yes'
                fh_out.write('%s %s %s %s\n' % (ent[1],ent[2],ent[3],round(ent[4],4)))
                clust_num+=1

##- - - - - - - - - - - - - - - - MAIN - - - - - - - - - - - - - -##

con = ConnectToDB()
bininfo=DBRange(con)
C_list,libraries=GetCTable(con,bininfo,lib_ids)
outputwrite(C_list,libraries)

##- - - - - - - - - - - - - - - - END - - - - - - - - - - - - - - ##




'''
Created on Jun 12, 2012

@author: atul
'''




'''
Created on Jul 22, 2012

@author: atul
'''
