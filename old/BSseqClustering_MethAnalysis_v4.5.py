##V3-Improvements
#1. Calculates cluster and records on the basis of libraries and not on Chromosomes
#2. Also records clusters with 0 value as multiple samples should show equal number of clusters for each sample

##v4- Improvements
#1. User inputs the library_ids, if not than all the libraries are used
#2. Variables like db, user, password and binsize are now given from command line
##4.5 Inprovements
#1 COntext based analysis and output file generation
#2. Mail sent on completion

##Script written by Atul:atulkakrana@gmail.com

import mysql.connector as sql
import sys
import re
import argparse
import smtplib
import datetime



##- - - - - - - - - - - - - GloBal Variables/Config- - - - - - - - - - - ##

#lib_ids = [1495,1496]

##- - - - - - - -- - - - - - Arguments - - - - - - - - - - - - - -##
parser = argparse.ArgumentParser()

#-db DATABASE -u USERNAME -p PASSWORD -size 20000
parser.add_argument("-host", "--hostname", dest = "hostname", default = "tarkan.dbi.udel.edu", help="Server name")
parser.add_argument("-t", "--table", dest = "table", default = "C_table", help="table to be analyzed")
parser.add_argument("-db", "--database", dest = "db", default = "MG_priv_BSseq", help="Database name")
parser.add_argument("-lib", "--libraries", dest = "lib", help="libraries to be analyzed")
parser.add_argument("-u", "--username",dest ="username", help="User name")
parser.add_argument("-p", "--password",dest = "password", default = "123456", help="Password")###Hard code pasword in place of 12345
parser.add_argument("-chr", "--chromosomes",dest = "chr", help="Choromosomes")
parser.add_argument("-size", "--binsize",dest = "binsize", default = 500, help="Size", type=int)
parser.add_argument("-mail", "--sendmail",dest = "mail", default = "N", help="Y- Sends email to specified account, N- Does not sends e-mail", type=str)
parser.add_argument("-to", "--mailto",dest = "to", default = "lipingchuan@gmail.com", help="Sends email to specified account", type=str)####Hard Code e-mail ID here

args = parser.parse_args()

print( "hostname {} table {} db {} lib {} user {} chr {} size {} ".format(
        args.hostname,
        args.table,
        args.db,
        args.lib,
        args.username,
        #args.password,
        args.chr,
        args.binsize
        ))


##- - - - - - - - - - - - - Modules - - - - - - - - - - - - - - - ##

#########################################################################################################################################################

def sendmail():
        print('Sending job complete mail')
        to = args.to
        gmail_user = 'daemon.blake.lab@gmail.com'
        gmail_pwd = 'WeLComE@MeYersLab!'
        smtpserver = smtplib.SMTP("smtp.gmail.com",587)
        smtpserver.ehlo()
        smtpserver.starttls()
        smtpserver.ehlo
        smtpserver.login(gmail_user, gmail_pwd)
        header = 'To:' + to + '\n' + 'From: ' + gmail_user + '\n' + 'Subject:Script run finished \n'
        #print (header)
        msg = (header + 'Master, \n\nYour BS seq clustering and methylation analysis script on %s using %s has just finished.\nStart Time:%s\nEnd Time:%s\n\nThanks.\n\nLab mail daemon\nServer room\nDBI' % (args.hostname,args.db,start_time,end_time))
        #msg = (header + '\n This is test msg from your server \n\n')
        smtpserver.sendmail(gmail_user, to, msg)
        print ('Mail Sent!')
        smtpserver.close()      
    
#############################################################################################################################################################

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

#############################################################################################################################################################

def DBRange(con, chromosomes):
    bininfo = []
    
    if chromosomes is None: 
        print('No chromosomes were specified from command line, using all the chromosomes as default\n')
        cur = con.cursor()
        cur.execute("SELECT distinct(chr_id) FROM %s.%s;" %(args.db, args.table))##Calculates the maximum number of chromosomes i.e total number of chromosomes
        maxchr= cur.fetchall()
        chromosomes = [] ## List to keep chromomosome ids
        for i in maxchr:
                chromosomes.append(i[0])
        #print(chromosomes)
                
        
    else:
        if args.chr.count(',') > 0: ## What if it is just one and no split required, so we check the number of inputs i.e if there is a comma do split
            chromosomes = args.chr.split(',')##List is automatically made by split
        else:
            chromosomes = []##Chromosmes has to be in form of list ven if single chromosome is input from commandline
            chromosomes.append(args.chr)
        
    for chromo in chromosomes:##For all the chromosomes
        cur = con.cursor()
        cur.execute("SELECT max(position) FROM %s.%s WHERE chr_id = %s;" % (args.db,args.table,chromo))##Get the length of each chromosome
        maxlen  = cur.fetchall()

  
        bins=int(maxlen[0][0]//args.binsize)##List index were used as the maxlen is a tuple
        lastbin=int(maxlen[0][0]%args.binsize)##List index were used as the maxlen is a tuple
        
        bininfo.append((chromo,bins,lastbin))##chr num, total number of clusters, last cluster remaining length

        
#    print(bininfo)

    return bininfo,chromosomes##List holdings number of bins and size of last/remainder bin

############################################################################################################################################################

def GetCTable(con, bininfo, libs):
    
    ##Check for library input, if empty use all libraries as default
    if libs is None:
        print('No libraries were specified, using all the libraries as default\n')    
        #Get the lib_ids
        cur=con.cursor()
        cur.execute("SELECT distinct(lib_id) FROM %s.%s;" % (args.db,args.table))
        libs=cur.fetchall()
        lib_ids = []
        for i in libs:
                lib_ids.append(i[0])    
        
    else:
        if args.lib.count(',')  > 0:
            lib_ids = args.lib.split(',')###Make list of lib_ids from the command line output  
        else:
            lib_ids = []## Even if the commandline input is a single library i.e 1495 it has to be in form of list for the next 'for' loop
            lib_ids.append(args.lib)##otherwise the single lib_id i.e 1495 will be split into 1,4,9,5 in the next 'for'loop
    #print (lib_ids)
    
    ##MASTER LISTS that will hold information for every bin-these are poulated from sub lists that are pupulated at entry in bin level
    C_list=[]##Big list holding information for every cluster as a tuple: lib_id, chrNo, cluster_start, cluster_end, value
    CHH_list = []##Big list holding CHH type information for every cluster as a tuple: lib_id, chrNo, cluster_start, cluster_end, value
    CG_list =[]##Big list holding CG type information for every cluster as a tuple: lib_id, chrNo, cluster_start, cluster_end, value
    CHG_list =[]##Big list holding CHG type information for every cluster as a tuple: lib_id, chrNo, cluster_start, cluster_end, value
    libraries =[]
    
    for lib_id in lib_ids:
        #lib_id = lib[0] ##Uncomment incase you are automatically getting lib_ids
        libraries.append(lib_id)
        print('\nThe library of id %s is being analyzed' %(lib_id))
        
        #print (bininfo[0])
        #break
        for ent in bininfo:## i.e number of chromosomes as each chromosome has a tuple in bininfo, revert to: for ent in bininfo :after testing
            print ('Caching the values of "cmethylated" and "totalc" in Chromosome %s' % (ent[0]))     
    
            pos=1##Variable to hold position on chromosome through the loop, set to 1 when next chromosome is analyzed
            cluster = 1
        
            ##Anlaysis at bin level
            for i in range(ent[1]):##ent[1] is number of bins i.e for every bin## use : for i in range(ent[1]) : after testing is complete
                cur = con.cursor()
                ##cur.execute("select sum(cmethylated),sum(totalc),count(*) from %s.%s where chr_id = %s and lib_id='%s' and position between %s and %s;" % (args.db,args.table,ent[0],lib_id,pos,pos+(args.binsize-1)))
                cur.execute("select cmethylated,totalc,type from %s.%s where chr_id = %s and lib_id='%s' and position between %s and %s;" % (args.db,args.table,ent[0],lib_id,pos,pos+(args.binsize-1)))
                bindata = cur.fetchall()
                ##print(bindata[0])
                #print(bindata)
                #for z in bindata:
                #        print(z)
                
                ## Temporary storage list, washed out after evry bin
                c_bin_list=[0]*3##Will keep information of every entry of any type in format: c_methylated, totalc, count
                chh_bin_list=[0]*3##Will keep information of CHH type in format: c_methylated, totalc, count
                cg_bin_list=[0]*3##Will keep information of CG  type in format: c_methylated, totalc, count
                chg_bin_list=[0]*3##Will keep information of CHG type in format: c_methylated, totalc, count
                #print(c_bin_list)
                
                ##Filter the data on the basis of type
                ##Example bindata
                ##(15, 81, 'CHH')
                ##(8, 77, 'CHH')
                ##(9, 68, 'CHH')
                ##(11, 259, 'CHH')
                
                ###Analysis at every methylation entry in a bin level- lowest level
                for met in bindata:##for every methylation entry in this bin, filter the methylation type and sum methylation, totalc
                        #print(met)
                        ##Overall methylation in bin including all types
                        c_bin_list[0]+=met[0]##C_methylated
                        c_bin_list[1]+=met[1]
                        c_bin_list[2]+=1##Counting total meth entries in bin
                        
                        ##Filter CHH
                        if met[2] == 'CHH':
                                #print('This is CHH:',met[2])##just to make sure that right type is being read
                                chh_bin_list[0]+=met[0]##C_methylated
                                chh_bin_list[1]+=met[1]##Totalc
                                chh_bin_list[2]+=1##Counting total CHH entries in bin
                        ##Filter CG
                        elif met[2] == 'CG':
                                cg_bin_list[0]+=met[0]##C_methylated
                                cg_bin_list[1]+=met[1]##Totalc
                                cg_bin_list[2]+=1##Counting total CHH entries in bin
                        ##Filter CHG
                        elif met[2] == 'CHG':
                                chg_bin_list[0]+=met[0]##C_methylated
                                chg_bin_list[1]+=met[1]##Totalc
                                chg_bin_list[2]+=1##Counting total CHH entries in bin
                #print(c_bin_list)
                #print(chh_bin_list)
                #print(cg_bin_list)
                #print(chg_bin_list)
                
                ###Add these lists to master lists: 
                
                ##C_list
                if c_bin_list[2] == 0:##If there is no entry found in bin i.e count == 0
                    C_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),0,'ALL' ))##Chr_num+bin number as unique identifier       
                else:
                    C_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),float( (int(bindata[0][0])*100)/int(bindata[0][1]) ),'ALL' ))##Chr_num+bin number as unique identifier, c_methylated*100/sum of total_C

                ##CHH_list
                if chh_bin_list[2] == 0:
                        CHH_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),0,'CHH' ))##Chr_num+bin number as unique identifier
                else:
                        CHH_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),float( (int(bindata[0][0])*100)/int(bindata[0][1]) ),'CHH' ))##Chr_num+bin number as unique identifier, c_methylated*100/sum of total_C   
                        
                ##CG_list
                if cg_bin_list[2] == 0:
                        CG_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),0,'CG' ))##Chr_num+bin number as unique identifier
                else:
                        CG_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),float( (int(bindata[0][0])*100)/int(bindata[0][1]) ),'CG' ))
                
                ##CHG_list       
                if chg_bin_list[2] == 0:
                        CHG_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),0,'CHG' ))##Chr_num+bin number as unique identifier
                else:
                        CHG_list.append((lib_id,'chr'+str(ent[0]),pos,pos+(args.binsize-1),float( (int(bindata[0][0])*100)/int(bindata[0][1]) ),'CHG' ))
                        
                #pass
                
                pos+=args.binsize##Change the position frame to next bin
                cluster+=1##change the cluster number 
            
            ##Final cluster of irregular length    
            finalbin_num = ent[1]+1
            finalbin_start = (ent[1]*500)+1##After the clusters of length 500 are analyzed,last cluster start position is calculated
            finalbin_end= (finalbin_start+ent[2])-1
    #        print(finalbin_start, finalbin_end)
            
            ##Final Bin query
            cur2 = con.cursor()
            cur2.execute("select sum(cmethylated),sum(totalc),count(*) from %s.%s where chr_id = %s and lib_id='%s' and position between %s and %s;" % (args.db,args.table, ent[0],lib_id,finalbin_start, finalbin_end))
            bindata_last = cur2.fetchall()
            
            ###All the list intialized for final cluster as well
            ## Temporary storage list, washed out after evry bin
            c_bin_list=[0]*3##Will keep information of every entry of any type in format: c_methylated, totalc, count
            chh_bin_list=[0]*3##Will keep information of CHH type in format: c_methylated, totalc, count
            cg_bin_list=[0]*3##Will keep information of CG  type in format: c_methylated, totalc, count
            chg_bin_list=[0]*3##Will keep information of CHG type in format: c_methylated, totalc, count
            #print(c_bin_list)
           
            ###Analysis at every methylation entry in a bin level- lowest level
            for met in bindata_last:##for every methylation entry in this bin, filter the methylation type and sum methylation, totalc
                    #print(met)
                    ##Overall methylation in bin including all types
                    c_bin_list[0]+=met[0]##C_methylated
                    c_bin_list[1]+=met[1]
                    c_bin_list[2]+=1##Counting total meth entries in bin
                    
                    ##Filter CHH
                    if met[2] == 'CHH':
                            #print('This is CHH:',met[2])##just to make sure that right type is being read
                            chh_bin_list[0]+=met[0]##C_methylated
                            chh_bin_list[1]+=met[1]##Totalc
                            chh_bin_list[2]+=1##Counting total CHH entries in bin
                    ##Filter CG
                    elif met[2] == 'CG':
                            cg_bin_list[0]+=met[0]##C_methylated
                            cg_bin_list[1]+=met[1]##Totalc
                            cg_bin_list[2]+=1##Counting total CHH entries in bin
                    ##Filter CHG
                    elif met[2] == 'CHG':
                            chg_bin_list[0]+=met[0]##C_methylated
                            chg_bin_list[1]+=met[1]##Totalc
                            chg_bin_list[2]+=1##Counting total CHH entries in bin
             ###Add these lists to master lists one last time: 
                
                ##C_list
            if c_bin_list[2] == 0:##If there is no entry found in bin i.e count == 0
                C_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,0,'ALL' ))##Chr_num+bin number as unique identifier       
            else:
                C_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,float( (int(bindata_last[0][0])*100)/int(bindata_last[0][1]) ),'ALL' ))##Chr_num+bin number as unique identifier, c_methylated*100/sum of total_C

            ##CHH_list
            if chh_bin_list[2] == 0:
                    CHH_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,0,'CHH' ))##Chr_num+bin number as unique identifier
            else:
                    CHH_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,float( (int(bindata_last[0][0])*100)/int(bindata_last[0][1]) ),'CHH' ))##Chr_num+bin number as unique identifier, c_methylated*100/sum of total_C   
                    
            ##CG_list
            if cg_bin_list[2] == 0:
                    CG_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,0,'CG' ))##Chr_num+bin number as unique identifier
            else:
                    CG_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,float( (int(bindata_last[0][0])*100)/int(bindata_last[0][1]) ),'CG' ))
            
            ##CHG_list       
            if chg_bin_list[2] == 0:
                    CHG_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,0,'CHG' ))##Chr_num+bin number as unique identifier
            else:
                    CHG_list.append((lib_id,'chr'+str(ent[0]),finalbin_start,finalbin_end,float( (int(bindata_last[0][0])*100)/int(bindata_last[0][1]) ),'CHG' ))    
    
    list_ring = (C_list,CHH_list,CG_list,CHG_list)##A key ring of all lists
    return list_ring,libraries###All master lists with library information, list frmat: lib_id, chrNo, cluster_start, cluster_end, value,Type

##################################################################################################################################################################

def outputwrite(alist,libraries,chromosomes):##alist format lib_id,chr+chrnumber,posstart,posend,value,type
    ##Detect context type i.e ALL,CHH,CG,CHG to be inculded in filename
    met_type = str(alist[0][5])
    chromo =''##Concatanate chromosomes name to be added in file name
    for i in chromosomes:##Lopp through chrmosomes in chromosomes list and append them to chromo string
        chromo += str('_')+str(i)
        
    #print (met_type)
    for lib in libraries:
        print ('\nThe results are being sorted and written for library %s' % (lib))

        ##OPen the file here as it made individually for every librray and every chromosome
        fh_out=open('lib%s_chr%s_%s_bin%s_C_data' % (lib,chromo,met_type,args.binsize), 'w')##Open one file for each chromosome
#               fh_out.write('Chromosome\tClust Start\tClust End\tConverted_value\n')
        #clust_num=1
        for ent in alist:## Go through each entry of context type list, format: lib_id, chrNo, cluster_start, cluster_end, value,Type
            if ent[0] == lib:##If entry has Chromosome number than write to file
                fh_out.write('%s %s %s %s\n' % (ent[1],ent[2],ent[3],round(ent[4],4)))
                #clust_num+=1
                                

##- - - - - - - - - - - - - - - - MAIN - - - - - - - - - - - - - -##

def main():
        ##Get start time
        start_time = datetime.datetime.now()
        ##Make a connection
        con = ConnectToDB()
        ##Make bins
        bininfo,chromosomes=DBRange(con, args.chr)
        ##Critical....do analysis
        list_ring,libraries=GetCTable(con,bininfo,args.lib)
        for alist in list_ring:
                ##Write results
                outputwrite(alist,libraries,chromosomes)
        
        ##End Time
        end_time = datetime.datetime.now()
        
        ##Print time info
        print('\nTime analysis started:%s\nTime analysis ended:%s' % (start_time,end_time))
        
        ##Send mail to lazy bones
        
        if args.mail == 'Y':  
            sendmail()
        else:
            print('\n~~~~~~~~~~~Script run finished~~~~~~~~~~~~~~\n')
        

if __name__ == "__main__":
	main()




'''
Created on Jun 12, 2012

@author: atul
'''
