##The function of this module is to calculate the concentration of 22mers and 21mers in the cluster

import argparse
import mysql.connector as sql

##Variables
clust_file = '50_perc_top_EPRV.list'

##- - - - - - - -- - - - - - Arguments - - - - - - - - - - - - - -##
parser = argparse.ArgumentParser()

#-db DATABSE -u USERNAME -p PASSWORD -size 20
parser.add_argument("-host", "--hostname", dest = "hostname", default = "taiji.dbi.udel.edu", help="Server name")
parser.add_argument("-db", "--database", dest = "db", help="Database name")##The db name with tag position summary
parser.add_argument("-u", "--username",dest ="username",default = "kakrana", help="User name")
parser.add_argument("-p", "--password",dest = "password", default = "livetheday", help="Password")

args = parser.parse_args()

print( "Hostname {} db {} User {}".format(
        args.hostname,
        args.db,
        args.username,
        args.password,
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

def ReadClust(filename):
    clust_list = []
    fh_in = open(clust_file, 'r')
    fh_out = open('miR_conc', 'w')
    
    for clust in fh_in:
        clust = clust.strip('\n')
        clust = clust.split(':')
    #    print(clust)
        chr = clust[1]
        loci = clust[2].split('-')
        start = loci[0]
        range = loci[1]
#        print(chr, start,range)
        clust_list.append((chr,start,range))
    return clust_list
    
def CalcmiRConc(con, clust,mer): ###'clust' comes from clust_list in format ((chr,start,range))
#    clust_abun = [] ##list to hold abundance values of 21mers and 22mers
#    for len in range (20,22):
    #    print(clust[0:])
    #    print(int(clust[1])+int(clust[2]))
        cur = con.cursor()
        cur.execute("select tag,norm_sum,hits from %s.tag_position_summary where len = %s and chr_id= %s and (position between %s and %s);" % (args.db, mer, clust[0], clust[1], int(clust[1])+int(clust[2])))
        clust_abun = cur.fetchall()
        return clust_abun
    #    print (mircount[0])



##MAIN
con = ConnectToDB()
clust_list = ReadClust(clust_file)
#print(clust_list)

##Writing results

fh_out = open('clust_abundances' , 'w')
fh_out.write('Cluster    \tLen\tAbundance\tTags\tLen\tAbundance\tTags\n')


clust_num = 0 ## To keep track of number of clusters analyzed in analysis
for clust in clust_list:
    
    abun_list = [] ## this resets after every cluster and fills values from next cluster
    clust_name = 'clust_'+clust[0]+':'+clust[1]+'-'+clust[2]### A name for cluster
    print ('%s is being analyzed' % (clust_name))
    
    for mer in range (21,23):## For different size of miRNAs                 
        mer_abun = CalcmiRConc(con, clust, mer)
        num_tags = len(mer_abun) ### This will be used to calculate average 22mer or 21mer abundance in cluster
        abun_sum = 0## Counter to add the number of 22mer or 21mer abundance
        for ent in mer_abun: ## Loop to add up abundances for all the tags and later divide by total tags (num_tags) for average abundance
            abun_sum += ent[1]
        abun_list.append((mer,abun_sum,num_tags))
    
    #print(abun_list)
    fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (clust_name,abun_list[0][0],abun_list[0][1],abun_list[0][2],abun_list[1][0],abun_list[1][1],abun_list[1][2]))
    clust_num += 1
    
print ('Analysis has completed successfully\n')
print ('Total number of clusters anlyzed in analysis: %s' % (clust_num))
    






'''
Created on Jul 17, 2012

@author: atul
'''
