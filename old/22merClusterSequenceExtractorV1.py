# The purpose of this program is to extract the 22mer clusters from a database and find the absolute start and end locations of these clusters.
# The clustering script is just static, meaning that it only clusters sRNAs within a 500 nt window, while some clusters could actually be merged
# if one ends at the end of a window, and another starts at the window. So, this program has the ability to merge clusters and return all of the
# locations of these clusters.
#
#import os
#import string
import sys

import MySQLdb as mdb

# Simple function to connec to a database while prompting the user for a username and password.
# Return: The connection to the database
# ***NOTE*** This is hardcoded to only connect to the database 'reza' on taiji. This must be changed if to be used on other databases.
def connectToDB():
    # Variable which will hold the host server name
    hostname = 'taiji.dbi.udel.edu'
#    hostname = '127.0.0.1'
    
    # Variable which will hold the database name
    database = 'reza'
    
    # Conncection variable that will be used whenever we wish to query the database
    con = None

    # Prompt the user for a username, and then store that in the 'username' variable
#    username = raw_input('What is your username?\n')
    username='kakrana'

    # This block of code will prompt the user for a password.
#    print "What is your password?"
    print "Executing analysis...please be patient"
    
    # Disable echo so that the user's password will not be displayed on the screen to ensure privacy.
#    os.system("stty -echo")
#    password = raw_input()
    password='livetheday'
    
    # Reenable echo after the user has input the password so that further output can be seen.
#    os.system("stty echo")

    # Tryt o connec tto the database
    try:
        con = mdb.connect(hostname, username, password, database)

    # If we cannot connect to the database, send an error to the user and exit the program.
    except mdb.Error, e:
        print "Error %d: %s" % (e.args[0],e.args[1])
        sys.exit(1)

    return con

# Simple function that will use the connection to obtain all clusters that have at least 1 hit to a 22mer.
# param: connection to database
# return: 22mer clusters
def get22merClusters(con):
    # Set cur to be the cursor so we can execute a query
    cur = con.cursor()
    
    # Execute the query to select the clusters with at least 1 22mer hit.
    cur.execute("select cluster_ID, chr_id, position, all_lib_sum_all from clusterID_based_clustering2 where all_lib_sum_all>0")
    
    # Variable clusters will contain all 22mers. This is stored as a tuple containing the cluster_ID, chr_id, position, and all_lib_sum_all
    clusters = cur.fetchall()
#    print clusters
    return clusters


# Function that will find all 22mers that are located within each cluster. These 22mers will be appended
# to the set of all clusters so that we can easily access the 22mers associated with each cluster. 
# param: connection to the database
# param: list of clusters provided by the function get22merClusters
# return: new set of 22mers with the 22mer positions appended to the end
def gather22mersByClusters(con,clusters):
    # Variable newClusters must be created so that we can mutate the clusters list. The tuples are immutable once set, so a new variable
    # is necessary to store the 22mers the way we want.
    newClusters = []

    # This loop will iterate through all clusters. Inside, we will find all 22mers found within the range of each cluster.
    for cluster in clusters: # change to 'for cluster in clusters:' and remove next line when running for real
        

        # Set cur to be the cursor so we can execute a query
        cur = con.cursor()

        # This query will find the positions of each 22mer between the start location of the cluster, and the end position of the cluster, which
        # is just the start location+499. Because the variables are stored as long ints, we must convert them to strings to perform this query.
        cur.execute("select position from tag_position_summary where len=22 and chr_id= %s and (position between %s and %s);" % (str(cluster[1]), str(cluster[2]), str(cluster[2]+499)))

        # Variable rows will contain all 22mers that are found in the range of the cluster being examined.
        rows = cur.fetchall()

        # Create a temporary list called temp and initialize as an empty list
        temp = []

        # Because of the way mysqldb handles the queries, it is necessary that we create a new temp variable and store the 22mer positions
        # in this variable so that we can easily work with those numbers later. The numbers are stored as long ints, but we can just convert
        # them to ints here.
        for row in rows:
            temp.append(int(row[0]))

        # Append this cluster to the variable 'newClusters'. 'newClusters' contains all of the information that 'clusters' had, but
        # we convert all long ints to just ints, and we append temp(start positions of the 22mers) to the end.
        newClusters.append((cluster[0], int(cluster[1]), int(cluster[2]), int(cluster[3]), temp))

    return newClusters

# This function will merge two clusters, and return the merged clusters as one to the calling function
# param: two separate clusters to be merged
# return: The new merged clusters
def mergeClusters(cluster1, cluster2):
    # Initialize the new cluster as a list containing the empty string in the first element
    newCluster = ['']

    # Temporary variable to read in the characters of the clusterID from the first cluster
    readCharacter = cluster1[0][0]

    # Initialize our counter to 1, as we already read the first character of the clusterID previously
    i = 1

    # Iterate through the cluster1's clusterID until we find a '-'
    while readCharacter != '-':
        # Begin naming the newCluster with the initial part of the clusterID
        newCluster[0] = newCluster[0] + readCharacter
        readCharacter = cluster1[0][i]
        i = i+1

    # For the proper format of a clusterID, we must add a hyphen to the clusterID of the new cluster we are creating
    newCluster[0] = newCluster[0]+'-'

    # The window length of the new cluster will be 500 + the length of cluster1.
    windowLength = int(cluster1[0][i::])
    windowLength = windowLength+500
    newCluster[0] = newCluster[0]+str(windowLength)

    # The chromosome number and the starting position of the new cluster will be that of cluster1's value.
    newCluster.append(cluster1[1])
    newCluster.append(cluster1[2])

    # The number of hits will be the sum of the hits in cluster1 and cluster2
    newCluster.append(cluster1[3]+cluster2[3])

    # Add the positions of the 22mers, so add the entire 22mer list from cluster1 to newCluster, then iterate through the
    # list of 22mers from cluster2 and append them to newCluster's list
    newCluster.append(cluster1[4])
    for i in range(len(cluster2[4])):
        newCluster[4].append(cluster2[4][i])

    # ***Remove print function when running non-testing version***
#    print newCluster

    return newCluster

# This function will attempt to merge the clusters. We will just try to merge with the cluster that is immediately to the right.
# It is not necessary to try to merge to the left as that will be handled in the n-1 step. The way this function works is by
# trying to find if 22mers is 2 adjacent clusters are within 20 nt of one another. If this is so, then we will merge the 2 clusters.
# param: list of clusters with the 22mer positions appended to the end
# return: new list of merged clusters
def findClustersToMerge(clusters):
    i = 0
    for cluster in clusters:
        if i < len(clusters)-2: 
            if(clusters[i][4][len(clusters[i][4])-1]+21+20 > clusters[i+1][4][0]):
                clusters[i] = mergeClusters(clusters[i],clusters[i+1])
            
                # After merging clusters[i] and clusters[i+1], we no longer need clusters[i+1], so we will delete it from clusters
                del clusters[i+1]
                i = i-1

            i = i + 1
    return clusters

def trimClusters(clusters):
    finalclusters=[]
    for acluster in clusters:
#	print "this is a cluster"
#        print acluster
        start22mers= acluster[4][0]
        end22mers = acluster[4][-1]+21
        finalclusters.append((acluster[0:-1],start22mers,end22mers))        

    
    return finalclusters

#This function output the results to a file so it can be analyzed in anyway. I have mentioned all the entities
# like clust_name, chr_id, clust_start and more so that person reading the script could understand whats going in out_file
def WriteOutput(finalclusters):
	
	
	out_file = '22merClusterSequenceExtractor_cutoff1_results'
	fh_out = open(out_file, 'w')
	fh_out.write ('clust_name\tchr_id\tclust_start\tall_lib_sum\tstart22mer\tend22mer\n')
	
	print " The results for analysis are being written in file: %s" % (out_file)
	for acluster in finalclusters:
		clust_name=acluster[0][0]
#        print clust_name
		chr_id = acluster[0][1]
		clust_start = acluster[0][2]
		all_lib_sum = acluster[0][3]
		start22mer = acluster[1]
		end22mer = acluster[2]
		fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (clust_name,chr_id,clust_start,all_lib_sum,start22mer,end22mer))
		
		
	return out_file
	fh_out.close()


###Function to extract FASTA files for the cluster. Only the start22mer and end22mer position were used as they are result of cluster trimming 
###from the outfile made in 'Writeoutput' function

def GetClustSeq(con,filename):
	fh_in = open(filename,'r')
	##Waste the header
	fh_in.readline()
	##Open a file to write results in FASTA format
	out_fasta = '22merClusterSequenceExtractor_fasta'
	fh_out2 = open(out_fasta, 'w')

	for ent in fh_in:
#		print ent
		seq = ''
		ent = ent.strip('\n')
		ent = ent.split('\t')
		clust_name = ent[0]
		start22mer = int(ent[4])
		end22mer = int(ent[5])
		length_required = (end22mer-start22mer)+1
		chr_id = ent[1]
		print "The start position on cluster is %s and length extracted is:%s" % (start22mer,length_required)
		
		cur = con.cursor()
		cur.execute("SELECT SUBSTRING(chromosome, %s, %s) FROM SLY_ITAG2_genome.chrom_sequence where chr_id = %s and strand = 'w'", (start22mer, length_required, chr_id))
		string=cur.fetchall()
		seq =string[0][0]
#		print seq
		fh_out2.write('>%s\n%s\n' % (clust_name,seq))
#    	break
    	
    	
	return out_fasta
	fh_in.close()
	fh_out2.close()
    			
		

    

def main():
    con = connectToDB()
    clusters = get22merClusters(con)
    clusters = gather22mersByClusters(con, clusters)
    clusters = findClustersToMerge(clusters)
    finalclusters = trimClusters(clusters)
    out_file= WriteOutput(finalclusters)
#    GetClustSeq(con,out_file)


if __name__ == "__main__":
    main()

