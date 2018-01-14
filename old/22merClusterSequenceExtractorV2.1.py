# The purpose of this program is to extract the 22mer clusters from a database and find the absolute start and end locations of these clusters.
# The clustering script is just static, meaning that it only clusters sRNAs within a 500 nt window, while some clusters could actually be merged
# if one ends at the end of a window, and another starts at the window. So, this program has the ability to merge clusters and return all of the
# locations of these clusters.
#
import os
import string
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
	database = 'AT_sbs_sRNA_v4'
	
	# Conncection variable that will be used whenever we wish to query the database
	con = None

	# Prompt the user for a username, and then store that in the 'username' variable
##        username = raw_input('What is your username?\n')
	username = 'kakrana'
	# This block of code will prompt the user for a password.
#        print "What is your password?"
	
	# Disable echo so that the user's password will not be displayed on the screen to ensure privacy.
#        os.system("stty -echo")
#        password = raw_input()
	
	# Reenable echo after the user has input the password so that further output can be seen.
#        os.system("stty echo")
	password = 'livetheday'
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
	cur.execute("select cluster_ID, chr_id, position, wt_mock_lib_sum_all from clusterID_based_clustering where (wt_mock_lib_sum_all >=20) and (Col10a_sum_all > 0 or Col10b_sum_all > 0 or Col10c_sum_all > 0 or Col7a_sum_all > 0 or Col7b_sum_all > 0 or Col7c_sum_all);")
	
	# Variable clusters will contain all 22mers. This is stored as a tuple containing the cluster_ID, chr_id, position, and pot_all_lib_sum_all
	clusters = cur.fetchall()

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
	#for cluster in clusters: # change to 'for cluster in clusters:' and remove next line when running for real
	for cluster in clusters:

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

	return newCluster

# This function will attempt to merge the clusters. We will just try to merge with the cluster that is immediately to the right.
# It is not necessary to try to merge to the left as that will be handled in the n-1 step. The way this function works is by
# trying to find if 22mers is 2 adjacent clusters are within 10 nt of one another. If this is so, then we will merge the 2 clusters.
# param: list of clusters with the 22mer positions appended to the end
# return: new list of merged clusters
def findClustersToMerge(clusters):
	i = 0
	for cluster in clusters:
		if i < len(clusters)-2: 
			if(clusters[i][4][len(clusters[i][4])-1]+21+10 > clusters[i+1][4][0]):
				clusters[i] = mergeClusters(clusters[i],clusters[i+1])
			
				# After merging clusters[i] and clusters[i+1], we no longer need clusters[i+1], so we will delete it from clusters
				del clusters[i+1]
				i = i-1

			i = i + 1
	return clusters

def trimClusters(clusters):
	finalclusters=[]

	for cluster in clusters:
		start22mers= cluster[4][0]
		end22mers = cluster[4][-1]+21
		finalclusters.append((cluster[0], cluster[1], cluster[2], cluster[3], cluster[4], start22mers, end22mers))
	
	return finalclusters

#This function output the results to a file so it can be analyzed in anyway. I have mentioned all the entities
# like clust_name, chr_id, clust_start and more so that person reading the script could understand whats going in out_file
def WriteOutput(finalclusters):	
	out_file = '22merClusterSequenceExtractor_results'
	fh_out = open('22merClusterSequenceExtractor_results', 'w')
	fh_out.write ('cluster_ID\tchr_id\tclust_start\tall_lib_sum\tstart22mer\tend22mer\n')
	
#	print " The results for analysis are being written in file: %s" % (out_file)
	for cluster in finalclusters:
		fh_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (cluster[0],cluster[1],cluster[2],cluster[3],cluster[5],cluster[6]))
		
	fh_out.close()
	return out_file

###Function to extract FASTA files for the cluster. Only the start22mer and end22mer position were used as they are result of cluster trimming 
###from the outfile made in 'Writeoutput' function
def GetClustSeq(con,clusters,filename):
	fh_out = open('22merClusterSequenceExtractor_fasta','w')
	for cluster in clusters:
		cur = con.cursor()
		cur.execute("select substring(chromosome, %s, %s) from AT_TAIR10_genome.chrom_sequence where chr_id = %s and strand = 'w'", (str(cluster[5]), str(cluster[6]-cluster[5]), str(cluster[1])))
		temp = cur.fetchall()
		fh_out.write('>%s\n%s\n' % (cluster[0],temp[0][0]))

	fh_out.close()

def main():
	con = connectToDB()
	print "Gathering clusters"
	clusters = get22merClusters(con)
	print "Adding 22mer locations to the clusters"
	clusters = gather22mersByClusters(con, clusters)
	print "Merging appropriate cluster windows"
	clusters = findClustersToMerge(clusters)
	print "Trimming clusters"
	finalclusters = trimClusters(clusters)
	print "Writing trimmed clusters to file"
	out_file= WriteOutput(finalclusters)
	print "Getting cluster sequences"
	GetClustSeq(con,finalclusters,out_file)

if __name__ == "__main__":
	main()

