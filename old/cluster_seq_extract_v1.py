import mysql.connector as sql

##I/O files
#fh_out = open('cluster_sequences', 'w')


#conn=sql.connect(host='taiji.udel.edu', user='kakrana', passwd='******')###Use to work directly from server
conn=sql.connect(host='127.0.0.1', user='kakrana', passwd='******')##Connecting via tunnel, start ssh tunnel first for required server
cur= conn.cursor() 

##query1 to find total number of entries in cluster table
cur.execute('select count(*) from reza.clusterID_based_clustering')##Get total number of clusters
number_ent = cur.fetchall()
entries = number_ent[0][0]
print('Total number of entries in Cluster table are: %s\n' % (entries))

##query 2 why not take all cluster IDs in a list - FOR TEST FIRST 1000 ENTRIES ARE USED, 'limit 1000' should be removed
cur.execute('select cluster_id from reza.clusterID_based_clustering limit 1000')##Get total number of clusters
cluster_ids = cur.fetchall()
#alist = cluster_ids.list()


## LOOP: Every single cluster_id is broken to take name, strat/stop position on chromosome
for i in cluster_ids:
#    print (i)
    name = i[0]
    ent = i[0].split(':')
    chromo = ent[1]
    start_pos = int(ent[2].split('-')[0])
    end_pos = int(start_pos+500)
    print('The name is %s, chromosome is %s and positions are %s-%s' % (name, chromo, start_pos,end_pos))
#    break
## Query 3: Extract the sequences from Genome DB
    cur.execute("SELECT SUBSTRING(chromosome, %s, %s) FROM SLY_ITAG2_genome.chrom_sequence where chr_id = %s and strand = 'w'", (start_pos, end_pos, chromo))
    string=cur.fetchall()
    seq =string[0][0]
#    fh_out.write('>%s\n%s\n'% (name,seq))## Writing the sequence to file in fasta format, name is same as cluster name
      
    print(seq)
    break

print(number_ent[0][0])

#fh_out.close()

'''
Created on Apr 24, 2012

@author: setu
'''
