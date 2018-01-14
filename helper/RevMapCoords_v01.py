#!/usr/local/bin/python3

##Written for Revmapping of target findr results for Tow
##Maintainer kakrana@udel.edu


import os
import sys

####Settings######

gene_coords = 'gene_coords'
res_file = 'Comb_targets_Genes_parse_v3.csv'


def RevMapCoord(gene_coords,res_file):
    
    ###Make dictionary of reverse coords file  
    global coord_dict_wat, coord_dict_crick
    coord_dict_wat = {} ## Dictionary of genes at watson strand, +
    coord_dict_crick = {}###Dictionary of genes at crick strand, - strand
    
    fh_coords = open(gene_coords, 'r')
    
    for line in fh_coords:
        #print (line)
        i = line.strip('\n').split(',')
        #print(i)
        strand = i[1]
        if strand == 'c':### if entry in reverse strand
            atuple = (i[0],i[4],i[5])
            coord_dict_crick[i[2]] = atuple###Gene name as key and chr_id,strand end and gene type as value
        else:
            atuple = (i[0],i[3],i[5])
            coord_dict_wat[i[2]] = atuple##Gene name as key and chr_id,strand start and gene type as value
    print ('**Strand dictionary made**')
    
    ###List of results
    fh_in = open(res_file, 'r')### Results to be reverse mapped
    fh_in.readline()##Remove header
    res_list= [] ##list of final results or parlist
    for res in fh_in:
        res_strp = res.strip('\n')
        ent =res_strp.split(',')
        #print(ent)
        res_list.append(ent)
    print ('**List from ScoInpExt ready to feed to RevMapCoord function**')
    
    
    ###Out results
    rev_map_file = res_file.split('.')[0] + '_Revmap.csv'
    fh_out = open(rev_map_file, 'w')
    fh_out.write('miRNA,Target,Chr,Strand,Bind Site,Score,miRNASeq,TarSeq')
    
    ###Rev Map
    print ('***Reverse mapping coordinates**')
    for ent in res_list:
        print(ent)
        gene_name = ent[1]
        ##cleave_site = int(ent[3])
        bind_site = ent[3].split('-')
        
        ###Check whether gene is from reverse or posiive strand by memebr ship test on dictionary
        if gene_name in coord_dict_wat:
            print ('\nEntry: %s in positive strand: %s' % (ent[0:4],coord_dict_wat[gene_name]))
            geno_start = coord_dict_wat[gene_name][1]###Use for reverse mapping of postive genes
    
            #print('Looking for chr_id')
            chr_id = coord_dict_wat[gene_name][0]
            #print('chr_id found')
            strand = 'w' ## AVlilable in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_wat[gene_name][2] ## Gene type
            ##new_cleave_site = (int(geno_start)-1)+int(cleave_site)###1 is reduced to give correct positions
            new_bind_site_start = (int(geno_start)-1)+int(bind_site[0])
            new_bind_site_end = (int(geno_start)-1)+int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)
        else:
            print ('\nEntry: %s in reverse strand: %s' % (ent[0:4],coord_dict_crick[gene_name]))
            geno_end = coord_dict_crick[gene_name][1]###use for reverse mapping of negative genes
            #print('Looking for chr_id')
            chr_id = coord_dict_crick[gene_name][0]
            #print('chr_id found')
            strand = 'c' ## Available in dictionary also coord_dict_crick[gene_name][1]
            gtype =  coord_dict_crick[gene_name][2] ##Gene type
            ##new_cleave_site = (int(geno_end)+1)-int(cleave_site)###1 is added to give correct positions
            new_bind_site_end = (int(geno_end)+1)-int(bind_site[0])###As the sequence was reversed before TF and CL, their binding start and end direction has also changed - Verified-OK
            new_bind_site_start = (int(geno_end)+1)-int(bind_site[1])
            new_bind_site = '%s-%s' % (new_bind_site_start,new_bind_site_end)   
    
        
        rev_mapped_entry = ('%s,%s,%s,%s,%s,%s,%s,%s\n' % (ent[0],gene_name,chr_id,strand,new_bind_site,ent[2],ent[4],ent[5]))
        print (rev_mapped_entry)
        fh_out.write(rev_mapped_entry)
        
    
    fh_in.close()
    fh_out.close()    
    
    return rev_map_file




def main(gene_coords,res_file):
    RevMapCoord(gene_coords,res_file)
    
    
    
if __name__ == '__main__':
    main(gene_coords,res_file)
    print('***\nScript finished without error')
    sys.exit()

    

