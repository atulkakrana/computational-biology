#!/usr/local/bin/python3

###Script to extract phasiRNA from a cluster file for specific library
###Used phased loci name as in PARE validation results or coords in csv format to extract the phasiRNAs
###Usage

import os,sys
import difflib
import time

##########Settings###################

coordsfile = 'LDMAR.txt' #### File with either phased loci ID 'Phas-100_w_7:6353898:6354260' or coords 'w,7,6353898,6354260'
phasedID = 'N'### If given than phas ID is used to extract coords, if not than coords specified by user are used

coordsep ='.' ### Seprator used in the file

####User specifed coords in input file
chrcol = '2' ### As in excel format
stratcol = '3' ### As in excel format
endcol =  '4' ### As in excel format

clustfile = '2598.txt.score_p1e-07_sRNA_21_out.cluster'
phasiLenFilter = 'Y'
phassize = '21'

startbuff = 25 ### WHile extracting sequence through coords2FASTA a buffer is added to start, strat position in phased ID has this buffer added, so minus this buffer for better matching with real loci
getAbun = 'Y' ## Get only the abundance for sRNAs in a phased loci and not seqeucnes

def getClust(clustfile):
    outfile = clustfile+'phasi.csv' ### 2437.txt.score_p1e-07_sRNA_21_out.cluster
    fh_out = open(outfile,'w')
    
    fh_in = open(coordsfile)
    #fh_in.readline()
    
    fh_in2 = open(clustfile,'r')
    clusters = fh_in2.read().split('>')
    #coordsdict = {}
    for ent in fh_in:#### Given an entry in coords file
        print('\n\nEntry from input file being matched: %s' % (ent))
        coords = ent.split(coordsep)
        #print('coords:',coords)
        if phasedID == 'Y': ## Extract coords from phased ID
            loci = coords[0].split('_')##Phas-10_w_3:27574117:27574772
            phasID = loci[0]
            fusedcoords = loci[2].split(':')###3:27574117:27574772
            get_chr_id = fusedcoords[0]
            get_start = int(fusedcoords[1])+startbuff ## It should be added to reduce length of loci
            get_end = int(fusedcoords[2])+1 ##1 added because when opening a range using start and end, end number is not included in range
        else: ## User specified coords
            get_chr_id = coords[1]
            get_start = int(int(coords[2])+startbuff)
            get_end = int(coords[3])+1 ##1 added because when opening a range using start and end, end number is not included in range
            phasID = '%s_%s_%s' % (get_chr_id,get_start,get_end) ### Name to be used in output file
        get_value  = (list(range(int(get_chr_id+str(get_start)),int(get_chr_id+str(get_end)))))
        #print (get_value) ### OK till here
        
        ###Find matching cluster in phasifile
        for aclust in clusters[1:]:
            aclust_splt = aclust.split('\n')
            #print (aclust_splt[1:-1])###Last entry is always empty
            header = aclust_splt[0].split()
            #print (header)
            #print ('This is one cluster:\n %s \n' % aclust)
            clust_id = header[2]
            chr_id = header[6]
            start = header[10]
            end = int(header[12])+1 ##1 added because when opening a range using start and end, end number is not included in range
            #print (chr_id,start,end,'\n',block)
            value = (list(range(int(chr_id+str(start)),int(chr_id+str(end)))))
            #print ('Cluster:', (value))
            
            sm=difflib.SequenceMatcher(None,get_value,value)
            if round(sm.ratio(),2) >= 0.50: ### Get phasiRNA from this cluster
                print ('Matching cluster found:%s' % ''.join(header))
                #head = '>phasID_' % ()
                #fh_out.write()
                if getAbun == 'Y': ## Fetch just the locus abundance
                    locusabun = 0 ##Abundances of each locus
                    print('Getting the abundances of the phased loci')
                    for i in aclust_splt[1:-1]:### In one locus - header was the first entry of block and not required here,###Last entry is always empty
                        print (i)
                        phasient = i.split('\t')
                        phasiname = phasient[4].replace("|","_")
                        print(phasient,phasiname)
                        phasilen = phasient[6]
                        
                        if phasiLenFilter == 'Y': ### If tags filter is ON
                            if phasilen == phassize: ### Size specified in settings
                                phasiabun = phasient[7]
                                locusabun += phasiabun
                                
                        else:
                            phasiabun = phasient[7]
                            locusabun += phasiabun
                    fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (phasID,clust_id,locusabun))

                
                else: ## Fetch phasiRNAs
                    print ('Getting the phasi RNAs')
                    for i in aclust_splt[1:-1]:###Because header was the first entry of block and not required here,###Last entry is always empty
                        print (i)
                        phasient = i.split('\t')
                        #print(phasient)
                        #print(phasient[4])
                        phasiname = phasient[4].replace("|","_")
                        phasiseq = phasient[5]
                        phasiabun = phasient[7]
                        phasilen = phasient[6]
                        #print(phasilen)
                        if phasiLenFilter == 'Y': ### If tags filter is ON
                            if phasilen == phassize: ### Size specified in settings
                                fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun))
                        else: ##output all the tags
                            fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun))###phasiabun added after sequence because filtering in csv file will retain seq after applying filter on abundance
            #else:
            #    print('No matching phased loci found for: %s' % (ent))
        
    fh_in.close()
    fh_in2.close()
    fh_out.close()
                        

def main():
    getClust(clustfile)
        
        
        
if __name__ == '__main__':
    main()
    print('\n**script finished sucessfully**\n')
    sys.exit()
    
    
## v01 -> V02
##Added get just the abundance of the locus in provided library file
    
    

    
    
    
