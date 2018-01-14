#!/usr/local/bin/python3

###Script to extract phasiRNA from a cluster file for specific library
###Used phased loci name as in PARE validation results or coords in csv format to extract the phasiRNAs
###Usage

import os,sys
import difflib
import time

##########Settings###################

coordsfile = 'Final_PHAS_Loci_ALL_v2_HiConf.csv' #### File with either phased loci ID 'Phas-100_w_7:6353898:6354260' or coords 'w,7,6353898,6354260'
phasedID = 'N'### If given than phas ID is used to extract coords, if not than coords specified by user are used
header = 'Y' ## Input coordsFile has header or not

coordsep =',' ### Seprator used in the file
matchType = 'A' ## A: Accurate match of coo-ordinates OR F: Fuzzy match based on ratio

####User specifed coords in input file
chrcol = 5 ### As in excel format
startcol = 6 ### As in excel format
endcol =  7 ### As in excel format

clustfile = 'ALL_p5e-05_sRNA_21_out.cluster' ### with suffix *.cluster
phasiLenFilter = 'Y'
phassize = '21'

startbuff = 0 ### WHile extracting sequence through coords2FASTA a buffer is added to start, strat position in phased ID has this buffer added, so minus this buffer for better matching with real loci
getAbun = 'N' ## Get only the abundance for sRNAs in a phased loci and not seqeucnes

def getClustFuzzy(coordsfile,clustfile):
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
            end = header[12]
            #print (chr_id,start,end,'\n',block)
            value = (list(range(int(chr_id+str(start)),int(chr_id+str(int(end)+1))))) ## +1 to end because when opening a range using start and end, end number is not included in range
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

def getClustAccurate(coordsfile,clustfile,header):
    
    outfile = clustfile+'phasi.csv' ### 2437.txt.score_p1e-07_sRNA_21_out.cluster
    fh_out = open(outfile,'w')
    
    fh_in = open(coordsfile)
    if header == 'Y':
        fh_in.readline()
    
    fh_in2 = open(clustfile,'r')
    clusters = fh_in2.read().split('>')
    
    #### For each entry in coordsFile i.e. phased loci
    for ent in fh_in:
        print('\n\nEntry from input file being matched: %s' % (ent.strip('\n')))
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
            get_chr_id = coords[chrcol-1] ## -1 for converting to python format
            get_start = coords[startcol-1]
            get_end = coords[endcol-1] 
            phasID = '%s_%s_%s' % (get_chr_id,get_start,get_end) ### Name to be used in output file
            print ('This is the phased ID being processed: %s\n' % (phasID))

        ### Find matching cluster in phasifile
        matchClust = {} ## Dictionary to store clusters from different libraries and than select one based on combined abundance
        matchAbun = {} ## Dictionary to store matched cluster number and abundance
        matchClustNum = 0 ## Numer all matched cluster and choose one with highest abundance
        for aclust in clusters[1:]:
            aclust_splt = aclust.split('\n')
            #print (aclust_splt[1:-1])###Last entry is always empty
            header = aclust_splt[0].split()
            #print (header)
            #print ('This is one cluster:\n %s \n' % aclust)
            clust_id = header[2] ## Of no use in all library combined file for accurate match
            chr_id = header[6].replace('chr','')
            start = header[10]
            end = header[12]
            #print ('Cluster coords:',chr_id,start,end,'\n')
            
            if chr_id == get_chr_id and start == get_start and end == get_end: ### Get phasiRNA from this cluster
                print ('Matching cluster found:%s' % ''.join(header))
                matchClustNum += 1 ## One cluster matched, this is used for controlling the different clusters from different library that match the same phased loci
               
                if getAbun == 'Y': ## Fetch just the locus abundance
                    print('Getting the abundances of the phased loci')
                    akey = 'Match%s' % (matchClustNum)
                    abunSum = 0 ## variable to hold total abundance of in phase siRNAs, will be used to choose best cluster from multiple libraries
                    avalue = []
                    for i in aclust_splt[1:-1]:### In one locus - header was the first entry of block and not required here,###Last entry is always empty
                        print (i)
                        phasient = i.split('\t')
                        phasiname = phasient[4].replace("|","_")
                        print(phasient,phasiname)
                        phasilen = phasient[6]
                        
                        if phasiLenFilter == 'Y': ### If tags filter is ON
                            if phasilen == phassize: ### Size specified in settings
                                phasiabun = phasient[7]
                                abunSum += phasiabun
                                avalue.append((phasID,clust_id,abunSum))
                                
                        else:
                            phasiabun = phasient[7]
                            abunSum += phasiabun
                            avalue.append((phasID,clust_id,abunSum))
                        
                    matchClust[akey] = avalue
                    matchAbun[akey] = abunSum
                
                else: ## Fetch phasiRNAs
                    print ('Getting the phasi RNAs')
                    akey = 'Match%s' % (matchClustNum)
                    abunSum = 0 ## variable to hold total abundance of in phase siRNAs, will be used to choose best cluster from multiple libraries
                    avalue = []
                    for i in aclust_splt[1:-1]:### Because header was the first entry of block and not required here,###Last entry is always empty
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
                                abunSum += int(phasiabun)
                                avalue.append((phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun))
                                #fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun))
                        else: ##output all the tags
                            abunSum += int(phasiabun)
                            avalue.append((phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun))
                            #fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun))###phasiabun added after sequence because filtering in csv file will retain seq after applying filter on abundance
                    #print ('writing matchClust')
                    matchClust[akey] = avalue
                    matchAbun[akey] = abunSum
        
        ##Find best matched locus i.e. with most abundance out of all matched ones and write reults -- OR SHOULD WE TAKE TAG NORM ABUNDANCE THAT IS DIVIDE THE abun with no of tags and than choose
        
        if matchAbun: ## While testing all entries might not have match
            bestLocus = max(matchAbun,key=matchAbun.get) ## With maximum abundance out of all matched ones from different library
            #print ('matchClust',matchClust)
            locusToWrite = matchClust[bestLocus]
            print('Most abundant locus with phasiAbun: %s is %s' % (abunSum,locusToWrite))
            if getAbun == 'Y':
                for i in locusToWrite:
                    fh_out.write('>%s_Clust%s_%s_%s,%s\n' % (i[0],i[1],i[2],i[2]))
            else:
                for i in locusToWrite:
                    fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (i[0],i[1],i[2],i[3],i[4],i[5]))
        else:
            print("****No match found***")
        
    fh_in.close()
    fh_in2.close()
    fh_out.close()


def main():
    if matchType == 'F': ## a.k.a fuzzy match
        getClustFuzzy(coordsfile,clustfile)
    elif matchType == 'A': ## Accurate match
        getClustAccurate(coordsfile,clustfile,header)
    else:
        print("Please select the correct match mode in variable 'matchType'")
        sys.exit()
        
        
        
if __name__ == '__main__':
    main()
    print('\n**script finished sucessfully**\n')
    sys.exit()
    
    
## v01 -> V02
##Added get just the abundance of the locus in provided library file

##v02->v03
## Instead of ratio based  (fuzzy) match - get exact (accurate) match - Use a concatanated filr of all phased loci for exact match of cluster
## The fuzzy method needs update to calculate key and value as done in PhasRedundant_v01
    
    

    
    
    
