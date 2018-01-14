#!/usr/local/bin/python3

###Script to extract phasiRNA from a cluster file for specific library
###Used phased loci name as in PARE validation results or coords in csv format to extract the phasiRNAs
###Usage

import os,sys
import difflib
import time


##########Settings###################

clustfile = 'AllLibs.score_p0.001_sRNA_24_out.cluster' ## Cluster file from phasing analysis, coul dbe for one lib or concatanated file fro all libs

#### User specifed coords in input file
coordsfile = 'Final_24PHAS_Loci_ALL_v4.txt'            ## File with either phased loci ID 'Phas-100_w_7:6353898:6354260' or coords 'w,7,6353898,6354260'
head = 'Y'
coordsep ='\t'                              ## Seprator used in the file
phasedID = 'N'                              ## If given than phas ID is used to extract coords, if not than coords specified by user are used
chrcol = 5                                  ## As in excel format
startcol = 6                                ## As in excel format
endcol =  7                                 ## As in excel format

phasiLenFilter = 'Y'                        ## 'Y' then tags of phase length will be written from the cluster | 'N' - All tags will be written
minAbun = 10                                ## Minimum abundance of tag to be written
matchThres = 0.95                           ## Ratio of match required by the cluster to phased loci
phase = 24
startbuff = 0                               ## While extracting sequence through coords2FASTA a buffer is added to start, start position in phased ID has this buffer added, 
                                            ## so minus this buffer for better matching with real loci

def getClust(clustfile,coordsfile):
    
    outfile = clustfile+'thres_%s_phasi.csv' % matchThres                       ## 2437.txt.score_p1e-07_sRNA_21_out.cluster
    fh_out = open(outfile,'w')
    
    fh_in = open(coordsfile)
    if head == 'Y':
        fh_in.readline()
    
    fh_in2 = open(clustfile,'r')
    clusters = fh_in2.read().split('>')
    #coordsdict = {}
    
    phasCount = 0                                                               ## Total phased loci in file
    uniqMatchCount = 0                                                          ## Atleast one cluster present for one phased loci
    allMatchCount  = 0                                                          ## Total number of matched cluster
                                                            
    for ent in fh_in:##                                                         ## Given an entry in coords file
        print('\n\nEntry from input file being matched: %s' % (ent))
        coords = ent.split(coordsep)
        #print('coords:',coords)
        if phasedID == 'Y':                                                     ## Extract coords from phased ID
            loci = coords[0].split('_')                                         ## Phas-10_w_3:27574117:27574772
            phasID = loci[0]
            fusedcoords = loci[2].split(':')                                    ## 3:27574117:27574772
            get_chr_id = fusedcoords[0].replace("chr","").replace("Chr","")
            get_start = int(fusedcoords[1])+startbuff                           ## It should be added to reduce length of loci
            get_end = int(fusedcoords[2])+1                                     ## 1 added because when opening a range using start and end, end number is not included in range
            phasCount += 1
        
        elif phasedID == 'N': ## User specified coords
            # print("File with user specified columns will be used")
            get_chr_id = coords[chrcol-1]
            get_start = int(int(coords[startcol-1])+startbuff)
            get_end = int(coords[endcol-1])+1                                   ## 1 added because when opening a range using start and end, end number is not included in range
            phasID = '%s_%s_%s' % (get_chr_id,get_start,get_end)                ## Name to be used in output file
            print("Phased Loci: %s #######################################################" % (phasID))
            phasCount += 1

        else:
            print("Please input correct value for 'phasedID' or check your file")

        get_value  = (list(range(int(str(get_start)),int(str(get_end)))))
        #print (get_value) ### OK till here


        ################################################ Matching ############################################
        ######################################################################################################

        ## Find matching cluster in phasifile
        matchCount = 0                                                          ## Total maching clusters for a phased loci - if same cluster in multiple libraries
        finalMatchList = []                                                     ## Holds best cluster, from multiple libraries
        
        for aclust in clusters[1:]:
            tempMatchList = [] ## To hold resulst current matching cluster
            aclust_splt = aclust.split('\n')
            header = aclust_splt[0].split()
            #print (aclust_splt[1:-1]) ## Last entry is always empty
            #print ('This is one cluster:\n %s \n' % aclust)
            clust_id = header[2]
            chr_id = header[6].replace("chr","").replace("Chr","")
            start = header[10]
            end = int(header[12])+1 ##1 added because when opening a range using start and end, end number is not included in range
            #print (chr_id,start,end,'\n',block)
            value = (list(range(int(str(start)),int(str(end)))))
            #print ('Cluster:', (value))
            
            if int(get_chr_id) == int(chr_id):
                sm=difflib.SequenceMatcher(None,get_value,value)
                if round(sm.ratio(),2) >= matchThres: ### Get phasiRNA from this cluster
                    print ('\nMatching cluster found:%s' % ''.join(header))
                    matchCount +=1
                    
                    phasiCyc = 0 ## Stores phasing cycles
                    phasiSig = 0 ## Stores total abundance of phased sRNAs
                    
                    for i in aclust_splt[1:-1]:## Because header was the first entry of block and not required here, Last entry is always empty
                        print (i)
                        phasient = i.split('\t')
                        phasiname = phasient[4].replace("|","_")
                        phasiseq = phasient[5]
                        phasilen = int(phasient[6])
                        phasiabun = int(phasient[7])
                        # print(phasiname,phasiabun,phasiseq,phasilen)
                        tempMatchList.append((phasiname,phasiabun,phasiseq,phasilen))

                        if int(phasilen) == phase:
                            phasiCyc +=1
                            phasiSig += phasiabun
                    
                    # print(phasID)
                    print("Current Cycles:%s | Current sig. strength:%s" % (phasiCyc,phasiSig))
                    tempMatchList.append((phasiCyc,phasiSig,phasID,clust_id))
                    # print("\nCluster list:",tempMatchList)

                    ## Decide the best and remove other from list ##############################
                    ############################################################################
                    if finalMatchList: ## There exists a previosly matched cluster
                        exist_phasiCyc = finalMatchList[-1][0]
                        exist_phasiSig = finalMatchList[-1][1]
                        print("Existing Cycles:%s | Existing sig. strength:%s" % (exist_phasiCyc,exist_phasiSig))

                        if phasiCyc > exist_phasiCyc: ## New cluster has more cycles
                            del finalMatchList[0:]
                            finalMatchList = list(tempMatchList)
                            print("--- New cluster selected ---")


                        elif phasiCyc == exist_phasiCyc: ## Both have same cycles
                            if phasiSig > exist_phasiSig: ## New one has more total abundance of phased siRNAs
                                del finalMatchList[0:]
                                finalMatchList = list(tempMatchList)
                                print("--- New cluster selected ---")

                            else: ## Existing/old one had more abundance of phasiRNAs
                                pass
                        else: ## Existing/old one was long i.e. had more cycles
                            pass

                    else: ## This is the first cluster
                        finalMatchList = list(tempMatchList)

                    # print("\nFinal Match List:",finalMatchList)


        ## Write from final list ###########################################################
        ##############################################################################
        if finalMatchList:
            if phasiLenFilter == 'Y': ### If tags filter is ON
                phasID = finalMatchList[-1][2]
                clustID = finalMatchList[-1][3]
                # print("phasID:%s | ClustID:%s" % (phasID,clustID))
                for i in finalMatchList[:-1]: ## Last entry holds other info - phasiCyc, phasiSig,phasID and Clust_id
                    if int(i[-1]) == int(phase) and int(i[1]) > minAbun: ### Size specified in settings
                        fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (phasID,clustID,i[0],i[1],i[2],i[1]) )  ## phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun
                
            else: ## output all the tags
                phasID,clust_id = finalMatchList[-1][1:-1]
                print(phasID,clustID)
                for i in finalMatchList[:-1]: ## Last entry holds other info - phasiCyc, phasiSig,phasID and Clust_id
                    if int(i[1]) > minAbun: ### Size specified in settings
                        fh_out.write('>%s_Clust%s_%s,%s\n%s,%s\n' % (phasID,clustID,i[0],i[1],i[2],i[1]) )  ## phasID,clust_id,phasiname,phasiabun,phasiseq,phasiabun


        allMatchCount += matchCount ## Add the mached cluster for each entry
        if matchCount > 0 :
            uniqMatchCount+=1
        else:
            print("####################################################################### NO MATCH FOUND\n")

    print ("\nSUMMARY: Phased loci in input file:%s | Loci match threshold: %s |Uniq matched cluster: %s | Total matched clusters found:%s" % (phasCount,matchThres,uniqMatchCount,allMatchCount))
    print("NOTE: If matched clusters more then phased loci that means same cluster was present in different libs\n")
    # print("NOTE: Don't forget to uniq the miRNAs")
    fh_in.close()
    fh_in2.close()
    fh_out.close()                 

def main():
    getClust(clustfile,coordsfile)

if __name__ == '__main__':
    main()
    print('\n**script finished sucessfully**\n')
    sys.exit()

#### v0.1 -> v0.3
## updated phas matching with chromosome not part of value instead matched before any ratio computation
## If clusters for all libs are concatanated and same phased locus matches clusters from different libs then the longest phase cycle is slected and is phase cycle are same most abundant is selcted