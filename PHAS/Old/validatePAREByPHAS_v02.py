#!/usr/local/bin/python3

import sys,time,os,glob
import mysql.connector as sql



#### Settings
genomeDB = 'RICE_MSU7_genome' ##Used to identify strand of target
dataserver = 'raichu.dbi.udel.edu'

PAREpval = 0.5
PAREresType = 'M' ### M for MPPP and C for CL3
PAREres = 'Rice_1153P_1154P_1155P_1156P_RevMap_ALL.csv' ## Has to be revmapped file with 14th column as cleave site
#PAREres='sco_inp_ext_geno'

phase = 24 ##Which phase you are genomeDBusing for double validation - Give corresponding file below
PHASresType = 'P' ##J: jixian and P:Ping list file ##PP - Phas redundant script generated *.ping format
PHASres = '1193.txt.cluster.boundary.without.PARE.validation.list'
pVal = '0.005'


##Read phased loci file and make dictionary
def CreatePHASdict(PHASres,phase,pVal):
    ''' Reads the phased result file
    and creates a dictionary of -5/+7 phased locations as value
    '''
    
    fh_in = open(PHASres,'r')
    #fh_out = open('%s_Dict' % (PHASres),'w')
    PHASdict = {}
    entries = fh_in.readlines()
    
    if phase == 21:
        phaseList = [-105,-84,-63,-42,-21,0,21,42,63,84,105,126,147]
    elif phase == 24:
        phaseList = [-120,-96,-72,-48,-24,0,24,48,72,96,120,144,168]
    else:
        print('What the phase for phased loci - Please specify correctly')
        sys.exit()
    
    if PHASresType == 'P':
        for i in entries:
            if i.strip(): ## Remove an empty line from end file
                ent_splt = i.strip('\n').split('=')
                #print(ent_splt[0].split('|'))
                pval,phase,trash = ent_splt[0].split('|')
                chromo_start,end = ent_splt[1].split('..')
                chromo,start = chromo_start.split(':')
                akey ='%s-%s-%s' % (chromo.strip(),start,end)
                aval = [sum(i) for i in zip([int(start)]*11,phaseList)]
                #print ('Key:%s | Value:%s' % (akey,aval))
                if pval == pVal:
                    PHASdict[akey] = aval
                    #print ('Key:%s | Value:%s' % (akey,aval))
                    #fh_out.write('%s\t%s\t%s\t%s\t%s\tNONE\tNONE\n' % (phase,pval,chromo.strip(),start,end))##Chromosome has space before it which later gives error while key matching
                    #fh_out.write('%s\t%s\n' % (str(akey),str(aval).strip('[]')))
    
    elif PHASresType == 'PP': ##For processes/merged scenario - use *.ping file
        for i in entries:
            if i.strip(): ## Remove an empty line from end file
                phase,pval,chromo,start,end,trash1,trash2 = i.strip('\n').split('\t')
                print(phase,pval,chromo,start,end,trash1,trash2)
                akey ='%s-%s-%s' % (chromo.strip(),start,end)
                aval = [sum(i) for i in zip([int(start)]*11,phaseList)]
                if pval == pVal:
                    PHASdict[akey] = aval
    
    elif PHASresType == 'J':
        for i in entries:
            if i.strip(): ## Remove an empty line from end file
                ent_splt = i.strip('\n').split('\t')
                chromo = ent_splt[1].split('_')[1]
                start,end = ent_splt[4].split('-')
                akey = '%s-%s-%s' % (chromo.strip(),start,end)
                aval = [sum(i) for i in zip([int(start)]*11,phaseList)]
                PHASdict[akey] = aval
                print ('Key:%s | Value:%s' % (akey,aval))
                #fh_out.write('%s\t%s\t%s\t%s\t%s\tNONE\tNONE\n' % (phase,pval,chromo.strip(),start,end))##Chromosome has space before it which later gives error while key matching
                #fh_out.write('%s\t%s\n' % (str(akey),str(aval).strip('[]')))
    else:
        print('Please enter correct PHASED result type')
        pass
    fh_in.close()
    #fh_out.close()

        
    return PHASdict

##Reads CL and MPPP file to return alist
def PAREreader(PAREres,PAREpval):
    
    fh_in = open(PAREres, 'r')
    header = fh_in.readline() ## Header
    entries = fh_in.readlines()
    
    resList = [] ## miRNA,Target,cleavesite,whole entry
    totalCount = 0
    if PAREresType == 'M':
        for i in entries:
            #print (i)
            ent_splt = i.split(',')
            #print(ent_splt[13])
            if float(ent_splt[13]) <= PAREpval:
                cleaveSite = int(ent_splt[15])
                tarName = ent_splt[1]
                mirName = ent_splt[0]
                resList.append((mirName,tarName,cleaveSite,i.strip('\n')))
                #print(mirName,tarName,cleaveSite)
                totalCount +=1
    
    elif PAREresType == 'C':
        for i in entries:
            #print (i)
            ent_splt = i.split(',')
            if float(ent_splt[9]) <= PAREpval:
                mirName = ent_splt[0]
                tarName = ent_splt[1]
                cleaveSite = int(ent_splt[5])
                resList.append((mirName,tarName,cleaveSite,i.strip('\n')))
                totalCount +=1
        header = ('miR,TarName,Chr,Strand,bindSite,cleaveSite,Score,mirSeq,tarSeq,p-val,Small,Large,Ratio' % ()) ## as missing from sco_inp_ext
    
    print('\n\nTotal entries passed p-val i.e in list:%s made from file:%s' % (totalCount,PAREres))
    
    return resList,header

##Identify if CleaveSite in phasedloci
def validatePHAS(resList,PHASdict,phase,con,header):
    ''
    ''
    fh_out = open('%s_PHASmatched.csv' % (PAREres),'w')
    fh_out.write('%s,PHASLoci,PHASindex\n' % (header.strip('\n')))
    totalCount = 0 ##Counter for total entries in reslist
    phasedCount =0 ##Counter for matched entries from the reslist
     
    for ent in resList: ##ent format: miRNA,Target,cleavesite,whole entry
        totalCount +=1
        #print (ent)
        mirName = ent[0]
        tarName = ent[1]
        cleaveSite = int(ent[2])
        #print(mirName,tarName,cleaveSite)
        
        cur= con.cursor()
        #print(genomeDB,tarName.split('_up')[0].split('_down')[0])
        cur.execute("SELECT strand FROM %s.gene_master where gene like '%s'" % (genomeDB,tarName.split('_up')[0].split('_down')[0]))### Convert intergenic to gene name so as to get strand
        tarStrand = cur.fetchall()
        #print (tarStrand)
    
        if tarStrand[0][0] == 'c': ##Add 3 to cleave site
            #print("\nCrick strand:%s" % (tarStrand[0][0]))
            cleaveSite += 3    
        else: ## Do nothing to cleave site
            #print("\nWatson strand:%s" % (tarStrand[0][0]))
            pass
        
        ##PHAS LOCI MATCH
        for akey in PHASdict.keys():
            aval = PHASdict[akey]
            if cleaveSite in aval:
                print('miRNA:%s | Target:%s |CleaveSite:%s | Strand:%s @ PHAS:%s - %s phase: %d\n' % (mirName,tarName,cleaveSite,tarStrand[0][0],akey,aval,(aval.index(cleaveSite)+1-6))) ## 1 added to correct the index in human format and 6 added because start is at 6th postion
                fh_out.write('%s,%s,%s\n' % (ent[3],akey,(aval.index(cleaveSite)+1-6)))
                phasedCount += 1
                
            else:
                pass
        
    print('Total entries in PARE list:%s and matched:%s' % (totalCount,phasedCount))
    
    pass

##Connect to DB
def ConnectToDB(server, infile):
    
    ##infile values are '0' when you dont want to pulaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE bioinfo_data.mir_page_results FIELDS TERMINATED BY ','")
    
    print ('\nTrying to connect to mySQL server on %s' % (server))
    # Try to connect to the database
    try:
        con=sql.connect(host= server, user='kakrana', passwd='livetheday')###local_infile = 1 not supported yet so a table has to be updated on row basis
        print ('Connection Established\n')

    # If we cannot connect to the database, send an error to the user and exit the program.
    except sql.Error:
        print ("Error %d: %s" % (sql.Error.args[0],sql.Error.args[1]))
        sys.exit(1)

    return con


##MAIN
def main():
    PHASdict = CreatePHASdict(PHASres,phase,pVal)
    resList,header = PAREreader(PAREres,PAREpval)
    
    con = ConnectToDB(dataserver,0)
    validatePHAS(resList,PHASdict,phase,con,header)

if __name__ == '__main__':
    main()
    sys.exit()
    
    
##V01 -> v02 ## Feb-9th
##Addition of category to MPP results let to change in PAREreader indexes and output - Should be used on new result format
##Fixed if '_down' is encountered while making mysql querty to get the strand
## header for CL result from PAREredaer should be removed as added to latest version of CL itself

    

    
    
