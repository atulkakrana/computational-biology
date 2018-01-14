#!/usr/local/bin/python3

import os,sys,operator,time,datetime,string,subprocess
import mysql.connector as sql
from collections import Counter
from multiprocessing import Process, Queue, Pool
import os.path

## This script is written for to fetch two max tag from miRNAs loci

#### USER SETTINGS ######

precursorFa = "MasterNormal.v3.21-22Priv.fa" ## File with precursor seqeunces, for foldback - Name and number will correspond to coordsF file
coordsF     = "MasterNormal.v3.21-22Priv.txt" ## File with precursor coordinates, and miRNA length - Name and number will correspond to precursorF
mirTagsF    = "aspa.mirTags.txt"  ## Tab separated file with miRNA names (w/o 5p or 3p) in first column matching above coords file, and seqeunces in second columns - This file is used to detect if the precursor had miRNA or not. This file could have two miRNAs of same name one for 5p and one for 3p, while coordsF and precursorF will have just one miRNA precusor entry.

mname   = 64        ## Excel format
mchr    = 67        ## Excel format
mstart  = 68        ## Excel format
mend    = 69        ## Excel format
mstrand = 71        ## Excel format
msize   = 12        ## Excel format
delim       = "\t"          ## Delimiter for coordsF
msizes      = [21,22]    ## miRNAs for these length to be used for analysis
head        = "Y"

libType     = 1             ## 0: Lib_ids (4518) | 1: lib_code ('leaf_1')
excludeLibs = []            ## If you wish to exclude a few libs then enter lib_id here, your tag position summary table shoul dhave libs ids and not lib codes

ntags       = 2            ## 1: report ratios based on most abundant tag 2: on basis of two topmost abundant tags
minAbun     = 50

### DEVELOPER SETTINGS ####
dataserver      = "tarkan.ddpsc.org"
msDB            = "kakrana"
tagPosTable     = "ASPARAGUS_privPHAS_sRNA_TagPosSummReNorm"
tagSizeBuffer   = 1 ## Tags of size +/-buffer to miRNA will be used in abundance calculation, because miRNA-3p are usally one nt short 
gap         = 12        ## Default:12 (from emboss website)
match       = 3         ## Default:3
mismatch    = -4        ## Default:-4
threshold   = 24        ## Default:50 | 24-selected after finding suxh cases in Maize
maxRepLen   = 2000      ## Default: 2000

#### Functions ####

def miRmax(coords):
    '''This function computes ratio of top one and top two sRNA tags
    from miRNA precursor, against all tags in that precursor'''

    ### GET THE BASIC INFO ######################
    con                 = ConnectToDB(dataserver)
    cur                 = con.cursor()
    lib_col             = prepareQuery(excludeLibs,cur)
    fastaList,fastaDict = fastaReader(precursorFa)
    mirtagsD            = mirTagsParse(mirTagsF)
    #############################################

    #### OutFile ################################
    outFile = "mirFetchMax.txt"
    fh_out  = open(outFile,'w')
    fh_out.write("miRNA\tmaxTag\tmaxtagAbun\tmaxtag2Abun\ttagsSum\tmaxratio\tmaxratio2\n")
    ##########################################################################

    ### Compute ratios of max tags in miRNA precursors ######################
    acount = 0 ### Counts precursors that didn't show foldback
    bcount = 0 ### COunts precursors that do not have any sRNA in their foldback
    ccount = 0 ### Counts precursors that ahve foldback and mapped sRNAs, but expected miRNA not found
    dcount = 0 ### All precursors tested
    for i in coords:
        aname,achr,astart,aend,astrand,alen = i
        print("\nEnt:",aname,achr,astart,aend,astrand,alen,"####################################")
        dcount+=1

        #### Foldback Coords #################################################
        aseq            = fastaDict[aname]
        outseq,outinv   = IRchecker(aseq,aname)
        IRcoords        = IRparser(outseq,outinv,aname)
        print("IRcoords",IRcoords)
        if IRcoords:
            IRname,score,matches,perc,gaps,alignLen,start5,end5,start3,end3,loop = IRcoords[0]
        else:
            acount +=1
            print("Precursor shows no foldback at given score thershold:%s - skipping precursor" % (threshold))
            continue

        #### Normalize IR coords to get genomic coords########################
        normstart5  = int(start5) + astart
        normend5    = int(end5)   + astart
        normstart3  = int(start3) + astart
        normend3    = int(end3)   + astart
        print("Normalized IR coords:",normstart5,normend5,normstart3,normend3)

        #### Query to get sRNAs from foldback #################################
        minlen = alen-tagSizeBuffer ## If miRNA is 5p then this length corresponds to 3p
        maxlen = alen+tagSizeBuffer ## If miRNA is 3p then this length corresponds to 5p
        cur.execute("select tag,%s from %s.%s where chr_id = %s and strand = '%s' and (len between %s and %s) and (position between %s and %s)" % (lib_col,msDB,tagPosTable,str(achr),str(astrand),str(minlen),str(maxlen),str(normstart5),str(normend5)))####
        temp5        = cur.fetchall()
        
        cur.execute("select tag,%s from %s.%s where chr_id = %s and strand = '%s' and (len between %s and %s) and (position between %s and %s)" % (lib_col,msDB,tagPosTable,str(achr),str(astrand),str(minlen),str(maxlen),str(normstart3),str(normend3)))####
        temp3        = cur.fetchall()
        
        # print("These are the tags from 5'arm:",temp5)
        # print("These are the tags from 3'arm:",temp3)

        temp = temp5+temp3
        # print("These are the tags from foldback:",temp)
        print("Total tags from foldback:",len(temp))

        #### Check - If miRNA was expressed in provided libraries ####################
        ##############################################################################
        if not temp:
            bcount+=1
            print("No tags found from foldback - skipping precursor")
            continue
        else:
            mirtags     = mirtagsD[aname] ### List of miRNAs related to this precuror i.e 3p, 5p 
            mirFound    = [] ### List to store miRNAs found in sRNAs from foldback
            for ent in temp:
                atag = ent[0]
                asum = sum(ent[1:-1])
                if (asum >= minAbun) and (atag in mirtags):
                    mirFound.append(atag)
                else:
                    pass

            if len(mirFound) == 0:
                print("No known miRNA identified from precursor foldback - skipiing precursor",mirFound)
                ccount+=1
                continue
                # sys.exit()
            else:
                print("Atleast one miRNA found at abundance threshold - processing precursor",mirFound)
                pass

        #### Collapse tag abundances tag from all libs to total ######################
        ###############################################################################
        precursorTags   = []    ## Temp list that stores merged abundances from all libs
        tagSet          = set() ## Maintain a uniq set, redundant tags found, problem with how tag position summary is computed
        allsum          = 0
        for ent in temp:
            atag = ent[0]
            if atag not in tagSet:
                asum = sum(ent[1:-1])
                # print(atag,asum)
                precursorTags.append((atag,int(asum)))
                tagSet.add(atag)
            else:
                print("Duplicate tag from same locus exists")
                # sys.exit()
                pass

        #### Get most abundant tags and total abundance of all tags in foldback ########
        ################################################################################
        precursorTags_s = sorted(precursorTags,key=operator.itemgetter(1),reverse=True) ## CHR_ID, strand and start 
        max1 = precursorTags_s[0]

        if len(precursorTags_s) > 1: ## In few cases there is no other tag of miRNA size, so max2 set to 0
            max2 = precursorTags_s[1]
        else:
            max2 = ["None",0]
        tagsSum = sum(i[1] for i in precursorTags_s)
        print("Max tag:",max1)
        print("Second max:",max2)
        print("Sum of all tags:",tagsSum)

        #### Compute rations and write results ########################################
        ############################################################################### 
        if ntags == 1:
            aratio1 = round(max1[1]/tagsSum,2)
            print('Res entry:',aname,max1[0],max1[1],max2[1],tagsSum,aratio1)
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (aname,max1[0],str(max1[1]),str(max2[1]),str(tagsSum),str(aratio1))) ### Name of miRNA entry, max tag, maxtag abun, maxtag2 abun, sum of all tags and aratio
        elif ntags == 2:
            aratio1 = round(max1[1]/tagsSum,2)
            aratio2 = round((max2[1]+max1[1])/tagsSum,2)
            print('Res entry:',aname,max1[0],max1[1],max2[1],tagsSum,aratio1,aratio2)
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (aname,max1[0],str(max1[1]),str(max2[1]),str(tagsSum),str(aratio1),str(aratio2))) ### Name of miRNA entry, max tag, maxtag abun, maxtag2 abun, sum of all tags and aratio
        else:
            print("Please input correct value for ntags")
            sys.exit()

    fh_out.close()
    print("\nPrecursor processed:%s | No foldback:%s | No sRNA:%s | Not expected miRNA:%s" % (dcount,acount,bcount,ccount))

    return outFile

def coordsParser(coordsF):
    '''parses input PHAS file'''

    print("\nFunction: coordsParser")

    fh_in   = open(coordsF,'r')
    if head == 'Y':
        fh_in.readline()
    coordsRead = fh_in.readlines()

    coordsList = []
    for i in coordsRead:
        print("Ent:",i)
        ent     = i.split(delim)
        aname   = ent[mname-1]
        achr    = int(ent[mchr-1].replace("chr",""))
        astart  = int(ent[mstart-1])
        aend    = int(ent[mend-1])
        asize   = int(ent[msize-1])
        astrand = ent[mstrand-1].translate(str.maketrans("+-","wc"))
        print("Phas Name %s | chr:%s | start:%s | end:%s | astrand:%s | Len:%s" % (aname,achr,astart,aend,astrand,asize))
        if asize in msizes:
            coordsList.append(((aname,achr,astart,aend,astrand,asize)))

    print("A list of query features prepared with entries:%s" % len(coordsList))
    print("Exiting function - coordsParser\n")
    
    return coordsList

def ConnectToDB(server):
    
    ##infile values are '0' when you dont want to pulaod data from local file and '1' when you wish to upload data by local file
    ##EX:con=sql.connect(host= server, user='kakrana', passwd='livetheday', local_infile = infile)
    ##Now later in script you can
    ##cur.execute("LOAD DATA LOCAL INFILE './scoring_input_extend2' INTO TABLE kakrana_data.mir_page_results FIELDS TERMINATED BY ','")
    
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

def fastaReader(fastaFile):
    
    '''Cleans FASTA file - multi-line fasta to single line, header clean, empty lines removal'''

    print("\nFUNCTION - fastaReader")
    ## Read seqeunce file
    print ('-Reading "%s" FASTA file' % (fastaFile))
    fh_in       = open(fastaFile, 'r')
    fasta       = fh_in.read()
    fasta_splt  = fasta.split('>')
    acount      = 0 ## count the number of entries
    empty_count = 0

    fastaList = [] ## Stores name and seq for fastFile
    fastaDict = {} ## Stores name as key and seq as value

    acount +=1
    for i in fasta_splt[1:]:
        acount  +=1
        ent     = i.split('\n')
        name    = ent[0].split()[0].strip()
        seq     = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        alen    = len(seq)
        fastaList.append((name,seq,alen))
        fastaDict[name] = seq

    print("-Total entries in phased fastaFile:%s" % (str(acount)))
    print("-fastaList generated with %s entries\n" % (str(len(fastaList))))
    print("-fastaDict generated with %s entries\n" % (str(len(fastaDict))))
    # print("-Length for %s:%s" % (name,alen))
    time.sleep(2)

    return fastaList,fastaDict

def IRchecker(combSeq,combName):
    '''
    This function takes a seqeunce, and runs einverted and reports back
    coordinates
    '''

    print("\nFUNCTION - IRchecker")
    seq         = combSeq
    name        = combName
    tempInput   = "tempSeq.fa"
    fh_out      = open(tempInput,'w')
    fh_out.write('>%s\n%s\n' % (name,seq))
    fh_out.close()

    outseq = "%s.fa.temp" % name
    outinv = "%s.inv.temp" % name

    retcode = subprocess.call(["einverted", "-sequence", tempInput, "-gap", str(gap), "-threshold", str(threshold), "-match",str(match),"-mismatch", str(mismatch), "-maxrepeat",str(maxRepLen), "-outfile",outinv, "-outseq",outseq ])

    if retcode == 0:## The bowtie mapping exit with status 0, all is well
        print('\n+einverted for %s\n' % (name) )
    else:
        print('Something wrong happened while running einverted for sequence: %s - - Debug for reason' % (name))
        sys.exit()

    ## Cleanup entry specifc FASTA file
    # if os.path.exists(tempInput):
    #     os.remove(tempInput)

    return outseq,outinv

def IRparser(outseq,outinv,combName):
    '''
    parse results and delete input files
    '''
    print("\nFUNCTION - IRparser")
    print ("+Parsing %s results" % (outinv))

    ## File to record results
    # resOut= "RES_%s.csv" % (combName)
    # fh_out = open(resOut,'w')
    IRcoords = [] ## Store coordinates for this inverted repeat
    # fh_out.write("EntryName,Score,Matches,Perc,Gaps,AlignLen,5'Start,5'End,3'start,3'end,Loop\n")

    ## Get list of files
    afile = outinv
    name = combName
    # print("Seqeunce being parsed:%s" % (afile))
    try:
        if os.stat(afile).st_size > 0:
            print ("+Sequence %s - Results !!!!" % (afile.split(".")[0]))
            fh_in = open(afile,'r')
            invs = fh_in.read().split("\n\n")
            # print("Empty line splitted:",invs)
            
            for i in invs:
                # print('\nInverted:',i.strip('\n'))
                invLines = i.strip('\n').split('\n')
                # print(invLines)
                
                resBlock_splt = invLines[0].split(":")
                # print(resBlock_splt)
                score = resBlock_splt[1].split()[1]
                # print(resBlock_splt[1].split())
                matchesInfo,gapsInfo = resBlock_splt[2].split(",")
                # print(matchesInfo,gapsInfo)
                gaps = gapsInfo.split()[0]
                
                # print(matchesInfo.strip().split())

                matches = matchesInfo.strip().split()[0]
                matched,total = matches.split("/")
                alignLen = int(total)*2
                perc = round(int(matched)/int(total),2)

                # matches,garbage1,perc,garbage2 = matchesInfo.strip().split()
                # alignLen = int(matches.split("/")[1])*2
                # print(score,matches,gaps)

                arm5 = invLines[1].strip()
                start5,seq5,end5 = arm5.split(" ")

                arm3 = invLines[3].strip()
                end3,seq3,start3 = arm3.split(" ")

                loop = int(start3)-int(end5)

                IRcoords.append((name,score,matches,str(perc),gaps,str(alignLen),start5,end5,start3,end3,loop))
                # fh_out.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (afile.split(".")[0],score,matches,str(perc),gaps,str(alignLen),start5,end5,start3,end3,loop)) ## name of file,score,match,gaps

            #### Clean up ### 
            garbage = [afile for afile in os.listdir('./') if afile.endswith (('.fa.temp','.inv.temp'))]
            for afile in garbage:
                print("+Deleting %s" % (afile))
                os.remove(afile)

        else:
            # print ("Sequence %s - No results" % (afile.split(".")[0]))
            pass

    except OSError:
        print("+No result file for sequence %s found - Please check" % (afile))
        print("+System will exit")
        sys.exit()


    return IRcoords

def mirTagsParse(miRtags):
    '''Reads miR tags file which has miRNA name matching fasta file and coords file
    in first column and tags in second column. One miRNA name can have multiple tags 
    i.e. from 3p and 5p, or from different precursors. All these tags should be recorded
    as value tuple with miRNA name as key.'''

    print("\nFunction: mirTagsParse")

    # mirtagsL = [] ### List to store results
    mirtagsD = {} ### DIct to store results

    fh_in = open(miRtags,'r')
    fh_in.readline()
    fileRead = fh_in.readlines()

    #### Get Uniq name list ####
    uniqmirs = set()
    for i in fileRead:
        ent     = i.strip('\n').split("\t")
        aname   = ent[0]
        uniqmirs.add(aname)

    for aname in uniqmirs:
        print("-caching entries for %s" % (aname))

        tempList = [] ### MiRNA tags
        for i in fileRead:
            ent     = i.strip('\n').split("\t")
            bname   = ent[0]
            atag    = ent[1]
            
            if bname == aname:
                tempList.append (atag)
            else:
                pass

        print("-%s tags found for %s" % (len(tempList),aname))
        mirtagsD[aname] = tempList

    print("\n+Tags dict prepared for %s miRNAs" % (len(uniqmirs)))

    print("Exiting function - mirTagsParse\n")

    return mirtagsD

def prepareQuery(excludeLibs,cur):

    ### Prepare query of libs #################

    ### get column names
    columns = [] ## empty list
    cur.execute("describe %s.%s" % (msDB,tagPosTable))
    tablefield = cur.fetchall()
    # print("\nTable fields:",tablefield)
    for i in tablefield:
        col = i[0]
        columns.append(col)

    libs    = columns[8:-2] ## norm_sum and max_norm are last two columns
    lastTwo = columns[-2:] ## norm_sum and max_norm
    print("\nLibs:",libs)

    if excludeLibs:
        print("\n\nLibs specifed in excludeLibs %s will be skipped\n\n" % (excludeLibs))
        selectLibs = [] ## Columns excluding the unwanted libraries
        excludeLibs_s = [str(i) for i in excludeLibs] ## Converting all entries in exclude list to string fro matching below
        
        for i in libs:

            ## Check if user made mistake in libType - as that will give results for all entries
            if type(i) is int and libType == 1:
                print("You seem to have input lib_id and chosen wrong libType")
                print("Check libType and excludeLibs match - Script will exit now")
                sys.exit()
            elif type(i) is str and libType == 0:
                print("You seem to have input lib_code and chosen wrong libType")
                print("Check libType and excludeLibs match - Script will exit now")
                sys.exit()
            else:
                print("All seems well")
                pass

            ## Get lib_if if libType == 0
            if libType == 0:
                lib = i.split('_')[1] ## Closed to match lib_codes too - Columns in tag position summary - 'lib_2499' or 'max_norm'
                print(lib)
            else:
                lib = i
                pass

            ### Filter libraries
            if str(lib) not in excludeLibs_s: ## Tested OK
                selectLibs.append(i)
            else:
                print("excluded:",i)
                # sys.exit()
                pass
        
        finalLibs = selectLibs+lastTwo ## Norm_Sum and max_norm are included

    else:
        finalLibs = libs+lastTwo ## Norm_Sum and max_norm are included

    # print("finalLibs:%s" % (finalLibs))
    
    lib_col =",".join(str(x) for x in finalLibs)### Manually mentioned strating lib column - should work on all tag position summary tables
    
    print("Library Columns:",lib_col)
    # queryLibs = 'SUM(%s)' % lib_col.replace(",","),SUM(")
    # sumLibs = "%s" % lib_col.replace(",","+") ## From phasiMax
    # queryLibs = "%s" % lib_col.replace(",",",") ## From phasiMax
    # print("\nThese are sumLibs:",sumLibs)
    # print("\nThis is query Libs:",queryLibs)
    # sys.exit()

    return lib_col


def main():

    coordsL = coordsParser(coordsF)
    resList = miRmax(coordsL)

    return None

if __name__ == '__main__':
    main()
    print("\nScript finished sucessfully")
    sys.exit()


### v01 #######
### Takes precursor coords, gathers sRNAs, and select top two tags to compute ratios against allsRNAs
### Only uniq tags used


### v01 -> v02 [major][stable]
### Checks for foldback of prursor using FASTA file
### Gets sRNAs from foldback region
### Checks sRNAs for mature miRNAs of same name


###
### 3p are usually int shorter than miRNAs, so query to get abundance need to be modified so that +/-1 nt of miRNAs length
### should be used in query