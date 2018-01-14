#!/usr/local/bin/python3

import os,sys,operator,time,datetime,string,subprocess
import mysql.connector as sql
from collections import Counter
from multiprocessing import Process, Queue, Pool
import os.path

###### USER SETTINGS ####################
miRBlastF   = "Final_Zm-Os-Ao.mature.mod.uniq.txt"      ## TXT file of filtered BLAST results
miRBlastFa  = "lilium.tagpos.tags.fa"                   ## FASTA file of all miRNA candidates, identified through miRAnno
transFa     = "Trinity.fasta"                           ## All tags, that need to be mapped to infer top two tags on precursors i.e. Trinity FASTA
 
indexStep  = 0                                        ## Index of transcriptome required to map Tags, mandatory for first run
includeLibs = ["Lilium_leaf","Lilium_4mm_an","Lilium_5mm_an","Lilium_6mm_an","Lilium_8mm_an","Lilium_10mm_an"]
##### ADVANCED SETTINGS ################
msDB        = "kakrana"
tagPosTable = "DAYLILY_priv_sRNA_TagPosSummNorm"
dataserver  = "tarkan.ddpsc.org"
numProc     = 24
gap         = 12        ## Default:12 (from emboss website)
match       = 3         ## Default:3
mismatch    = -4        ## Default:-4
threshold   = 24        ## Default:50 | 24-selected after finding suxh cases in Maize
maxRepLen   = 2000      ## Default: 2000
miRsize     = [20,21,22]## candidate miRNAs of this length will be considered for analysis
sRNAsize    = [20,21,22,23,24] ## sRNAs of these sizes will be used to compute ratios
ntags       = 2
########################################

def inferCandidates(miRBlastF,miRBlastFa):
    '''This module takes the filtered BLAST results of miRNA query against all tags of species, 
    maps seqeunces for these results to transcripts and gives back candiates'''

    ### Get the miRNA candidates from miRblast
    miRcandsL,miRcandsD = blastReader(miRBlastF)

    ### Fetch seqeunces of miRNA candidates from FASTA file
    fastaL,fastaD   = fastaReader(miRBlastFa)
    candsFa         = "tempTags.fa" ## Write canidate miRNA tags for mapping to transcriptome
    fh_out          = open(candsFa,'w')
    
    for cand in miRcandsL:
        candName = cand[1]
        print("Fetching tag for miRcand:%s" % (candName))
        candTag  = fastaD[candName]
        print("-seqeunce found for %s" % (candName))
        fh_out.write(">%s\n%s\n" % (candName,candTag))
    fh_out.close()

    ### Map candidate miRNAs to transcritome so as to get precursors
    fetchedLibs = ["tempTags"] ### OPEN after test
    mapLibs(transFa,fetchedLibs,indexStep) ## OPEN after test
    alib            = fetchedLibs[0]        ## For comaptibility with map2Dict functions
    mapDict         = map2Dict(alib)        ## Fetches mapping results with transcript as key
   
    ### Get precursor candidates from mapped miRNA cand file, and their foldback results
    candsFastaL,candsFastaD     = fastaReader(transFa)  ## Read transcripts file
    
    ## For every precursor get it's foldback, filter on foldback length and miRNA position on stem
    precFa      = "precursors.fa"        ## Output file for miRNA candidates
    fh_out      = open(precFa,'w')    
    precSumm    = "precursorsSumm.txt" 
    fh_out2      = open(precSumm,'w')
    header      = "miRcand\tputativeName\tmiRstrand\tmiRseq\tmiRlen\tmiRhits\tmappos\tacand\tprecursorSeq\tprecursorLen\tstatusFlag\tprecursorName\tscore\tmatches\tstr(perc)\tgaps\talignLen\tstart5\tend5\tstart3\tend3\tloop"
    fh_out2.write("%s\n" % (header))
    
    precSummL = [] ### List to store all results
    # precSummD = {} ### Dict to store all results, Unique tag is key and entry is value
    for acand in mapDict.keys():     ## Values are bascially miRNA candidate and mapping details
        atransSeq = candsFastaD[acand]
        candRes,IRcoords = predictPrecursor(acand,atransSeq,mapDict,miRcandsD)

        if candRes:
            for miRcand in candRes:
                finalStatus = miRcand[-1]

                if finalStatus == "Y":
                    fh_out.write(">%s\n%s\n" % (acand,atransSeq))
                    fh_out2.write("%s\t%s\n" % ('\t'.join(str(x) for x in miRcand),'\t'.join(str(x) for x in IRcoords[0])))
                    tempResL = miRcand + IRcoords[0]
                    precSummL.append(tempResL)
                else:
                    print("Candidate miRNA and precursor shows no evidence of being a miRNA")
                    pass
        else:
            print("Precursor does not foldback")
            pass
    
    fh_out.close()
    fh_out2.close()

    print("Snippet of precSummL",precSummL[1:5])
    return precSummL,precFa,header

def ratioComputer(precSummL,precFa,header):
    '''
    This function maps all sRNAs to final precursors, and computes ratio of two one or two tags 
    based on tags mapped in foldback region
    '''
    con = ConnectToDB(dataserver)
    cur = con.cursor()
    outFile = "mirFetchMaxLocal.txt"
    fh_out  = open(outFile,'w')
    fh_out.write("miRNAPrec\tmaxTag\tmaxtagAbun\tmaxtag2Abun\ttagsSum\tmaxratio\tmaxratio2\n")
    
    ## Make a dict of IR coords
    IRcoordsD   = {}
    tempSet     = set()
    for i in precSummL:
        aprec = i[7]
        if aprec not in tempSet:
            IRcoordsD[aprec] = i
            tempSet.add(aprec)
    print("Dictionary of IRcoords prepared with entries:%s" % (len(IRcoordsD)))


    ### Map sRNAs to precursor fasta
    fastaIn   = "%s" % precFa.rpartition(".")[0]
    print("These is the precursor name for index generations:%s" % (fastaIn))
    fetchedLibs = [miRBlastFa.rpartition(".")[0]] ### OPEN after test
    mapLibs(precFa,fetchedLibs,1) ## OPEN after test
    alib            = fetchedLibs[0]        ## For comaptibility with map2Dict functions
    mapDict         = map2Dict(alib)        ## Fetches mapping results with transcript as key

    ### For every precursor (key), use IR coords to compute ratio
    for aprec in mapDict.keys():
        print("\n\nComputing ratio for:%s" % (aprec))
        coords  = IRcoordsD[aprec]
        start5  = int(coords[17])
        end5    = int(coords[18])
        start3  = int(coords[19])
        end3    = int(coords[20])
        print("Foldback coords:",start5,end5,start3,end3)
        
        mapTags = mapDict[aprec]
        lib_col = ','.join(x for x in includeLibs)
        maxLen  = max(sRNAsize)
        minLen  = min(sRNAsize)
        print("Libs",lib_col)

        ###### GET TAG ABUNDANCES ###############
        #########################################

        validTags   =[] ### List of tags from foldback
        tagSet      = set() ## Temporarily store tags to unique them
        for i in mapTags:
            print("Mapped tag ent:",i)
            mappos = i[-1]
            maptag = i[2]
            
            if len(maptag) in miRsize: ## Considering only tags of specified length, similar to phasiMAx where just 21 or 24-nt tags are used
                if mappos >= start5 and mappos <= end5:
                    print("-Position on 5' arm")
                    cur.execute("select tag,%s from %s.%s where tag = '%s' and len between %s and %s" % (lib_col,msDB,tagPosTable,maptag,minLen,maxLen))
                    temp        = cur.fetchall() ## Many hits of same tag, select just one
                    print("-These are the tags:",temp[0])
                    if maptag not in tagSet:
                        validTags.append(temp[0])
                        tagSet.add(maptag)

                    ### Fetch abundance, and add tag and abun to list
                elif mappos >= start3 and mappos <= end3:
                    print("-position on 3' arm")
                    cur.execute("select tag,%s from %s.%s where tag = '%s' and len between %s and %s" % (lib_col,msDB,tagPosTable,maptag,minLen,maxLen))
                    temp        = cur.fetchall()
                    print("-These are the tags:",temp[0])
                    if maptag not in tagSet:
                        validTags.append(temp[0])
                        tagSet.add(maptag)
                else:
                    print("-position not in stem")
                    aflag = 0 
                    pass
            else:
                ## To make the ratios comaprable with fetchMax script, only tags of expected sizes used
                ## i.e miRNA query is 20,21, and 22nt, considering this will cover 5p and 3p, plus noise
                print("The tag size is incomaptible with current analysis")


            ### Check to see if foldback has any tags mapped to its stem
            if len(validTags) == 0:
                print("No tags of desired length mapped to stem of precursor %s" % (aprec))
                continue
                # sys.exit()
            else:
                pass

        ####### COMPUTE RATIOS ###############
        ######################################
        ### Sum abundances from multiple libraries for each tag, sort on abundances, and compute ratios
        precursorTags = [] ## Store tags with their summed abundnaces
        for tag in validTags:
            atag = tag[0]
            asum = sum(tag[1:-1]) ## Sum of abundnace from all libraries
            precursorTags.append((atag,asum))
        precursorTags_s = sorted(precursorTags,key=operator.itemgetter(1),reverse=True) ## CHR_ID, strand and start

        ### Get max tags, and total sum of all tags on this precursor
        max1 = precursorTags_s[0]
        if len(precursorTags_s) > 1: ## In few cases there is no other tag of miRNA size, so max2 set to 0
            max2 = precursorTags_s[1]
        else:
            max2 = ["None",0]
        tagsSum = sum(i[1] for i in precursorTags_s)
        
        print("Max tag:",max1)
        print("Second max:",max2)
        print("Sum of all tags:",tagsSum)

        ## Compute rations and write to file
        if ntags == 1:
            aratio1 = round(max1[1]/tagsSum,2)
            print(aprec,max1[0],max1[1],max2[1],tagsSum,aratio1)
            fh_out.write("%s\t%s\t%s\t%s\t%s\n" % (aprec,max1[0],str(max1[1]),str(tagsSum),str(aratio1))) ### Name of miRNA entry, max tag, maxtag abun, maxtag2 abun, sum of all tags and aratio
        elif ntags == 2:
            aratio1 = round(max1[1]/tagsSum,2)
            aratio2 = round((max2[1]+max1[1])/tagsSum,2)
            print(aprec,max1[0],max1[1],max2[1],tagsSum,aratio1,aratio2)
            fh_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (aprec,max1[0],str(max1[1]),str(max2[1]),str(tagsSum),str(aratio1),str(aratio2))) ### Name of miRNA entry, max tag, maxtag abun, maxtag2 abun, sum of all tags and aratio
        else:
            print("Please input correct value for ntags")
            sys.exit()


        ## Use valid tags to find top two tags and ratios


    return None

def blastReader(miRBlastF):
    '''Read all three blast results and give back a transcript'''

    fh_in       = open(miRBlastF,'r')
    header      = fh_in.readline().strip("\n")
    fileRead    = fh_in.readlines()
    fh_in.close()

    miRresL = [] ## Empty list to store passed miRNA candidates
    acount  = 0 ## Total entries read
    bcount  = 0 ## Cands found
    for i in fileRead:
        ent = i.strip("\n").split("\t")
        # print("-miR cand ent:",ent)
        putativeName = ent[0]
        candName     = ent[1]
        candLen      = int(ent[7])
        totalVar     = ent[-2]
        # print(putativeName,candName,candLen)
        acount       +=1
        if ent[-1] == "pass":
            if candLen in miRsize:
                # print("-candidate found")
                print("-putative: %s | cand: %s | totalVar: %s" % (putativeName,candName,totalVar))
                miRresL.append((ent))
                bcount+=1

    miRresL_se  = sorted(miRresL,key=lambda x: float(x[13]),reverse=True) ## Sorted on e-val
    miRresL_s   = sorted(miRresL_se,key=lambda x: float(x[14]),reverse=True) ## then Sorted on bitscore
    tagSet      = set() ## An empty set to do uniq later
    miRcandsL   = [] ## Store unique results
    miRcandsD   = {} ## Stores uniq name as key and rest as values
    for i in miRresL_s:
        sname = i[1]
        if sname not in tagSet:
            miRcandsL.append(i)
            miRcandsD[sname] = i
            tagSet.add(sname)
        else:
            pass

    print("-Total miRBLAST entries:%s | cands:%s" % (acount,bcount))
    print("-Total uniq miRBLAST cands:%s" % (len(miRcandsL)))

    return miRcandsL,miRcandsD

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

def mapLibs(fastaFile,fetchedLibs,indexStep):
    
    '''
    This module prepares bowtie maps for all the libraries
    '''

    print("\nFUNCTION: mapLibs")
    print("-Inputs",fastaFile,fetchedLibs,indexStep)
    ## Prepare - Make index, and filenames
    # fastaIndex = "fasta.index"
    fastaIndex = "%s.index" % fastaFile
    
    if indexStep == 1: 
        retcode = subprocess.call(["bowtie-build",fastaFile,fastaIndex])
        if retcode == 0:
            print("+Index file prepared for %s" % (fastaFile))
        else:
            print("-There is some problem in preparing index")
            sys.exit()
    elif indexStep == 0:
        print("Assuming that index of name '%s' exists in this directory" % (fastaIndex))
        pass
    else:
        print("Please choose correct input for setting 'indexStep'")
        sys.exit()

    for alib in fetchedLibs:
        print("+Mapping library:%s" % (alib))
        inFile = '%s.fa' % (alib)
        # print ('+Processing %s for mapping to genome' % (inFile))
        # fastaFile = tagCount2FASTA(inFile,'N') ## Unique reads to FASTA format 
        fastaFile = inFile

        mapFile = ('./%s.map' % (alib))
        print("-",fastaIndex,inFile,fastaFile,mapFile)

        ## Start mapping 
        print ('-Mapping %s processed file to genome' % (alib))
        nproc2 = str(nproc)
        mismat = str(0)

        retcode = subprocess.call(["bowtie","-f","-n",mismat,"-p", nproc2,"-t" ,fastaIndex, fastaFile, mapFile])
        
        if retcode == 0:## The bowtie mapping exit with status 0, all is well
            print('+Bowtie mapping for %s complete' % (inFile) )
        else:
            print ("-There is some problem with mapping of '%s' to cDNA/genomic index - Debug for reason" % (inFile))
            print ("Script exiting.......")
            sys.exit()

    return None

def map2Dict(alib):
    '''
    parse the bowtie map, it assumes that map file have been generated already
    '''

    print("\nFUNCTION: map2Dict")
    mapFile = './%s.map' % (alib)
    fh_in   = open(mapFile,'r')
    mapRead = fh_in.readlines()
    fh_in.close()

    transSet = set() ## To be used as key later
    for i in mapRead:
        ent = i.strip('\n').split("\t")
        # print(ent)
        atrans = ent[2].strip()
        if atrans not in transSet:
            transSet.add(atrans)
        else:
            # print("-transcript %s already added to set once" % (atrans))
            pass


    print("-Preparing dictionaries of sRNAs")
    srnaDict = {}
    for atrans in transSet:
        value = [] ## Store values for key
        # print("-Caching sRNAs for trans:%s" % (atrans))
        for i in mapRead:
            ent = i.strip('\n').split("\t")
            # print(ent)
            phasiname,astrand,trans,apos,phasiseq,trash1,hits,trash2, = ent ## Last column is Comma-separated list of mismatch descriptors. 
            aname               = phasiname
            phasistrand         = str(astrand).translate(str.maketrans("+-","wc"))
            phasilen            = len(phasiseq.strip())
            phasihits           = int(hits)+1 ## This is not a good poxy of sRNA hits, more processing required to get hits, so this number can't be trusted
            phasiflag           = 'S' ## sRNA from bowtie map file
            phasipos            = int(apos)+1 ## Convert 0-based offset to 1-based offset like Pingchuans cluster file
            # print(aname,phasistrand,phasiseq,phasilen,phasihits,phasiflag,phasipos)

            if atrans == trans  and phasistrand == 'w':
                value.append((aname,phasistrand,phasiseq,phasilen,phasihits,phasiflag,phasipos))
            else:
                # sys.exit()
                pass

        # print("Key:%s | values:%s" % (atrans,value))
        srnaDict[atrans] = value

    print("-Total unique transcipts recorded from map file:%s" % (str(len(transSet))))
    print("-Total mappings recorded in dictionary:%s" % (str(len(srnaDict)))) ## This inlcudes multiple mappings on same transcripts


    return srnaDict

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

def predictPrecursor(acand,atransSeq,mapDict,miRcandsD):
    '''Takes a precursor transcipt candidate, fetched FASTA, checks for foldback, 
    checks for location of miR on foldback, apply filters to resturn a falg of status'''

    print("Fetching transcript for candidate precursor:%s" % (acand))
    
    aflag      = "" ## Flag for arm of mapped miRNA
    statusFlag = "" ## Flag for final status of precursor - Y or N 

    ## Foldback and IR coords
    outseq,outinv   = IRchecker(atransSeq,acand)
    IRcoords        = IRparser(outseq,outinv,acand)
    print("Foldback coords:",IRcoords)

    candRes = [] ## List with results using all tags that mapped to this precursor
    if IRcoords:
        start5  = int(IRcoords[0][6])
        end5    = int(IRcoords[0][7])
        start3  = int(IRcoords[0][8])
        end3    = int(IRcoords[0][9])
        ## Is miRsite within foldback i.e. on stem
        avalues = mapDict[acand]
        print("These are miRNA candidates mapped to this precursor",avalues)

        for ent in avalues: ### Multiple tags mapped to same transcript, see how many validates and resturn all
            miRcand = ent[0]
            miRstrand = ent[1]
            miRseq  = ent[2]
            miRlen  = ent[3]
            miRhits = ent[4]
            mappos  = int(ent[6])
            
            if mappos >= start5 and mappos <= end5:
                print("-Position on 5' arm")
                aflag = 5
            elif mappos >= start3 and mappos <= end3:
                print("-position on 3' arm")
                aflag = 3
            else:
                print("-position not in stem")
                aflag = 0 
                pass

            ### Other filters
            if aflag == 3 or aflag == 5:
                statusFlag = "Y"
            else:
                statusFlag = "N"

            ### Add miR candidate wise results
            blastEnt = miRcandsD[miRcand]
            putativeName = blastEnt[0]
            candRes.append((miRcand,putativeName,miRstrand,miRseq,miRlen,miRhits,mappos,acand,atransSeq,len(atransSeq),statusFlag)) ## acand is precursor name and atransSeq is its seqeunce


    else:
        print("-Precursor transcripts doesn't folds")
        pass


    # sys.exit()

    print("-Final status for precursor %s - %s" % (acand,statusFlag))

    return candRes,IRcoords

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
                print("You seem to have input lib_id and chosen wrong libType")
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

#######################################

def main():
    precSumm,precFa,header = inferCandidates(miRBlastF,miRBlastFa)
    ratioComputer(precSumm,precFa,header)


if __name__ == '__main__':
    #### Assign Cores
    if numProc == 0:
        nproc = int(multiprocessing.cpu_count()*0.95)
    else:
        nproc = int(numProc)
    main()
    print("Script finished sucessfully")
    sys.exit()



#### v.01 -> v02
#### miRNA blast results are sorted on e-val, bitscore and uniq on tags