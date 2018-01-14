#!/usr/local/bin/python3

## Written by ATUL to create plots for phased trsncripts and text outputs for R
## In text outputs crick strand positions add added with three to match watson strand positions
import sys,os,operator,time,datetime,string,subprocess,difflib,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

candfile  = "cand.list" ## A single column list
clustfile = "24PHAS_p1e-05_clust.mod.select.txt" ## Typical cluster file
phase     = 24


### Developer Settings

## Plots
minNoiseLen = 18        ## Min noise length included in chart-input, consensus-input and consensusus file, ## Recommended: 18
maxNoiseLen = 28        ## Max noise length included in chart-input, consensus-input and consensusus file, ## Recommended: 24
minPlotLen  = 18        ## Min length considered for plotting
maxPlotLen  = 28        ## Max length used for plotting

def plottrans(candsL,phasedDict):

    '''
    Fetches candidates from list, fetches entries from clustfile, and plots
    '''

    acount = 0
    for aname in candsL:
        aclust      = phasedDict[aname]
        # print("Cluster:",aclust)
        aclust_s    = sorted(aclust, key=operator.itemgetter(7),reverse=False)
        anormclust  = strandmerge(aname,aclust_s)
        # print("\nCluster Sorted:",anormclust)
        chart(aname,anormclust)
        # sys.exit()

    return None

def strandmerge(aname,aclust):
    '''
    Merges abudnaces from same strand
    '''

    print("\nFn: strandmerge ####################")

    posclust        = []
    negclust        = []

    ### Prepare strand-scpecific phasiRNAs list
    ###
    for ent in aclust:
        if ent[1] == '+':
            posclust.append(ent)
        elif ent[1] == '-':
            ## Transform to positive strand for heart-beat chart
            # print("Negative Ent:",ent)
            # newstrand   = '+'
            newpos      = ent[7]+3
            negclust.append((ent[0],ent[1],ent[2],ent[3],ent[4],ent[5],ent[6],newpos))
        else:
            pass

    # print("Negative clust transformed",negclust)
    print("Ent in posclust:%s | Ent in negclust:%s" % (len(posclust),len(negclust)))


    ### Write results for plotting in R - Needs to be merged in excel by user
    ###
    outfile     = "%s.txt" % (aname)
    fh_out      = open(outfile,'w')
    posnegclust = posclust+negclust
    posnegclust_s = sorted(posnegclust, key=operator.itemgetter(7),reverse=False)
    for phasi in posnegclust_s:
        fh_out.write("%s\n" % ('\t'.join(str(i) for i in phasi )))
    fh_out.close()


    ### Merge abudnaces from negative strand to positive strand
    ###
    anormclust  = [] ## Stores a merged cluster
    negset      = set() ## Stores tag fron negative strand foro which abudnace has been added to positive position 
    for aent in posclust:
        aname       = aent[0]
        astrand     = aent[1]
        aabun       = aent[2]
        aseq        = aent[3]
        asize       = aent[4]
        ahits       = aent[5]
        aflag       = aent[6]
        apos        = aent[7]
        ## Find matching entry in negclust based on pistion and merge abundance
        for bent in negclust:
            babun   = bent[2]
            bseq    = bent[3]
            bsize   = bent[4]
            bpos    = bent[7]

            if bpos == apos and asize == bsize:
                # print("\naent:",aent)
                # print("bent:",bent)
                aabun+=babun
                # print("New abun:",aabun)
                negset.add(bseq)
            else:
                pass

        anormclust.append((aname,astrand,aabun,aseq,asize,ahits,aflag,apos))

    ### Add unique positions from negative strand 
    for bent in negclust:
        if bent[3] not in negset:
            anormclust.append(bent)
        else:
            pass

    ### Sort after transforming positions fron negative strand
    anormclust_s    = sorted(anormclust, key=operator.itemgetter(7),reverse=False)

    # print("Original clust",aclust)
    # print("\nMerged clust",anormclust)

    return anormclust_s

def parseList(candsfile):
    '''
    parses single column list
    '''

    print("\nFn: clust2Dict #####################")
    fh_in       = open(candsfile)
    fileread    = fh_in.readlines()
    fh_in.close()

    candsL = set()
    acount = 0
    for i in fileread:
        ent = i.strip('\n').strip()
        candsL.add(ent)
        acount+=1

    print("Candidate entries:%s | Unqiue cached:%s" % (acount,len(candsL)))

    return candsL

def clust2Dict(clustfile):
    '''
    This module reads cluster file and gives back a dictionary of coords. In case of 
    directIRs only one phased cluster, that had valid coords is expected with phasiRNAs
    uniuqe on seq and positions.
    '''

    print("\nFn: clust2Dict #####################")
    # clustFile = 'phased.clust'
    fh_in       = open(clustfile,'r')
    clustRead   = fh_in.read().split('>')
    fh_in.close()

    
    ## Make a set of transcript names
    transSet = set() ## Empty set for phased transcripts
    for i in clustRead[1:]:
        # print(i)
        aclust      = i.split("\n")
        transname   = aclust[0].strip('\n').split('|')[1].split('=')[1].strip().split()[0].strip()
        transSet.add(transname)
        # print(transname)
    print("Trascripts set prepared with %s entries" % (len(transSet)))
    # sys.exit()

    ## Search clusters with transcript names, and append their phasiRNAs
    phasedDict  = {} ## Transcripts wise phasiRNAs
    for gettrans in transSet:
        value       = [] ### List to hold phaseRNAs for a transcript
        print("Caching phasiRNAs for transcript:%s" % (gettrans))
        # phasiset    = set() ## Set of uniq phasiRNAs for this transcript
        
        for ablock in clustRead [1:]:
            aclust      = ablock.split("\n")
            transname   = aclust[0].strip('\n').split('|')[1].split('=')[1].strip().split()[0].strip()
            # print(gettrans,"   -   ",transname)


            if transname == gettrans:
                # print("-Matching block found with transcript:%s" % (transname))
                # print(ablock)
                
                for x in aclust[1:-1]: ## First entry was header and last is always empty
                    ent = x.split("\t")
                    # print("-ent:%s" % (ent))
                    # sys.exit()
                    if ent:
                        # print(ent)
                        phasiname   = ent[1]
                        phasistrand = ent[2]
                        phasipos    = int(ent[3])
                        phasiseq    = ent[5]
                        phasilen    = ent[6]
                        phasiabun   = int(ent[7])
                        phasihits   = int(ent[10].split('=')[1])
                        phasiflag   = 'P' ## Phased from cluster file
                        if phasistrand == "+" or phasistrand == "-": ## Only one strand is used because these are transcripts and not genome, and in case of IR-based, 'c' strand maps to other arm so kind of redundant
                            value.append((phasiname,phasistrand,phasiabun,phasiseq,phasilen,phasihits,phasiflag,phasipos))
                            # phasiset.add(phasiseq)
            else:
                # print("This cluster corresponds to different transcipt %s" % (transname))
                pass

        # print("-%s phasiRNA cached for transcript:%s" % (len(value),gettrans))
        ## Results from multiple libraries will add phasiRNAs, 
        ## Filter most abundant uniq phasiRNAs. and add to dictionary
        # value_s = sorted(value, key=operator.itemgetter(2),reverse=True) ## Sorted on Abundance from high to low
        value_s = sorted(value, key=operator.itemgetter(1),reverse=True) ## Sorted on strand from w to c, we will fetch abundances later
        # print("\n-Snippet of abundance sorted phasiRNA list:",value_s)
        
        checklist   = [] ## List to check if most abundant instance of this phasiRNA have been recorded
        valueuniq   = [] ## List to hold uniq phasiRNAs from multiple libraries, sorted on abundance
        for anent in value_s:
            aphasiseq = anent[3]
            aphasilen = int(anent[4])
            # print("-phasiRNA beind added",anent)
            if aphasiseq not in checklist:
                checklist.append(aphasiseq) ## add to checklist, which means considered once
                valueuniq.append(anent)
            else:
                # print("This phasiRNA %s has been recorded before" % (aphasiseq))
                pass
        print("-%s uniq phasiRNA cached for transcript:%s" % (len(valueuniq),gettrans))

        ## Finally add unique values to dict        
        # print("-trans name:%s | value snippet:%s" % (gettrans,valueuniq[0:3]))
        phasedDict[gettrans] = valueuniq


    print("-phasedDict with %s elements prepared\n" % (str(len(phasedDict))))
    # print("Length of dictionary for 'TRINITY_DN279522_c0_g1_i2' is %s and values:" % (len(phasedDict['TRINITY_DN279522_c0_g1_i2'])))
    # print("And values:",phasedDict['TRINITY_DN279522_c0_g1_i2'])
    # sys.exit()

    return phasedDict

def chart(aname,anormclust):
    '''
    This prepares line plot for pair
    '''

    print("\nFn: Chart")
    # print("Cluster:",anormclust)

    ## Foldback coords
    start1      = int(anormclust[0][7])
    end1        = int(anormclust[-1][7])
    print("PHAS Start:%s | PHAS end:%s" % (start1,end1))
    # sys.exit()
    
    ## Prepare #####
    plotFile        = "%s" % (anormclust[0][0]) 
    # anormclust_s    = sorted(anormclust,key=lambda x: float(x[7]),reverse=False)
    # bnormclust_s    = sorted(bnormclust,key=lambda x: float(x[7]),reverse=False)
    ### Prepare abundance and positions list for A
    #############################################
    print("Prepare list for 'A' for plotting")
    a21x        = [] ## Record size specific positions
    a22x        = []
    a23x        = []
    a24x        = []
    aOthersx    = [] ## Other sizes in allowed in plotLen
    aALLx       = [] ## Collect psoitions for all sizes (sorted)
    a21y        = [] ## Record size specific abundances
    a22y        = []
    a23y        = []
    a24y        = []
    aOthersy    = [] ## Other sizes in allowed in plotLen
    aALLy       = [] ## Collect abundances for all sizes (sorted on position)
    aphasposS   = set() ### Records all phased psoition to define x-axis
    for aphasi in anormclust:
        # print("Ent",aphasi) ### name,strand,abun,seq,size,hits,flag,pos
        phasipos    = int(aphasi[7])
        phasiabun   = round(math.log(int(aphasi[2]),2),2)
        phasilen    = int(aphasi[4])
        phasiflag   = aphasi[6]
        # print("Phasi Abundance:%s" % (phasiabun))
        if phasiflag == "P": ## Record the psoition to define x-axis
            aphasposS.add(phasipos)

        aALLx.append(phasipos)
        aALLy.append(phasiabun)

        if phasilen     == 21:
            a21x.append(phasipos)
            a21y.append(phasiabun)
        elif phasilen   == 22:
            a22x.append(phasipos)
            a22y.append(phasiabun)
        elif phasilen   == 23:
            a23x.append(phasipos)
            a23y.append(phasiabun)
        elif phasilen   == 24:
            a24x.append(phasipos)
            a24y.append(phasiabun)
        elif phasilen >= minPlotLen and phasilen <= maxPlotLen:
            aOthersx.append(phasipos)
            aOthersy.append(phasiabun)
        else:
            print("The sRNA length not allowed for plotting - will be trashed")
            pass
    print("Length of 21 list:%s | 22 list:%s | 23 list:%s | 24 list:%s" % (len(a21x),len(a22y),len(a23x),len(a24x)))
    # print("phased positions on A:",aphasposS)


    #### Prepare chart ###############
    ##################################

    a21 = plt.scatter(a21x, a21y, s=[n**2 for n in a21y], c='#27408B', lw = 0, alpha=0.5)
    a22 = plt.scatter(a22x, a22y, s=[n**2 for n in a22y], c='#2E8B57', lw = 0, alpha=0.5)
    a23 = plt.scatter(a23x, a23y, s=[n**2 for n in a23y], c='#7D26CD', lw = 0, alpha=0.5)
    a24 = plt.scatter(a24x, a24y, s=[n**2 for n in a24y], c='#EE9A00', lw = 0, alpha=0.5)
    aOthers = plt.scatter(aOthersx, aOthersy, s=[n**2 for n in aOthersy], c='#383838', lw = 0, alpha=0.5)
    aAll = plt.plot(aALLx,aALLy, c='#828282', lw = 0.5)

    ##### Add phasing lines #############
    aphasmin        = min(aphasposS)
    aphasmax        = max(aphasposS)
    aphasx_marker   = list(np.arange(aphasmin,aphasmax+phase,step=phase)) ## Phased position from first till last phase, plus an extra marker added
    aphasy_marker   = [max(aALLy)]*len(aphasx_marker) ## Most abundant value on y-axis multipled the times of x marker to get list of same size as aphasx_marker

    for i,j in zip(aphasx_marker,aphasy_marker):
        # plt.axvline(i,linewidth=0.5, color='#454545')
        plt.plot((i,i),(0,j),linewidth=0.3, color='#454545',linestyle='--') ## plt.plot((x1, x2), (y1, y2), 'k-')

    atitle      = "%s" % (anormclust) 
    plt.axhline(linewidth=0.5, color='#454545')
    plt.ylabel('Abundance', fontproperties=font_manager.FontProperties(size=8))
    plt.xlabel('Fold-back normalized positions [phas-lines computed on first phasiRNA - for indicative purpose only]', fontproperties=font_manager.FontProperties(size=7))
    plt.title("First phase: %s(5') | Last phase: %s(5') | Total phase: %s(5')" % (aphasmin,aphasmax,len(aphasx_marker)-1),fontproperties=font_manager.FontProperties(size=7)) ## In figures 1 or 2 phases might be extra, these are noise sRNAs showing up in phass

    plt.savefig(plotFile, format=None , facecolor='w', edgecolor='w', orientation='portrait', papertype=None,  transparent=False, bbox_inches=None, pad_inches=0)
    plt.clf() ## Clear figure so that next plot is different file

    print("\nFn: Exiting Chart")

    return None

def main():
    candsL      = parseList(candfile)
    phasedDict  = clust2Dict(clustfile)
    plottrans(candsL,phasedDict)


if __name__ == '__main__':
    main()
    print("Done!!")
    sys.exit()

