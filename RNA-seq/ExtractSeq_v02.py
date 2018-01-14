#!/usr/bin/python3

'''Script to parse cuffmerge output and extrct sequences from the genome used 
for the experiment - Written by atulkakra@udel.edu'''

import os,sys,operator
import itertools as it


### User settings
coordFile   = 'merged.gtf' ## Mandatory Input 
transFile   = 'isoforms.attr_table' ## Required for extracting transcript
genoFile    = '/alldata/Genomic/Asparagus/UGA/v2/GenomicSEQ/Asparagus.V2.0.genome.stable_modifed.fa'
mode        = 1 ## 1: Extract CDS [use merged GTF] | 2: Extract Gene transcript [use genes attribute file] | 3: Extract Isoforms [use isofor attribute file] 4: Extract Promoters [use gene attribute file]

def parseCoordFile(coordFile):

	'''Parse Rocket merged GTF file and prepare
	a list of coords'''

	print("\nModule:parseCoordFile")
	print("Parsing GTF file and preparing coords dictionary\n")

	coordsList = [] ## List to store coords for sequence extraction
	coordsDict = {} ## Dictionary with unique gene_id ans key and other info as values
	transList = [] ##List of transcript coords - Generated only in mode 2

	fh_in = open(coordFile,'r')
	
	if (mode == 1) or (mode == 2):
		fileRead = fh_in.readlines()
		faultCount = 0

		for ent in fileRead:
			# print(ent)
			chr_id,trash1,trash2,start,stop,dot1,strand,dot2,info = ent.strip('\n').split('\t')
			# print(chr_id,start,stop,strand,info)
			info_splt = info.split(';')[:4]
			gene_id = info_splt[0].split('"')[1].replace('"','')
			transcript_id = info_splt[1].split('"')[1].replace('"','')
			exon = info_splt[2].split('"')[1].replace('"','')
			gene_name = info_splt[3].split('"')[1].replace('"','')
			# print(chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name)
			
			## Change strand format
			if strand == '+':
				strand = 'w'
			elif strand == '-':
				strand = 'c'
			else:
				# print('Wrong strand encountered while splitting: %s' % strand)
				# print(ent)
				faultCount += 1

			## Record
			coordsList.append((chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name))
			coordsDict[gene_id] = ((chr_id,strand,gene_name))
	## Extract entries and add gene strand from coordsList
	
	if mode == 2:
		print("Preparing coords for transcript extraction")
		fh_in = open(transFile,'r')
		fh_in.readline() ## Remove header
		readFile = fh_in.readlines()

		for line in readFile:
			ent = line.strip('\n').split('\t')
			# print(ent)
			gene_id = ent[0]
			gene_name = ent[4]
			coords = ent[6]
			# print(coords)
			chr_id = coords.split(':')[0]
			start2,stop2 = coords.split(':')[1].split('-')

			## Get strand from the coordsDict
			chr_id,strand,gene_name = coordsDict[gene_id]
			transList.append((chr_id,start2,stop2,strand,gene_id,gene_name))



	##Sort the list on chromosome,strand and stop
	sorted(coordsList, key=operator.itemgetter(0,3,1))
	sorted(transList, key=operator.itemgetter(0,3,1))

	# print(coordsList[:3])
	print('Total entries recorded:%s | Entries with problem:%s\n'% (len(coordsList),faultCount))
	fh_in.close()

	return coordsList,transList

def makeGenomeDict(genoFile):

	''' Prepare a dictionary of genome with chr/scaffold as key
	and seq as value'''

	print("\nModule:makeGenomeDict")
	print("Preparing genome dictionary\n")

	fh_in = open(genoFile,'r')
	readFile = fh_in.read().split('>')
	genomeDict = {} 

	for i in readFile[1:]:
		# print(i.split('\n'))
		chromo = i.split('\n')
		head = chromo[0].split(' ')[0] ## Shorten the header to be used as key
		chrSeq = '' ## Empty string for sequence
		for i in chromo[1:]:
			# print("seq:%s" % (i))
			chrSeq+=i
		# print(">%s\n%s" % (head,chrSeq))
		
		## Add to dict
		genomeDict[head] = chrSeq

	print("Entries in the genomeDict:%s\n" % (len(genomeDict)))
	return genomeDict

def extractCDS(coordsList,genomeDict):

	''' Module extracts CDS sequences from the merged gtf fileby stiching togther the exons for every transcript
	This includes isoform for same genes - I am lazy, could have used mySQLite'''

	print("\nModule:extractCDS")
	print("Constructing CDS from Exons")

	featureSeq = '' ## Intialize empty sequence and then empty after every feature (below)
	seqList = [] ## List to store results
	chrLoaded = '' ## Will keep track of running chromosome, since coordsList is sorted this will reduce runtime
	exonCount  = 0 ## Keep counts of exons for a gene
	
	# for ents in zip(coordsList[0:],coordsList[1:]):
	for i in range(0,int(len(coordsList))-1):
		ent1 = coordsList[i]
		ent2 = coordsList[i+1]
		chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name = ent1
		print("This is the ent:",chr_id,start,stop,strand,gene_id,transcript_id,exon,gene_name)
		## To accomodate missing one nucleotide at start
		if start == 0:
			pass ## Bug-1 fixed
		else:
			start = int(start)-1

		## Load Chromosome 
		if chrLoaded:
			## This is not the first entry and chromosome has been loaded
			if chr_id == chrLoaded:
				##No need to fetch chromosme
				pass
			else:
				## The chromosome is different from last entry - Load and update chrLoaded
				chromo = genomeDict[ent1[0]] ## Use chr_id to fetch sequence
				chrLoaded = chr_id ##update chrLoaded
		else:
			## This is the first entry - Load chromosme and update chrLoaded
			chromo = genomeDict[ent1[0]] ## Use chr_id to fetch sequence
			chrLoaded = chr_id ## update chrLoaded

		## Fetch Sequences
		if ent2:
			## Not EOF
			if transcript_id == ent2[5]:
				## Exons belong to same isoform
				# chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
				exonSeq = chromo[int(start):int(stop)]
				# print(exonSeq)
				featureSeq+=exonSeq
				exonCount +=1
			
			else:
				if exonCount >= 1:
				## This is last exon of multi exon gene
					# chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
					exonSeq = chromo[int(start):int(stop)]
					print("This is final exon of multi-exon transcript:%s" % (exonSeq))
					featureSeq+=exonSeq
					# print("This is final transcript:%s" % (featureSeq))
					if strand == 'w':
						seqList.append((gene_id,transcript_id,featureSeq))
					elif strand == 'c':
						seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
					else:
						## Strand is uknown and must be '.' in file - So record both
						print("GeneID:%s - Unknown strand" % (gene_id))
						seqList.append((gene_id,transcript_id,featureSeq))
						seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
					## Current transcript is complete - reintialize transcript-specific variables
					featureSeq = '' ## Empty for next gene as current transcript is complete
					exonCount = 0 ## Set end of gene

				else:
					## This is the first Exon of single exon transcript
					featureSeq = '' ## Just to make sure
					# chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
					exonSeq = chromo[int(start):int(stop)]
					print("This is single exon transcript:%s\%s" % (ent1,exonSeq))
					featureSeq += exonSeq
					if strand == 'w':
						seqList.append((gene_id,transcript_id,featureSeq))
					elif strand == 'c':
						seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
					else:
						## Strand is uknown and must be '.' in file - So record both
						print("GeneID:%s - Unknown strand" % (gene_id))
						seqList.append((gene_id,transcript_id,featureSeq))
						seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))

					## Current transcript is complete - reintialize transcript-specific variables
					featureSeq = '' ## Empty for next gene as current transcript is complete
					exonCount = 0 ## Set end of gene
		else:
			## This is EOF - Final entry - Wrap up the transcript
			print ("This is the final entry")
			if exonCount >= 1:
			## This is last exon of multi exon gene
				# chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
				exonSeq = chromo[int(start):int(stop)]
				print("This is final exon of multi-exon transcript:%s" % (exonSeq))
				featureSeq+=exonSeq
				# print("This is final transcript:%s" % (featureSeq))
				if strand == 'w':
					seqList.append((gene_id,transcript_id,featureSeq))
				elif strand == 'c':
					seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
				else:
					## Strand is uknown and must be '.' in file - So record both
					print("GeneID:%s - Unknown strand" % (gene_id))
					seqList.append((gene_id,transcript_id,featureSeq))
					seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))


				## Current transcript is complete - reintialize transcript-specific variables
				featureSeq = '' ## Empty for next gene as current transcript is complete
				exonCount = 0 ## Set end of gene
				break

			else:
				## This is the first Exon of single exon transcript
				featureSeq = '' ## Just to make sure
				# chromo = genomeDict[ent1[0]] ## Use chr to fetch sequence
				exonSeq = chromo[int(start):int(stop)]
				print("This is single exon transcript:%s\%s" % (ent1,exonSeq))
				featureSeq += exonSeq
				if strand == 'w':
					seqList.append((gene_id,transcript_id,featureSeq))
				elif strand == 'c':
					seqList.append((gene_id,transcript_id,featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))
				else:
					## Strand is uknown and must be '.' in file - So record both
					print("GeneID:%s - Unknown strand" % (gene_id))
					seqList.append((gene_id,transcript_id,featureSeq))
					seqList.append((gene_id+'_rc',transcript_id+'_rc',featureSeq[::-1].translate(str.maketrans("TACG","ATGC"))))

				## Current transcript is complete - reintialize transcript-specific variables
				featureSeq = '' ## Empty for next gene as current transcript is complete
				exonCount = 0 ## Set end of gene
				break
	
	print("Total CDS extracted: %s" % (len(seqList)))
	return seqList

def extractTrans(genomeDict,transList):
	''' Module takes gene attributes file from Rocket
	to extract transcripts '''

	print("\nModule:extractTrans")
	print("Constructing CDS from Exons")
	
	chrLoaded = '' ## Will keep track of running chromosome, since coordsList is sorted this will reduce runtime
	transSeqList = [] ## List to record results

	for i in transList:
		chr_id,start,stop,strand,gene_id,gene_name = i
		print("This is the ent:",i)
		## To accomodate missing one nucleotide at start
		if start == '0':
			# print("Start is Zero")
			start = int(start) ## Bug-1 fixed
		else:
			start = int(start)-1

		## Load Chromosome 
		if chrLoaded:
			## This is not the first entry and chromosome has been loaded
			if chr_id == chrLoaded:
				##No need to fetch chromosme
				pass
			else:
				## The chromosome is different from last entry - Load and update chrLoaded
				chromo = genomeDict[chr_id] ## Use chr_id to fetch sequence
				chrLoaded = chr_id ##update chrLoaded
		else:
			## This is the first entry - Load chromosme and update chrLoaded
			chromo = genomeDict[chr_id] ## Use chr_id to fetch sequence
			chrLoaded = chr_id ## update chrLoaded

		# print("This is the chromosome:",chr_id,chromo)

		## Extract Transcript - Both strand in case strand is unknown
		if strand == 'w':
			featureSeq = chromo[int(start):int(stop)]
			# print("FeatureSeq:",featureSeq)
			transSeqList.append((gene_id,gene_name,featureSeq))
			
		elif strand == 'c':
			featureSeq_rc = chromo[int(start):int(stop)][::-1].translate(str.maketrans("TAGC","ATGC"))
			# print("FeatureSeq:",featureSeq_rc)
			transSeqList.append((gene_id,gene_name,featureSeq_rc))
			
		elif strand == '.':
			## Strand is uknown and must be '.' in file - So record both
			print("GeneID:%s - Unknown strand" % (gene_id))
			featureSeq = chromo[int(start):int(stop)]
			featureSeq_rc = chromo[int(start):int(stop)][::-1].translate(str.maketrans("TAGC","ATGC"))
			# print("FeatureSeq:",featureSeq)
			# print("FeatureSeq:",featureSeq_rc)
			
			transSeqList.append((gene_id,gene_name,featureSeq))
			transSeqList.append((gene_id+'_rc',gene_name+'_rc',featureSeq_rc))

	return transSeqList

def extractProm(coords,genoFile):
	pass
def writer(resList):
	'''Write the CDS, transcripts or promoter identified above
	resList has three entries - ID1,ID2 and Seq'''

	print("\nModule:writer")
	print("Writing CDS, transcripts or promoter")

	if mode == 1:
		outFile = '%s.fa' % (coordFile.rpartition('.')[0])
		fh_out = open(outFile,'w')

		for i in resList:
			# print(i)
			gene_id,transcript_id,featureSeq = i
			fh_out.write('>%s_%s\n%s\n' % (gene_id,transcript_id,featureSeq))
		pass

	elif mode == 2:
		outFile = '%s.fa' % (transFile)
		fh_out = open(outFile,'w')

		for i in resList:
			# print(i)
			gene_id,gene_name,featureSeq = i
			fh_out.write('>%s_%s\n%s\n' % (gene_id,gene_name,featureSeq))
		pass

def main():
	coordsList,transList = parseCoordFile(coordFile)
	genomeDict = makeGenomeDict(genoFile)
	if mode == 1:
		resList = extractCDS(coordsList,genomeDict)
	elif mode == 2:
		resList = extractTrans(genomeDict,transList)
	else:
		print("Please choose correct mode\n")
		sys.exit()
	
	writer(resList)

if __name__ == '__main__':
	main()
	sys.exit()


###### LOG ########
## v01
## One last entry is missed in main for loop of extractCDS
## Fixed - Bug-1 - Some coordinates had 0 as start and subtracting 1 was giving empty seqiences

## v01 -> V02
## Fixed major bug in reverse complement C and G were not being translated
