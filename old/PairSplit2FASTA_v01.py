#!/usr/local/bin/python3
##Usage: ./PairSplitFASTA_v01.py -l PAIRED_END_LIB	-b POSITION_TO_SPLIT
##Written by kakrana@gmail.com

import os
import sys
import re
import argparse

########Parameters################

# lib = 'test.fasta'
# breakpt = '75'

####################################################################################################################



def SplitPaired(lib,breakpt):

	breakpt = int(breakpt)##make sure its int

	fh_in = open(lib,'r')
	fasta_file = lib.split('.')[0] + '.fa'
	fh_out = open(fasta_file, 'w')

	regex = re.compile('^@[A-Za-z]')
	blacklist = re.compile('[$%"/&+#]')
	print ("\n***Converting paired end '%s' file to Single end FASTA file '%s'\n" % (lib,fasta_file))
	for i in fh_in:
		if re.search(regex,i):##Possible header
			if re.search(blacklist,i): ## Possibly contaminated header
				# print (i)
				pass

			else: ##genuine header
				header = i.strip('\n').split()[0][1:]
				seq = fh_in.readline().strip('\n')##Simplest possible implementation, use 'iterzip' in future
		
				read5p = seq[:breakpt]
				read3p = seq[breakpt:]

				# print (">%s_5p\n%s\n>%s_3p\n%s\n" % (header,read5p,header,read3p))
				fh_out.write(">%s_5p\n%s\n>%s_3p\n%s\n" % (header,read5p,header,read3p))

	
	fh_in.close()
	fh_out.close()

	return fasta_file




def main(lib,breakpt):
	fasta_file = SplitPaired(lib,breakpt)


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	##Command line options
	parser.add_argument("-l", "--lib", dest = "lib", help="Paired end FASTQ file")
	parser.add_argument("-b", "--break", dest = "breakpt", help="Position to split")

	args = parser.parse_args()

	print( "lib {} breakpt {}".format(
	        args.lib,
	        args.breakpt
	        ))

	main(args.lib,args.breakpt)
	print('\n**Script finished sucessfully**\n')
	sys.exit()

############AA####################TT###############################UU###########################LL###################