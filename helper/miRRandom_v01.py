#!/usr/local/bin/python3

### To generate random miRNAs

import os
import glob
import sys
import difflib
import time
import string
import random
import datetime

size = 21
number = 438### Number of random miRs required

def RandGenerate(size,number):
    
    ###Ouput
    year = datetime.date.today().strftime("%Y")
    month = datetime.date.today().strftime("%B")
    date = datetime.date.today().strftime("%d")
    outfile = 'miRRandom_%s_%s_%s_%s.fa' % (size,date,month,number)
    fhout = open(outfile, 'w')
    
    #chars = string.ascii_uppercase
    chars = 'ATGC'
    mirset = set()
    anumber = number+500 ### Adding extra number because few randomly generated miRNAs could be redundant
    for mir in range (anumber):
        rand_miR = ''.join(random.choice(chars) for x in range(size))
        mirset.add(rand_miR)
        if len(mirset) == number:
            print('Number of required random miRNAs generated')
            break
    
    namecount = 1
    for mir in mirset:
        fhout.write('>%s_%s_%s\n%s\n' % (namecount,date,month[:3],mir))
        namecount+=1
    
    fhout.close()
    return outfile 

def main():
    RandGenerate(size,number)
    
    

if __name__ == '__main__':
    main()
    print('\nScript finished sucessfully\n')
    sys.exit()
    
    
    
####
####Get similar complexity as original miRNAs
    
    
