#!/usr/local/bin/python3

##This module cleans the FASTA file header
##Trunchates the header till first whitespace and pipe, you can change the nature below

##Usage: python3 clean_fasta_header_v1.py ENTER_FASTA_FILE

import re,sys

## Config ####
inp_file_name = str(sys.argv[1])
mode  = int(sys.argv[2])
# mode = 0           ## 0: Clean Fasta as is 1: clean and reverse 2: clean and rev comp | 3: Sequnce complemented only | 4: Convert DNA to RNA (T->U) | 5: Convert RNA to DNA (U->T)



def FASTAClean(filename,mode):
    
    '''Cleans FASTA file - multi-line fasta to single line, header clean, empty lines removal'''

    ## Read seqeunce file
    fh_in       = open(filename, 'r')
    print ('\nCleaning "%s" FASTA file' % (filename))
    
    ## Write file
    if mode == 0:
        out_file1 = ('%s.clean.fa' %            (filename.rpartition('.')[0]))
    elif mode == 1:
        out_file1 = ('%s.clean.rev.fa' %        (filename.rpartition('.')[0]))
    elif mode == 2:
        out_file1 = ('%s.clean.revcomp.fa' %    (filename.rpartition('.')[0]))
    elif mode == 3:
        out_file1 = ('%s.clean.comp.fa' %       (filename.rpartition('.')[0]))
    elif mode == 4:
        out_file1 = ('%s.clean.rna.fa' %        (filename.rpartition('.')[0]))
    elif mode == 5:
        out_file1 = ('%s.clean.dna.fa' %        (filename.rpartition('.')[0]))
    else:
        print("Input correct mode- 0: Normal | 1: Seqeunces reversed | 2: Seqeunces reverse complemented | 3: Seqeunces complemented only")
        print("USAGE: cleanFasta.v.x.x.py FASTAFILE MODE")
        sys.exit()

    fh_out1     = open(out_file1, 'w')

    out_file2   = ('%s.summ.txt' % (filename.split('.')[0]))
    fh_out2     = open(out_file2, 'w')
    fh_out2.write("Name\tLen\n")
    
    print ('\nCleaning "%s" FASTA file\n' % (filename))
    
    fasta = fh_in.read()
    fasta_splt = fasta.split('>')
    acount = 0 ## count the number of entries
    empty_count = 0
    for i in fasta_splt[1:]:
        ent = i.split('\n')
        name = ent[0].split()[0].strip()
        seq = ''.join(x.strip() for x in ent[1:]) ## Sequence in multiple lines
        alen    = len(seq)

        if mode == 0:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 1:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq[::-1]))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 2:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq[::-1].translate(str.maketrans("TAGC","ATCG"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 3:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq.translate(str.maketrans("TAGC","ATCG"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 4:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq.translate(str.maketrans("TAGC","UAGC"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        elif mode == 5:
            if seq:
                fh_out1.write('>%s\n%s\n' % (name,seq.translate(str.maketrans("UAGC","TAGC"))))
                fh_out2.write('%s\t%s\n' % (name,alen))
                acount+=1
            else:
                empty_count+=1
                pass
        else:
            print("Please enter correct mode")
            pass


        acount+=1
    
    fh_in.close()
    fh_out1.close()
    fh_out2.close() 

    print('\nThe fasta file with reduced header: "%s" with total entries %s has been prepared\n' % (out_file1, acount))
    print('**There were %s entries found with empty sequences and were removed**' % (empty_count))
    return None

def fastatolist(alib):
    '''
    New FASTA reader
    '''

    ### Sanity check
    try:
        f = file(alib)
    except IOError:                    
        print ("The file, %s, does not exist" % (alib))
        return None

    order       = []
    sequences   = {}

    for line in f:
        if line.startswith('>'):
          name = line[1:].rstrip('\n')
          name = name.replace('_', ' ')
          order.append(name)
          sequences[name] = ''
        else:
          sequences[name] += line.rstrip('\n').rstrip('*')
    
    print ("%d sequences found" % (len(order)))

    return order, sequences

def main():
    FASTAClean(inp_file_name,mode)

if __name__ == '__main__':

    main()
    sys.exit()


'''
Created on Apr 23, 2012

@author: atul
'''
## v1 -> v1.1
##Added function to remove gene entries with empty sequences
##Changed the output file name by removing '.fa' in between of name

## v1.1 -> v1.2
## Removes space at the end of multiline fasta

##v1.2 -> v1.3
## Name or output file changed
## Modes added fro reverse and reverse compelmenting FASTA file 

## v1.3 -> v1.6
## Added two more modes to convert DNA to RNA and RNA to DNA

## v1.6 -> v1.7
## Fixed a file renaming bug

