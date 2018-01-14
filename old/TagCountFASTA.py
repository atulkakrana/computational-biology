#!/usr/local/bin/python3


import sys


######### CONFIG ############
tagcount = str(sys.argv[1])###File with IDs
tag_len = (21,22,23,24,25)
############################

def TagtoFASTA(tagcount):
    fh_in = open(tagcount, 'r')
    
    out_file = '%s.fa' % (tagcount.split('.')[0])
    fh_out = open(out_file, 'w')
    
    tag_num = 1 ### for tag naming that will be used in FASTA format
    
    print ('Converting %s to FASTA' % (tagcount) )
    for entry in fh_in:##All the entries of the library
        
        ent =entry.strip('\n').split('\t')
        if len(ent[0]) in tag_len: ## Length of tag specified by the user
            name = ('%s_%s' % (tag_num,ent[1]))
            fh_out.write('>%s\n%s\n' % (name, ent[0]))
            tag_num += 1
        else:
            #print ('Length is not 20nt')
            pass
    print ('The %s has been converted to FASTA file\n' % (tagcount))
    
    fh_out.close()
    
def main():
    TagtoFASTA(tagcount)
    

if __name__ == '__main__':
    main()
    
    sys.exit()


