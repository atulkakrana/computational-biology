#!/usr/local/bin/python3


import sys

seq_file = str(sys.argv[1])###File with positive strand
rev_comp = '%s_RevComp.fa' % (seq_file.split('.')[0] )

def RevComp(seq_file,rev_comp):

    fh_in = open(seq_file, 'r')
    fh_out = open(rev_comp, 'w')##Result file
    main_file = fh_in.read().split('>')[1:]
    print('Preparing reverse complement of FASTA file')
    for i in main_file:
        ent = i.split('\n')
        seq = ''.join(ent[1:])##Sequence in multiple lines
        new_name = '%s_reverse' % (ent[0])
        #seq_r = seq[::-1]
        seq_rc = seq.translate(str.maketrans("atgcATGC","tacgTACG"))[::-1]

        fh_out.write('>%s\n%s\n' % (new_name, seq_rc))
    fh_out.close()
    fh_in.close()
    return fh_out


def main():
    RevComp(seq_file,rev_comp)



if __name__ == '__main__':
    main()
    print('The revese complement of FASTA file is ready')
    sys.exit()
    



    