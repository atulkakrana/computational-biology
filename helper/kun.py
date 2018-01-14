#!/usr/local/bin/python3


import sys

inp_file = str(sys.argv[1])###File with IDs
seq_file = str(sys.argv[2])###File from where entries need to ne extracted, should have clean header


#inp_file = './names'
#seq_file = './hairpin.fa'
match_file = './matched_entries'

def MatchHead():

    fh_in = open(inp_file, 'r')
    fh_in2 = open(seq_file, 'r')
    fh_out = open(match_file, 'w')##Result file
    main_file = fh_in2.read().split('>')[1:]
    adict  ={}
    for i in main_file:
        ent = i.split('\n')
        seq = ''.join(ent[1:])##Sequence in multiple lines
        name = ent[0].strip().lower().split()##Tuple as key holds acomplete header as list
        #print(name)
        atuple = (seq,name[1:])##Hold sequence and rest of header
        #print (name[0])
        adict[name[0]] = atuple

    for i in fh_in:
        name = i.strip('\n').split('-')
        #print (name)
        new_name = '%s-%s' % (name[0],name[1])
        #print (new_name)
        
        ##Fetch dfrom dictionary
        try:
            seq = adict[new_name.lower()][0]
            header_list = adict[new_name.lower()][1:][0]
            header = ' '.join(header_list)
            #print (header)
            #print(header[0][0])
            fh_out.write('>%s %s\n%s\n' % (new_name, header,seq))
            
        except KeyError:
            pass
        
    fh_in.close()
    fh_in2.close()
    fh_out.close()
    
    return match_file
        


def main():
    matched = MatchHead()
    
        
if __name__ == '__main__':
    main()
    print ('\n Job Done..Cheers !_!')

