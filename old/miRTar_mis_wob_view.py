#!/usr/local/bin/python3

###Find mismatches in 10 and 11th position

mis_pos = 11 ###look for mismatches at these positions, keep same if there is one position
mis_pos2 = 12 ###look for mismatches at these positions, keep same if there is one position

def MisPos(afile):
    fh_in = open(afile, 'r')
    mismatch_file = 'mismatch_%s_%s.csv' % (mis_pos,mis_pos2)
    fh_out = open(mismatch_file, 'w')
    wobble_file = 'wobble_%s_%s.csv' % (mis_pos,mis_pos2)
    fh_out2 = open(wobble_file, 'w')
    for i in fh_in: ###For every result entry
        print ('\n%s' % (i.strip('\n')))
        anent = i.split(',')
        tar = anent[6].strip()###There is a whitespace in file just before target start
        mir_r = anent[5][::-1]###miR reversed to read from left to right (5-3') for simplicity
        tar_r = tar[::-1]
        tar_rc = tar.translate(str.maketrans("AUGC","UACG"))[::-1]###translated and reversed for simplicity (3-5')
        
        print(mir_r,tar_rc)
        nt_count = 1 ###keep track of actual position, strats from position 1 in mir and tar and not '0' as in python
        pos_count = 0
        
        for bp in range (len( mir_r)):
            print (mir_r[bp], tar_rc[bp])
            if mir_r[bp] == '-' or tar_rc[bp] == '-':##Skip gap, don't count in in nt_count if there is a gap or bulge in miR
                #print(tar[::-1][aling_count])
                print ('Gap/Bulge found at %s pos: %s %s' % (nt_count,mir_r[bp],tar_rc[bp]))
                pos_count += 1
                
                pass
            else: ####If no gap in miRNA               
                #print(mir_r[bp],tar_rc[bp])
                if nt_count == mis_pos or nt_count == mis_pos2:###when at specified position
                    if str(mir_r[bp]) == 'G' and str(tar_r[bp]) == 'U' or str(mir_r[bp]) == 'U' and str(tar_r[bp]) == 'G': ##Wobble testing
                        print ("Wobble at %sth: %s - %s" % (nt_count,mir_r[bp],tar_r[bp]))
                        fh_out2.write(i)
                        
                    elif str(mir_r[bp]) is not str(tar_rc[bp]): ### Mismatch testing
                        print ("**mismatch at %sth: %s - %s**" % (nt_count,mir_r[bp],tar_r[bp]))
                        fh_out.write(i)
                nt_count += 1#Count if no gap
                pos_count+=1#Counted in all cases
    fh_out.close()
    fh_out2.close()
    return mismatch_file, wobble_file


def main():
    res_file = MisPos(afile)
    

if __name__ == '__main__':
    afile = 'scoring_inp_22'
    main()
    