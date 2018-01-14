## This script reads an expression matirx or final combined file with counts, fpkm and anno and expands 
## entries like PB1.1,PB1.2,PB1.3 which got clubbed by cuff-norm

## Example: XLOC_000001  -   -   XLOC_000001 PB.1,PB.3   TSS1,TSS2,TSS3,TSS4,TSS5
## In above entry coulmn 5 needs to be expanded making seprate entries for PB.1 and PB.3

import sys,os

afile   = "anno_fpkm_counts.v2.genes.txt"           ## File 
acol    = 5                                         ## The column which needs to expanded, i.e. written into multiple lines - Excel format
sep     = "\t"                                      ## delimiter used in file



def expander(afile):

    fh_in   = open(afile,'r')
    header  = fh_in.readline()
    aread   = fh_in.readlines()

    outfile = '%s_expanded.txt' % (afile.rpartition('.')[0])
    fh_out  = open(outfile,'w')
    fh_out.write('%s\n' % header.strip('\n'))

    for line in aread:
        ent         = line.strip('\n').split(sep)
        # print("\nEntry",ent)
        fusedvals   = ent[acol-1]
        ntrans      = fusedvals.split(',')
        # print("nTrans",ntrans)
        if len(ntrans) > 1:
            print("This entry corresponds to multiple trans:",ntrans)
            for i in ntrans:
                # print('\t'.join(x for x in ent[:acol-1]),i,'\t'.join(x for x in ent[acol:]))
                fh_out.write('%s\t%s\t%s\n' % (('\t'.join(x for x in ent[:acol-1])),i,('\t'.join(x for x in ent[acol:]))))
                # sys.exit()

        else:
            fh_out.write('%s\n' % ('\t'.join(x for x in ent)))
            pass

    fh_in.close()
    fh_out.close()

    return outfile


def main():
    outfile = expander(afile)


if __name__ == '__main__':
    main()
