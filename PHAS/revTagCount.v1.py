#!/usr/local/bin/python3

import os,sys

fileEnd = "chopped.txt"

def revTagCount(fileEnd):
    fileList = [file for file in os.listdir('./') if file.endswith (('%s' % fileEnd))]
    print("These are all the files:",fileList)

    for afile in fileList:

        print("\nReverse complementing file:%s" % (afile))
        fh_in = open(afile, 'r')
        aread = fh_in.readlines()

        outfile = "%s.revcomp.txt" % afile.rpartition(".")[0]
        fh_out = open(outfile,'w')

        for i in aread:
            tag,count = i.strip("\n").split("\t")
            # tagrev = tag.strip()[::-1]
            tagrevcomp = tag.strip()[::-1].translate(str.maketrans("TAGC","ATCG"))
            fh_out.write("%s\t%s\n" % (tagrevcomp,count))

        fh_out.close()
        fh_in.close()

def main():
    revTagCount(fileEnd)

if __name__ == '__main__':
    print("\n\nAyush- You could have written this script!!!\n\n")
    main()
    print("\nJob done succesfully - One beer to atul !!")
    sys.exit()