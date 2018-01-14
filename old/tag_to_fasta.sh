#!/bin/bash
#
#
clear
#
#PLEASE CHANGE INPUT FILENAME 'tag_file' before proceeding
#
echo -e "Today is $(date)"
echo -e "\n"
echo -e "Convert tag count file to FASTA format\n"
echo -e "The tag counts will be named as the sequence only\n"
#echo "Please change the input filename 'tag_file' and output filename 'result_file"
#
echo -e "Processing time depends upon the size of your input file\n"
echo -e "Please wait......\n"
#

$file='tomato_pare.txt'

#Select the first column from tag_count file
#PLEASE RENAME INPUT FILE DOWN HERE
awk '{print $1}' $file > file1
#
#copy each line of file after same line separted with tab
#
paste file1 file1 > file2
#
#insert > before every line
sed 's/^/>/' file2 > file3
#
#substitute tab with new line
#PLEASE RENAME  OUTPUT FILE DOWN HERE
sed 's/\t/\n/' file3 > $file.fa
#
#
echo " Today is $(date)"
#Cleaning intermediate files
rm file1 file2 file3
#
#print signal 'DONE'
echo -e "Success\n"
exit 0




