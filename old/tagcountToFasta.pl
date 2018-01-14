#! /usr/bin/perl -w

my $inputfile = "ve_tag_count.txt";

my $outputfile = "ve.fasta";
open(FILE, "<$inputfile");
open(OUT, ">$outputfile");

while(<FILE>)
{
	chomp;
	my @data = split(/\t/, $_);
	print OUT ">$data[0]\n$data[0]\n";
}
close(FILE);
close(OUT);
