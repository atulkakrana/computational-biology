#usr/bin/perl
#usage: perl extract_matched_header.pl contigs.fa test3 > test3.fa

my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];

open(FILE1, "<$inputfile1");

my %data = ();
my $sub;
my $text = "";

while(<FILE1>)
{
	chomp;
	if(m/^>/)
	{
		if($text ne "")
		{
			$data{$sub} = $text;
		}
		$sub = substr($_, 1);
		$text = "";
	}
	else
	{
		$text = $text . $_;
	}
}
if($text ne "")
{
	$data{$sub} = $text;
}
close(FILE1);

open(FILE2, "<$inputfile2");

my @data2 = ();
while(<FILE2>)
{
	chomp;
	$sub = substr($_, 1);
	push(@data2, $sub);
}
close(FILE2);

foreach my $k (keys %data)
{
	if((grep $_ eq $k, @data2))
	{
		print ">$k\n$data{$k}\n";
	}
}
