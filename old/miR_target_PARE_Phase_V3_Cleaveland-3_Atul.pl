#!/usr/bin/perl
# written by Jixian Zhai, zhai@dbi.udel.edu
use DBI;

#input required

#1. list of phased loci to extract all phased siRNAs
#2. PARE DB (DB name)
#3. small RNA summary2 table
#4. Target sequence (could be whole genome, in fasta format)

############################################# General Input #################################################
#tasi_loci_sequence (subject)
$seq_file = "$ARGV[0]"; ###Output from cleavland file
#>miR1507        CCTCGTTCCATACATCATCTAG  tasi_NO22_chr3_10490780_10491094        715     6       -
$range = 500; #range used in extracting tasi loci sequence
$phase = 21;
$cycle = 10;#	$s = $start - 1 -$range;#	$length = 2*$range + $phase*$cycle + ($end - $start + 1);#	$subseq = substr($chromosome{$chr}, $s, $length);
$target_range_large = 15; #check cleavage in PARE library
$target_range_small = 2;
#read PARE DB
my $user_name = "kakrana";
my $user_password = 'livetheday';
my $database_name = "SOY_priv_PARE";###same as below line
my $PARE_database_name = "SOY_priv_PARE";
my $run_master_table = "run_master";##No need to change
my $tag_info_table = "tag_info";##No need to change
my $summary_table = "kakrana_data.SOY_priv_PARE_summary_table";## from summary script
my $tag_position_table = "tag_position";##No need to change
#Gene Master table
my $gene_master_DB = "SOY_Gmax101_genome";###Genome database

#Small RNA summary2  table
my $sRNA_DB = "SOY_pub_sRNA_v4";###Use tomato sRNA database name, same as used in summary table
my $sRNA_summary_table = "kakrana_data.SOY_priv_sRNA_v4_summary_table ";## if summary table created for sRNA feed here

##phasing score table
#my $phasing_table = "jxzhai.Med_merge_phase_21_8libs";
#my $phase_range = 5000;


#number of tasi tags from each loci
my $tasi_tag_num = 2;

########################################## Main ##############################################################

# connect to PARE DB
&DBI_connection ($database_name,$user_name,$user_password);
print "MySQL database $database_name connected\n";
print "PARE Summary table: $summary_table\n";
print "Small RNA Summary table: $sRNA_summary_table\n";
my $sth = &SQL_row("SELECT * from $summary_table limit 1");  # Pass the SQL statement to the server
my $aref = $sth->{NAME};
my @Column_name = @$aref;
my $library_num = SQL_number("SELECT count(*) from $target_DB.library_info");
print "PARE RNA libraries included ($library_num)\n\n";

open (OUT, ">Output_tasi_target\_$seq_file");

print OUT "miR_name\tmiR_seq\tmiR_len\tmiR_abun\tSubject\tScore\tStrand\tGene\tGene_title\tSub_chr\tCleavage_site\tPhase_start\tPhase_score\tPhase_check"; #\tCleavage_site\tCleavage_strand\tPhasing_check
foreach my $i (5..$#Column_name)
{
	my $code = $Column_name[$i];
	print OUT "\t$code\_s\t$code\_l";
	$abun_filter_stm .= "$code >= $min_abun_each OR ";
}
print OUT "\n";

# print "File opened!\n";

#analyze Cleaveland result (in LINE format)

open (SGOUT, $seq_file);

while (<SGOUT>) {
	
	chomp;
	next if ($_ =~ m/\#/);
	my ($miR,$tag,$subject,$csite,$score,$strand) = split/\t/;
	$miR =~ s/\>//;
	$miR =~ s/\r//;
	$miR =~ s/\n//;
	$subject =~ s/\>//;
	$subject =~ s/\r//;
	$subject =~ s/\n//;
#print "$miR\n";
	if ($subject =~ /Chr(\d+)/i) {	
		$chr = $1;
#		$g_start = $2;
#		$g_end = $3;
	}	
#	$start = $g_start -$range;
#	$length = 2*$range + $phase*$cycle + ($g_end - $g_start + 1);
#	$end = $start + $length - 1;
#	$csite = $start + $csite -1;
#print "$chr\t$csite\t$start\t$start\t$csite\n";
#die;
	$miR_length = length($tag);
	$miR_abun = 0;
	$miR_abun = SQL_number("SELECT acc_abundance_norm from $sRNA_summary_table where tag like \"$tag\" ");
	
#Calculate cleavage site

	$gene = SQL_number("SELECT gene from $gene_master_DB\.gene_master where chr_id = $chr and $csite between start and end limit 1");
	$gene_title = SQL_number("SELECT title from $gene_master_DB\.gene_master where chr_id = $chr and $csite between start and end limit 1");

#	print OUT "$miR\t$tag\t$miR_length\t$miR_abun\t$subject\t$score\t$strand\t$gene\t$gene_title";
#	print "$miR\t$tag\t$miR_length\t$miR_abun\t$subject\t$score\t$strand\t$gene\t$gene_title";
# miR map to "+" strand, cleave on "-" strand
	if ($strand eq "+") {
		$s_start = $subject_start;
		$s_end = $subject_end;
#		print OUT "\t$chr\t$csite";
#		print  "\t$chr\t$csite";
	#5' start position of the PARE tag
		$pare_start = $csite; #on the minus strand
#	#check wether the cleavage is in phase with sRNA
#		$range1 = $csite - $phase_range;
#		$range2 = $csite + $phase_range;
#		$phase_start = &SQL_number ("select Position from $phasing_table where Chr_id = $chr and Position between $range1 and $range2 order by P_score desc limit 1");
#		$phase_score = &SQL_number ("select P_score from $phasing_table where Chr_id = $chr and Position between $range1 and $range2 order by P_score desc limit 1");
#		$phase_check = ($phase_start - ($csite +3 ))%21;
#		print OUT "\t$phase_check";
		print OUT "$miR\t$tag\t$miR_length\t$miR_abun\t$subject\t$score\t$strand\t$gene\t$gene_title\t$chr\t$csite";
#		print "$miR\t$tag\t$miR_length\t$miR_abun\t$subject\t$score\t$strand\t$gene\t$gene_title\t$chr\t$csite";
		print OUT "\t$phase_start\t$phase_score\t$phase_check";
	#confirm cleavage in PARE DB
	$target_start_l = $csite - $target_range_large;
	$target_end_l = $csite + $target_range_large;
	$target_start_s = $csite - $target_range_small;
	$target_end_s = $csite + $target_range_small;
	my $tag_position_ref= &SQL_array("SELECT tag,chr_id, position, strand from $PARE_database_name\.$tag_position_table where chr_id = $chr and position between $target_start_l and $target_end_l and strand like 'c' group by tag, position,strand"); 
	my $target_abun_large = 0;
	my $target_abun_small = 0;
	my %abun_s;
	my %abun_l;
	foreach my $row (@$tag_position_ref)
	{
		my ($tag,$chr_id, $position, $strand) = @$row;
#		print "$tag\t$chr_id\t$position\t$strand";
#print "$tag\t$chr_id\t$position\t$strand\n";

		foreach my $i (5..$#Column_name){
			my $code = $Column_name[$i];
			$abun_s{$code} = 0;
			$abun_l{$code} = 0;
		}

		$abun_ref = &SQL_array("select * from $summary_table where tag = \'$tag\'");
		foreach my $row (@$abun_ref)
		{
			@summary = @$row;
			foreach my $i (5..$#summary)
			{
#				print "\t$summary[$i]";
				$abun_s{$i} += $summary[$i] if ($target_start_s <= $position and $position <= $target_end_s);
				$abun_l{$i} += $summary[$i];
			}
		}
#		print "\n";
	}
	
#	print "$seq_name\t$query_seq\t$target_seq\t$score\t$chr_id\t$start\t$end\t$gene\t\"$gene_func\"";
	foreach my $i (5..$#Column_name) {
		$abun_s{$i} = 0 if (!$abun_s{$i});
		$abun_l{$i} = 0 if (!$abun_l{$i});
		print OUT "\t$abun_s{$i}\t$abun_l{$i}";
#		print  "\t$abun_s{$i}\t$abun_l{$i}";
	}
	print OUT "\n";
#	print  "\n";
}

# miR map to "-" strand, cleave on "+" strand
	if ($strand eq "-") {
		$s_start = $subject_start;
		$s_end = $subject_end;
#		print OUT "\t$chr\t$s_start\t$s_end";
	#5' start position of the PARE tag
		$pare_start = $csite; #on the minus strand
#	#check wether the cleavage is in phase with sRNA
#		$phase_start = $pare_start; #3'overhang, converted to "+" strand, should be the same as PARE start because all phased sRNAs are converted to top strand
#		$phase_check = ($phase_start - $g_start)%21;
#		print OUT "\t$pare_start\t+\t$phase_check";
#		$range1 = $csite - $phase_range;
#		$range2 = $csite + $phase_range;
#		$phase_start = &SQL_number ("select Position from $phasing_table where Chr_id = $chr and Position between $range1 and $range2 order by P_score desc limit 1");
#		$phase_score = &SQL_number ("select P_score from $phasing_table where Chr_id = $chr and Position between $range1 and $range2 order by P_score desc limit 1");
#		$phase_check = ($phase_start - $csite)%21;
		print OUT "$miR\t$tag\t$miR_length\t$miR_abun\t$subject\t$score\t$strand\t$gene\t$gene_title\t$chr\t$csite";
		print OUT "\t$phase_start\t$phase_score\t$phase_check";
	#confirm cleavage in PARE DB
	$target_start_l = $csite - $target_range_large;
	$target_end_l = $csite + $target_range_large;
	$target_start_s = $csite - $target_range_small;
	$target_end_s = $csite + $target_range_small;
	my $tag_position_ref= &SQL_array("SELECT tag,chr_id, position, strand from $tag_position_table where chr_id = $chr and position between $target_start_l and $target_end_l and strand like 'w' group by tag, position,strand"); 
	my $target_abun_large = 0;
	my $target_abun_small = 0;
	my %abun_s;
	my %abun_l;
	foreach my $row (@$tag_position_ref)
	{
		my ($tag,$chr_id, $position, $strand) = @$row;
#		print "$tag\t$chr_id\t$position\t$strand";

		foreach my $i (5..$#Column_name){
			my $code = $Column_name[$i];
			$abun_s{$code} = 0;
			$abun_l{$code} = 0;
		}

		$abun_ref = &SQL_array("select * from $summary_table where tag = \'$tag\'");
		foreach my $row (@$abun_ref)
		{
			@summary = @$row;
			foreach my $i (5..$#summary)
			{
#				print "\t$summary[$i]";
				$abun_s{$i} += $summary[$i] if ($target_start_s <= $position and $position <= $target_end_s);
				$abun_l{$i} += $summary[$i];
			}
		}
#		print "\n";
	}
	
#	print "$seq_name\t$query_seq\t$target_seq\t$score\t$chr_id\t$start\t$end\t$gene\t\"$gene_func\"";
	foreach my $i (5..$#Column_name) {
		$abun_s{$i} = 0 if (!$abun_s{$i});
		$abun_l{$i} = 0 if (!$abun_l{$i});
		print OUT "\t$abun_s{$i}\t$abun_l{$i}";
	}
	print OUT "\n";

}

}

#close TEMP;
#$/ = ">";


close OUT;
##tasi_NO1	chr_1	phased_1	29	4433752-4433752
#	my ($tasi_num, $chr, $phase_num, $abun, $pst) = split /\t/;
#	$chr =~ s/chr\_//i;
#	my ($start,$end) =split/-/,$pst;
#	$s = $start - 1 -$range;
#	$length = 2*$range + $phase*$cycle + ($end - $start + 1);
#
#>ath-miR390b	4	15	0.0001	1	17	17	2	2	0	+	510	528	>tasi_NO4_chr2_11490602_11490623	M1;D1;M15;S1;D1;M1;S1;
#>ath-miR390b	4	8	0.2471	10	21	18	2	1	1	+	737	757	>tasi_NO4_chr2_11490602_11490623	M8;I1;S1;M3;D1;M5;S1;M2;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub SQL_row
{
	my $SQL_command = $_[0];
	my $SQL_select = $dbh->prepare($SQL_command) or die "Could not prepare SQL statement.\n";
	$SQL_select->execute or die "Could not execute SQL statement.\n";
	return ($SQL_select);
	#while ( my @row= $SQL_select->fetchrow_array())
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub SQL_array
{	
	my $SQL_command = $_[0];
	my $SQL_select = $dbh->prepare($SQL_command) or die "Could not prepare SQL statement.\n";
	$SQL_select->execute or die "Could not execute SQL statement.\n";
	my $array_ref = $SQL_select->fetchall_arrayref(); 
	return ($array_ref);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub SQL_number
{	
	my $SQL_command = $_[0];
	my $SQL_select = $dbh->prepare($SQL_command) or die "Could not prepare SQL statement  .\n";
	$SQL_select->execute or die "Could not execute SQL statement. $SQL_command\n";
	my $number;
	while(my @array = $SQL_select->fetchrow_array()) 
	{
		$number = $array[0];
	}
	return $number;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub DBI_connection
{	($database_name,$user_name,$user_password) = @_;
	$dbh = DBI->connect("dbi:mysql:$database_name:host=localhost:port=3306", "$user_name", "$user_password") or die "Could not connect to database $database_name.\n";
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub TIME
{   my $t0 = $_[0];
	my $t1 = new Benchmark;
    my $td = timediff($t1, $t0);
    print "\n(Time for code execution :",timestr($td),")\n";
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub MAX
{	# find max value
	# call &MAX(default value, other values . . . )
	my ($m,@l) = @_;
    foreach my $x (@l) { $m = $x if ($x > $m); }
    return $m;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
