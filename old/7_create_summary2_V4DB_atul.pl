#!/usr/bin/perl

#### Last modified by Guna on Nov 20th, 2009 ####

$USAGE = "\nUSAGE: perl 7_create_summary2_v4DB.pl\n";

use Time::HiRes qw(gettimeofday);
use DBI;


### Modify Hardcoded input parameters (for each run):
$database_name = "$ARGV[0]";##PARE database name/sRNA database or any expression database
$user_name = "kakrana";
$user_password = 'livetheday';
$output_table = "kakrana_data.$ARGV[0]\_summary_table";###keep in your database 'jxzhai_raichu' is a database name. Get your own database
$server = "pikachu.dbi.udel.edu";###Change the server to pikachu and get all the databases required on pikachu, PARE and genome database
$port = 3306; # 3307 for secondary MySQL server
			  # 3306 for default primary MySQL
### End custom parameters


$dbh = DBI->connect("dbi:mysql:$database_name;host=$server;port=$port", "$user_name", "$user_password") or die "Could not connect to database.\n";


($date, $time) = &time_stamp();
print "Creating $output_table table... ($date, $time)\n";

$select_library_codes = $dbh->prepare( "SELECT code, lib_id FROM library_info" ) or die "Could not prepare SQL statement.\n"; # where lib_id IN (154,155,156,157)
$select_library_codes->execute or die "Could not execute SQL statement.\n";

$num_libs = 0;

while(@array = $select_library_codes->fetchrow_array()) 
{
	$code = $array[0];
	$code =~ s/-/_/g;
	$code =~ s/\./_/g;
	$lib_id = $array[1];
	#$code = $code."_".$lib_id;
	if($lib_id == 0 or $code eq '') { next; }
	$codes_by_lib{$lib_id} = $code;
	$codes[$num_libs] = $code;
	$lib_ids[$num_libs] = $lib_id;
	$num_libs++;
}

$stm = "CREATE TABLE $output_table ( tag varchar(50) NOT NULL, len int(19) NOT NULL DEFAULT 0, hits int(10) NOT NULL DEFAULT 0, rna_type SET('miRNA', 'natRNA', 'pRNA', 'piRNA', 'pre_tRNA', 'qiRNA', 'rRNA', 'rasiRNA', 'siRNA', 'slRNA', 'snRNA', 'snoRNA', 'tRNA', 'tasiRNA', 'misc_RNA', 'other_RNA') not null default '', acc_abundance_norm int(50) NOT NULL DEFAULT 0";


for($i = 0; $i < $num_libs; $i++) 
{
	$stm .= ", $codes[$i] int(10) NOT NULL DEFAULT 0";
}

$stm .= ")";

#print "$stm\n";


$drop_summary2 = $dbh->prepare( "DROP TABLE IF EXISTS $output_table" ) or die "Could not prepare SQL statement.\n";
$drop_summary2->execute or die "Could not execute SQL statement.\n";

$create_summary2 = $dbh->prepare( "$stm" ) or die "Could not prepare SQL statement.\n";
$create_summary2->execute or die "Could not execute SQL statement.\n";

$select_tags = $dbh->prepare( "SELECT tag, len, hits, rna_type+0 FROM tag_info" ) or die "Could not prepare SQL statement.\n";
$select_tags->execute or die "Could not execute SQL statement.\n";


$count_tags = 0;
while(@array = $select_tags->fetchrow_array()) 
{
	$t0 = gettimeofday();
	
	$tag = $array[0];
	if($tag eq '') { next; }
	$len = $array[1];
	$hits = $array[2];
	$rna_type = $array[3];
	#raw_value should be norm
	$select_norms = $dbh->prepare( "SELECT lib_id, norm FROM run_master WHERE tag = \'$tag\'" ) or die "Could not prepare SQL statement.\n";
	$select_norms->execute or die "Could not execute SQL statement.\n";
	
	$stm = "INSERT INTO $output_table (tag, len, hits, rna_type";
	
	$stm2 = ") VALUES (\'$tag\', $len, $hits, $rna_type";
	
	$acc_abundance_norm = 0;
	while(@array2 = $select_norms->fetchrow_array() ) 
	{
		$lib_id = $array2[0];
		$norm = $array2[1];
		$acc_abundance_norm = $acc_abundance_norm + $norm;
		$stm .= ", $codes_by_lib{$lib_id}";
		$stm2 .= ", $norm";
	#	print "$lib_id\t$norm\t$acc_abundance_norm\n";
	}
	
	$stm = $stm . ", acc_abundance_norm" . $stm2 . ", $acc_abundance_norm)";
	
	#print "$stm\n\n";

	$insert_into_summary2 = $dbh->prepare( "$stm" ) or die "Could not prepare SQL statement.\n";
	$insert_into_summary2->execute or die "Could not execute SQL statement.\n";
	$count_tags++;
	
	$t1 = gettimeofday();
	$elapsed_time = $t1 - $t0;
	
	
}

($date, $time) = &time_stamp();
print "$output_table table populated with $count_tags entries... ($date, $time)\n";
print "Creating indexes... ($date, $time)\n";

$create_index = $dbh->prepare( "ALTER TABLE $output_table ADD PRIMARY KEY (tag)" ) or die "Could not prepare SQL statement.\n";
$create_index->execute or die "Could not execute SQL statement.\n";

$create_index = $dbh->prepare( "CREATE INDEX len ON $output_table (len)" ) or die "Could not prepare SQL statement.\n";
$create_index->execute or die "Could not execute SQL statement.\n";

$create_index = $dbh->prepare( "CREATE INDEX hits ON $output_table (hits)" ) or die "Could not prepare SQL statement.\n";
$create_index->execute or die "Could not execute SQL statement.\n";

$create_index = $dbh->prepare( "CREATE INDEX rna_type ON $output_table (rna_type)" ) or die "Could not prepare SQL statement.\n";
$create_index->execute or die "Could not execute SQL statement.\n";

$create_index = $dbh->prepare( "CREATE INDEX acc_abundance_norm ON $output_table (acc_abundance_norm)" ) or die "Could not prepare SQL statement.\n";
$create_index->execute or die "Could not execute SQL statement.\n";

#$insert_into_db_log = $dbh->prepare( "INSERT INTO db_log (updates) VALUES ('Executed create_summary2.pl using $ARGV[0]')" ) or die "Could not prepare SQL statement.\n";
#$insert_into_db_log->execute or die "Could not execute SQL statement.\n";

($date, $time) = &time_stamp();
print "Done! ($date, $time)\n";


sub time_stamp 
{
	my ($d,$t);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	$year += 1900;
	$mon++;
	$d = sprintf("%4d-%2.2d-%2.2d",$year,$mon,$mday);
	$t = sprintf("%2.2d:%2.2d:%2.2d",$hour,$min,$sec);
	return($d,$t);
}
