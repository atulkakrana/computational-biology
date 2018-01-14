#!usr/bin/perl -w
###Jixian original script can be found in dropbox:00-miRNA:soyafolder

############# By Jixian, Modified to do the Genome-wide analysis, mapping all small RNAs together, turned off random shuffles to increase speed, added calculation for matching scores on the "+" strand
use Getopt::Std;
### CRITICAL:  Must set the location of the SHRiMP rmapper-ls executable here ... will differ for each system###Find for pikachu
##$shrimp = "/gpfs/home/mja18/work/SHRiMP_1_3_1/bin/rmapper-ls";
$shrimp = "/usr/local/shrimp/bin/rmapper-ls";###Get it installed and change the address here
###
################## -c is the genome index file
$usage = "-q : The FASTA file of query sequences: All seqs must be 20-25 in length and all characters acgutACGUT\n" . 
    "-c : FASTA file of contigs: usuall transcriptome only revcomp hits will be search\n" .
    "-s : Maximum score cutoff\.  DEFAULT: 7 \(all alignments will be kept\)\n" . 
    "-r : Number of random shuffles to also predict\.  Full alignments not retained for these\.  DEFAULT: 100\n";

getopt('qcsr');
## check required files
unless(-r $opt_c) {
#    die "FATAL:  Could not read the contig file $opt_c\nUSAGE:\n$usage\n";
}
unless(-x $shrimp) {
    die "FATAL:  The SHRiMP rmapper-ls $shrimp was not executable\.  Open this script and edit \$shrimp on line 5\nUSAGE:\n$usage\n";
}
(open(QFILE, "$opt_q")) || die "FATAL:  Could not open the query file $opt_q check option -q\n$usage\n";
@q_fa_data = <QFILE>;
close QFILE;

## set options

if($opt_s) {
    $max_score = $opt_s;
} else {
    $max_score = 7;
}

if($opt_r) {
    $n_shuffles = $opt_r;
} else {
    $n_shuffles = 100;
}

## go through the query FASTA file one query at a time
open(TQ, ">TEMP_query\_$opt_q");
foreach $fa_line (@q_fa_data) {
    chomp $fa_line;
    if($fa_line =~ /^>/) {
	$qname = $fa_line;
	$qname =~ s/\|/:/g;  ## get rid of pipe metacharacters in the name
	$qname =~ s/\s/--/g;  ## get rid of whitespace characters in the name
	print TQ "$qname\n";
    } else {
	$qseq = uc $fa_line;
        ## verify query sequence and its name

	$qseq =~ s/U/T/g;  ## convert to DNA form
	unless(($qseq =~ /^[ACTG]+$/) and
	       ((length $qseq) >= 18) and
	       ((length $qseq) <= 25)) {
	    die "Query seqeunce $qseq is not acceptable\n$usage\n";
	}
	unless($qname =~ /^>/) {
	    $qname = ">" . "$qname";
	}
	print TQ "$qseq\n";
	$HASH_miR{$qname}=$qseq;
	}
}
        ## open outfiles using a standardized naming format.  Fails if the files already exist
	$real_out = "$opt_q" . "_$opt_c" . "_REAL_targets_$max_score";
#	$ran_out = "$opt_q" . "_$opt_c" . "_RAN$n_shuffles" . "_targets_$max_score";
	$real_out =~ s/>//g;
#	$ran_out =~ s/>//g;
	
#	if((-e $real_out) or
#	   (-e $ran_out)) {
#	    die "output file $real_out and:or $ran_out already exists FATAL\n$usage\n";
#	} else {
	    open(REALOUT, ">$real_out");
#	    open(RANOUT, ">>$ran_out");
		open (LINEOUT, ">LINE\_$real_out");
		open (LOG, ">LOG\_$real_out");
#	}
	
	print LOG "axtell_targetfinder_2_batchmode\.pl\n";
	print LOG "Date: ";
	print LOG `date`;
	print LOG "Hostname: ";
	print LOG `hostname`;
	print LOG "Working Directory: ";
	print LOG `pwd`;
	print LOG "SHRiMP rmapper-ls program: $shrimp\n";
	print LOG "-q Query fasta file: $opt_q\n";
#	print LOG "Query \(DNA form\): $qseq\n";
#	print LOG "Query name: $qname\n";
	print LOG "-c Contigs\/Transcriptome: $opt_c\n";
	print LOG "-s Maximum score allowed: $max_score\n"; 
	print LOG "-r Number of random queries examined: $n_shuffles\n";
	print LOG "Query results file: $real_out\n";
	#print LOG "Randomized queries results file: $ran_out\n";

        ## write the query to a temporary fasta file for use in calling rmapper-ls
#	$qtempfile = "$qname" . "_temp";
#	$qtempfile =~ s/>//g;  ## don't want those symbols on file names
#	open(TQ, ">$qtempfile");
#	print TQ "$qname\n$qseq\n";
	close TQ;
	
        # NEW 3/4/10
	$LOGtemp = "$opt_q" . "_$opt_c". "_LOGtemp";
	$shrimpout = "$opt_q" . "_$opt_c" . "_shrimptemp";
	print "$shrimpout\n";
	#print LOG "Running $shrimp on $qtempfile ...";
	#system "$shrimp -s 11100111,11001100011,11001001001001,10010001000010011 -n 1  -g -16 -q -16 -e -10 -f -10 -m 15 -i -10 -h 33% -w 175% -o 10000 -N 12 -K 1 -P  TEMP_query\_$opt_q $opt_c > $shrimpout 2> $LOGtemp";
	########################## Use Pre-Indexed genome
	system "$shrimp -s 11100111,11001100011,11001001001001,10010001000010011 -n 1 -t -1 -g -16 -q -16 -e -10 -f -10 -m 15 -i -10 -h 33% -w 175% -o 10000 -P   TEMP_query\_$opt_q $opt_c  > $shrimpout 2> $LOGtemp";
#	system "rm -f $LOGtemp";
	#print $shrimpout;
	#die;
	(open(SH, "$shrimpout")) || die "FATAL could not open shrimp output file $shrimpout\n";
	@sh_data = <SH>;
	close SH;
	#system "rm -f $shrimpout";
	
	print LOG "Done\nParsing data...";
	
        ## parse the data, and print the results to REALOUT
	$full = 1;  ## indicated to the sub-routine whether or not to print full alignments or just one line information
	($lineresults, $real_results) = shrimptar_parse (\@sh_data,\$max_score,\$full);
	print LINEOUT "$lineresults\n";
	print REALOUT "$real_results\n";
#	print  "$real_results\n";

        ## examine the real hits and summarize in LOG...
	@real_res = split ("\n", $real_results);
	%rl = ();
	%final = ();
	$running_tally = 0;
	## [0] : real number
	## [1] : cumulative real number
	## [2] : random mean number
	## [3] : random stdev
	## [4] : random cumulative mean number
	## [5] : random cumulative stdev number
	
	for($i = 0; $i <= $max_score; $i += 0.5) {
	    $rl{$i} = 0;
	}
	foreach $rr (@real_res) {
	    chomp $rr;
	    if($rr =~ /^>/) {
		@fields = split ("\t", $rr);
		++$rl{$fields[4]};
	    }
	}
	for($i = 0; $i <= $max_score; $i += 0.5) {
	    $running_tally += $rl{$i};
	    push (@{$final{$i}}, $rl{$i});
	    push (@{$final{$i}}, $running_tally);
	}
	print "$running_tally\n";
	print "Done\n";

	
	

#     ## now, make the randomized versions
#     ## print to temp files
#     ## do in blocks of 100
#     ## also, store the random hit frequencies in a hash
#	%sh_hash = ();
#	$x = 0;
#	$tempfa = "ran_" . "$qtempfile";
#	$LOGtemp = "$tempfa" . "_LOGtemp";
#	$shrimpout = "$tempfa" . "_shrimptemp";
#	open(TEMPFA, ">$tempfa");
#	for($i = 1; $i <= $n_shuffles; ++$i) {
#	    ## make a random permutation of the query and write to fasta file
#	    $random = '';
#	    @letters = split ('', $qseq);
#	    while (@letters) {
#		$size = scalar @letters;
#		$index = int(rand($size));
#		$random .= $letters[$index];
#		splice (@letters, $index, 1);
#	    }
#	    print TEMPFA "$qname" . "_ran$i\n$random\n";
#	    ++$x;
#	    if($x == 100) {
#		print LOG "Searching against a set of 100 randomized queries...\n";
#		$x = 0;
#		close TEMPFA;
#		
#		### run the search
#		#@sh_data = shrimptar ($shrimp, $tempfa, $opt_c);
#		system "$shrimp -s 11100111,11001100011,11001001001001,10010001000010011 -n 1  -g -16 -q -16 -e -10 -f -10 -m 15 -i -10 -h 33% -w 175% -o 10000 -N 8 -P $tempfa $opt_c > $shrimpout 2> $LOGtemp";
#		system "rm -f $LOGtemp";
#		(open(SH, "$shrimpout")) || die "FATAL could not open shrimp output file $shrimpout\n";
#		@sh_data = <SH>;
#		close SH;
#		system "rm -f $shrimpout";
#		print LOG "Done\nParsing the results...";
#
#		### parse, sort, and output  the results
#		$full = 0;
#		$unsorted_ran = shrimptar_parse (\@sh_data, \$max_score, \$full);
#		@us_randata = split ("\n", $unsorted_ran);
#		@randata = sort @us_randata;
#		@us_randata = (); ## just to save some memory
#		$sorted_ran = join ("\n", @randata);
#		print RANOUT "$sorted_ran\n";
#		
#		## parse the shuffle data for the sh_hash
#		parseran(\@randata, \%sh_hash, \$max_score);
#		
#		### delete the old tempfile, and open a new one
#		system "rm -f $tempfa";
#		open(TEMPFA, ">$tempfa");
#		print LOG "Done\n";
#	    }
#	}
	
## run the search if there are any in the tempfa file
#	if($x) {
#	    print LOG "Search against final set of randomized queries...\n";
#	    ## run the search
#	    system "$shrimp -s 11100111,11001100011,11001001001001,10010001000010011 -n 1 -t -1 -g -16 -q -16 -e -10 -f -10 -m 15 -i -10 -h 33% -w 175% -o 10000 -C -P $tempfa $opt_c > $shrimpout 2> $LOGtemp";
#	    system "rm -f $LOGtemp";
#
#	    (open(SH, "$shrimpout")) || die "FATAL could not open shrimp output file $shrimpout\n";
#	    @sh_data = <SH>;
#	    close SH;
#	    system "rm -f $shrimpout";
#	    print LOG "Done\nParsing the results...";
#
#	    ### parse, sort, and output  the results
#	    $full = 0;
#	    $unsorted_ran = shrimptar_parse (\@sh_data, \$max_score, \$full);
#	    @us_randata = split ("\n", $unsorted_ran);
#	    @randata = sort @us_randata;
#	    @us_randata = (); ## just to save some memory
#	    $sorted_ran = join ("\n", @randata);
#	    print RANOUT "$sorted_ran\n";
#	    
#	    ## parse the shuffle data for the sh_hash
#	    parseran(\@randata, \%sh_hash, \$max_score);
#	    
#	}
### delete the old tempfile
#	system "rm -f $tempfa";
#	print LOG "Done\n";
### report on the noise    
## first get the simple means and SDs for each score
#	for($i = 0; $i <= $max_score; $i += 0.5) {
#	    $mean = get_mean (@{$sh_hash{$i}});
#	    $sd = get_sd (@{$sh_hash{$i}});
#	    push (@{$final{$i}}, $mean);
#	    push (@{$final{$i}}, $sd);
#	}
#	
### now calculate the cumulative noise at each score
#	%cumul = ();
#	for ($x = 0; $x < $n_shuffles; ++$x) {
#	    $running_tally = 0; 
#	    for($i = 0; $i <= $max_score; $i += 0.5) {
#		if(${$sh_hash{$i}}[$x]) {
#		    $running_tally += ${$sh_hash{$i}}[$x];
#		}
#		push (@{$cumul{$i}}, $running_tally);
#	    }
#	}
#	
#	for ($i = 0; $i <= $max_score; $i += 0.5) {
#	    $mean = get_mean (@{$cumul{$i}});
#	    $sd = get_sd (@{$cumul{$i}});
#	    push (@{$final{$i}}, $mean);
#	    push (@{$final{$i}}, $sd);
#	}
#	
## finally, output the summary to LOG
	print LOG "\n";
	print LOG "SCORE\tQUERY\tQUERY_CUMULATIVE\tRANDOM_MEAN\tRANDOM_STDEV\tRANDOM_CUMULATIVE_MEAN\tRANDOM_CUMULATIVE_STDEV\n";
	for ($i = 0; $i <= $max_score; $i += 0.5) {
	    print LOG "$i\t";
	    $out = join ("\t", @{$final{$i}});
	    print LOG "$out\n";
	}
	
# delete the last tempfile
#	system "rm -f $qtempfile";
#   }
#}

### sub-routines

sub shrimptar_parse {
    my($sh_data, $max_score, $full) = @_;  ## passed by reference .. sh_data is an array, $max_score is a string
    ## initialize some variables
    my $sh_line;
    my $dataline;
    my @fields = ();
    my $q_name;
    my $hit_name;
    my $strand;
    my $contigstart;
    my $contigstop;
    my $xcontigstart;
    my $xcontigstop;
    my $querystart;
    my $querystop;
    my $querylen;
    my $contigseq_al;
    my $qseq_al;
    my $contigseq_for;
    my $qseq_for;
    my $nfivep_gaps;
    my $fivetoadd;
    my $fivetoadd_rev;
    my $nthree_gaps;
    my $threetoadd_rev;
    my $threetoadd;
    my $al_line;
    my @qal = ();
    my @cal = ();
    my $pos_query;
    my $score;
    my $cgaps_1stten;
    my $qgaps_1stten;
    my $i;
    my $qbase;
    my $cbase;
    my $this_score;
    my $this_score_f;
    my $al_line_string;
    my $magic_pos;
    my $final_result;
    my @G = ();
    my @R = ();
    my $g;
    my $r;
    my $first_letter;
    my $last_letter;
    my $gcount;
    
    # initialize a few simple hashes for scoring matrices
    my %wc = (   ## Watson-Crick matches for RNA
		 'A' => 'U',
		 'G' => 'C',
		 'C' => 'G',
		 'U' => 'A',
	);

    my %gu = (  ## GU wobbles
		'G' => 'U',
		'U' => 'G',
	);
    
    # here is the parsing section
    foreach $sh_line (@$sh_data) {
	chomp $sh_line;
#	print "$sh_line\n";
######################by Jixian, need to deal with matching on both strands
	if($sh_line =~ /^>/) {  ## this is the shrimp data line.  all queries should begin with the > fasta header character
	    $dataline = $sh_line;
	    @fields = split ("\t", $sh_line);
	    $q_name = $fields[0];
	    $hit_name = $fields[1];
	    $strand = $fields[2];
#	    unless($strand eq "-") {
#		die "Strand $strand is not - at entry $sh_line  \n";
#	    }
	    $xcontigstart = $fields[3];  ## major bug in earlier versions .. the SHRiMP numbers here are NOT the all of the contig bases shown!
	    $xcontigstop = $fields[4];  ## as above
	    $querystart = $fields[5];
	    $querystop = $fields[6];
	    $querylen = $fields[7];
	    
	} elsif ($sh_line =~ /^G:.*\s([ATCGURYMKWSBDHVN-]+)\s/) {  ## ACCOUNTS FOR AMB CHARACTERS
	    $contigseq_al = $1;
	} elsif ($sh_line =~ /^R:.*\s([ATCG-]+)\s/) {  ## NO AMB CHARACTERS EXPECTED IN QUERIES
	    $qseq_al = $1;
	    
	    ## get the real contig start and stop numbers
	    @G = split ('', $contigseq_al);
	    @R = split ('', $qseq_al);
	    $gcount = 0;
	    $last_letter = 0;
	    $first_letter = 0;
	    foreach $g (@G) {
		$r = shift @R;
		if($g =~ /[ACTG]/) {
		    ++$gcount;
		    if($r eq $g) {
			unless($first_letter) {
			    $first_letter = $gcount;
			}
			$last_letter = $gcount;
		    }
		}
	    }
		################################## Modified by Jixian ##############################################
	    if ($strand eq "-") {
		$contigstop = $xcontigstop + $first_letter - 1;
	    $contigstart = $xcontigstart - ($gcount - $last_letter);
	    } elsif ($strand eq "+") {
	    $contigstart = $xcontigstart - $first_letter + 1;
		$contigstop = $xcontigstop + ($gcount - $last_letter);
		}
	    #print "$contigstop\t$contigstart\n";
	    
	    ## convert to RNA-RNA alignment
	    ## for the contig, simply get the reverse complement of the aligned bases, then convert Ts to Us
	    $contigseq_for = reverse $contigseq_al;
	    $contigseq_for =~ tr/ACTGRYMKWSBDHVNU/UGACYRKMWSVHDBNA/;  ## FOR AMB CODES IN CONTIGS
	    ## for the query, reverse it and convert it to RNA, but NOT the complement
	    $qseq_for = reverse $qseq_al;
	    $qseq_for =~ s/T/U/g;  ## 3' to 5' direction, RNA version
	    ## fill in any gaps on the 5' end of the query (now at the end of the string)
	    $nfivep_gaps = $querystart - 1;
	    if($nfivep_gaps > 0) {
		$fivetoadd_rev = substr($HASH_miR{$q_name},0,$nfivep_gaps);
		$fivetoadd = reverse $fivetoadd_rev;
		$fivetoadd =~ s/T/U/g;
		$qseq_for =~ s/-+$/$fivetoadd/;
	    }
	    ## fill in any gaps on the 3' end of the query (now at the beginning of the string
	    $nthreep_gaps = $querylen - $querystop;
	    if($nthreep_gaps > 0 ) {
		$threetoadd_rev = substr($HASH_miR{$q_name},$querystop,$nthreep_gaps);
		$threetoadd = reverse $threetoadd_rev;
		$threetoadd =~ s/T/U/g;
		$qseq_for =~ s/^-+/$threetoadd/;
	    }
	    ## score the RNA-RNA alignment, using Allen et al rules
	    $al_line = ();  ## | for W/C match, o for G:U wobble
	    @qal = split ('', $qseq_for);
	    @cal = split ('', $contigseq_for);
	    unless ((scalar @qal) == (scalar @cal)) {
		die "unequal sizes of qseq_for $qseq_for and contigseq_for $contigseq_for\n";
	    }
	    $pos_query = 0;
	    $score = 0;
	    $cgaps_1stten = 0;
	    $qgaps_1stten = 0;
	    for($i = ((scalar @cal) - 1); $i >= 0; --$i) {
		$qbase = $qal[$i];
		$cbase = $cal[$i];
		if($qbase =~ /[ACUG]/) {
		    ++$pos_query;
		}
		
		if(($qbase =~ /[^ACUG]/) or ($cbase =~ /[^ACUG]/)) {  ## one or both is a gap or ambig.
		    $this_score = 1;
		    $al_line .= " ";
		    
		} elsif ((exists($gu{$qbase})) and ($gu{$qbase} eq $cbase)) { ## it is a G:U wobble
		    $this_score = 0.5;
		    $al_line .= "o";
		} elsif ($wc{$qbase} eq $cbase) { ## it is a W:C match
		    $this_score = 0;
		    $al_line .= "\|";
		} else {   ## it is a mismatch
		    $this_score = 1;
		    $al_line .= " ";
		}
		## penalty doubled at query positions 2-13
		if(($pos_query >= 2) and ($pos_query <= 13)) {
		    $this_score_f = 2 * $this_score;
		} else {
		    $this_score_f = $this_score;
		}
		## track gaps in the first ten positions
		if($pos_query <= 10) {
		    if($qbase eq "-") {
			++$qgaps_1stten;
		    }
		    if($cbase eq "-") {
			++$cgaps_1stten;
		    }
		}
		$score += $this_score_f;
	    }
	    $al_line_string = reverse $al_line;
	    
	    ## find the position of the contig across from the tenth base of the query
		############################# modified by Jixian ################################
		if ($strand eq "-") {
			$magic_pos = $contigstop - 9 + $cgaps_1stten - $qgaps_1stten;
		}elsif ($strand eq "+") {
			$magic_pos = $contigstart + 9 - $cgaps_1stten + $qgaps_1stten;
		}

	    
	    ## mark the cleavage position
		################## modified by Jixian
	    $mark_line = '';
		if ($strand eq "-") {
	    $cpos = $contigstart - 1 ;
	    foreach $cal_character (@cal) {
		if($cal_character =~ /[AUCG]/) {
		    ++$cpos;
		}
		if(($cpos == $magic_pos) and ($cal_character =~ /[AUCG]/)) {
		    $mark_line .= "\*";
		} else {
		    $mark_line .= " ";
		}
	    }
		}
	    elsif ($strand eq "+") {
		$cpos = $contigstop + 1 ;
	    foreach $cal_character (@cal) {
		if($cal_character =~ /[AUCG]/) {
		    --$cpos;
		}
		if(($cpos == $magic_pos) and ($cal_character =~ /[AUCG]/)) {
		    $mark_line .= "\*";
		} else {
		    $mark_line .= " ";
		}
		}
		}
    
	    ## output the results, for the real query
	    if($score <= $$max_score) {
		$line_result .= "$q_name\t$HASH_miR{$q_name}\t$hit_name\t$magic_pos\t$score\t$strand\n";
		$final_result .= "$q_name\t$HASH_miR{$q_name}\t$hit_name\t$magic_pos\t$score\t$strand\n";
		if($$full) {
		    $final_result .= "$mark_line\n";
		    $final_result .= "$contigseq_for\t$hit_name 5\'-3\' $contigstart to $contigstop CLEAVAGE SITE: $magic_pos\n";
		    $final_result .= "$al_line_string\tSCORE: $score\n";
		    $final_result .= "$qseq_for\t$q_name 3\'-5\'\n\n";
		}
	    }
	}
    }
#	print "$final_result\n";
    return ($line_result,$final_result);
}

sub parseran {
    my($randata,$sh_hash,$max_score) = @_;  ## passed by reference.  randata=@, sh_hash=%, max_score=$
    my %single = ();
    my %already = ();
    my $dataline;
    my @fields = ();
    my $i;
    my $last = "null";
    
    
    for($i = 0; $i <= $$max_score; $i += 0.5) {
	$single{$i} = 0;
    }
    foreach $dataline (@$randata) {
	chomp $dataline;
	@fields = split ("\t", $dataline);
	if($fields[0] ne $last) {
	    
	    ## check to ensure that the new query has not been encountered before.  if it has, die, because the input data was not properly sorted
	    if(exists($already{$fields[0]})) {
		die "FATAL in sub-routine parseran: data not sorted properly: $fields[0] encountered more than once\n";
	    }
	    
	    ## now make the switches
	    $last = $fields[0];
	    $already{$last} = 1;
	    
	    ## parse the previous data
	    
	    for($i = 0; $i <= $$max_score; $i += 0.5) {
		push (@{$sh_hash{$i}}, $single{$i});
		$single{$i} = 0;  ## clear the hash once the tally is recorded
	    }
	    
	}
	++$single{$fields[3]};
    }
}

sub get_mean {
    my(@info) = @_;
    my $denom = scalar @info;
    my $sum = 0;
    my $datapoint;
    foreach $datapoint (@info) {
	$sum += $datapoint;
    }
    my $mean = $sum / $denom;
    my $mean3 = sprintf ("%.3f", $mean);
    return $mean3;
}

sub get_sd {
    my(@info) = @_;
    my $denom = scalar @info;
    my $sum = 0;
    my $datapoint;
    foreach $datapoint (@info) {
	$sum += $datapoint;
    }
    my $mean = $sum / $denom;
    my $total2 = 0;
    foreach $datapoint (@info) {
	$total2 += ($mean - $datapoint)**2;
    }
    my $mean2 = $total2 / $denom;
    my $stdev = sqrt($mean2);
    my $stdev3 = sprintf ("%.3f", $stdev);
    return $stdev3;
}
    
