#!/usr/bin/perl

#package SrnaTools::Module::HpTool ;
#use base SrnaTools::Module ;
use strict ;
use warnings;
use File::Temp qw( tempfile tempdir );
use IO::String ;
use Cwd;
use File::Copy;
use List::Util 'shuffle';


sub run{
# Commands 
#my $rnaplot_ps_macro = '1 8 omark ' ;
#my $RNAfold_cmd = 'RNAfold --MEA -d2 -p' ;
my $RNAfold_cmd = 'RNAfold' ; 
my $RNAplot_cmd = 'RNAplot' ;
my $ps2pdf_cmd = 'ps2pdf' ; 
my $ps2eps_cmd = 'ps2eps' ;
my $working_dir = getcwd;
my $bcount = 0; ### Count matched sRNAs
# my $allcount = 0; ### count all sRNAs 
#print $working_dir;
#my $long_seq_param = $self->param('longSeq');
#my $short_seqs_param = $self->param('shortSeqs');
my $cwd = getcwd;
($cwd) =($cwd=~/^(.+)$/);
#print $cwd;
# rgb colors for highlighting short seqs
my @rgb_colors = ('255 40 40', 
                  '255 0 255', 
                  '0 0 255', 
                  '0 255 255', 
                  '0 255 0', 
                  '255 255 0', 
                  '255 127 0',
                  '255 184 235',
                  '181 255 184',
                  '199 255 176',
                  '127 62 62',
                  '255 227 99',
                  '127 0 0', 
                  '127 0 127', 
                  '0 0 127',
                  '0 127 127',
                  '0 127 0',
                  '130 127 0',
                  '153 153 153',
                  '204 204 204',
                  '62 130 153', 
                  '0 62 127',
                  '62 127 127',
                  '130 127 62',
                  '130 127 0',
                  '62 127 153',
                  '123 0 204',
                  '155 284 235',
                  '81 225 215',
                  '129 155 126',
                  '227 62 162',
                  '105 87 199',
                  '127 200 100', 
                  '227 105 157', 
                  '110 205 217',
                  '190 227 127',
                  '210 217 205',
                  '230 92 160',
                  '180 190 53',
                  '104 232 104',
                  '93 78 223',
                  '42 187 67',
                  '68 23 53',
                  '219 92 104',
                  '254 129 253', 
                  '89 167 27',
                  '176 28 227',
                  '180 67 162',
                  '234 192 109',
                  '32 156 253',
                  '66 65 214',
                  '198 244 25',
                  '172 55 225',
                  '97 255 226',
                  '183 38 62',
                  '8 87 19',
                  '0 219 119', 
                  '201 105 127', 
                  '109 173 117',
                  '164 227 27',
                  '233 21 25',
                  '164 92 60',
                  '22 37 153',) ;
  
  my $long_seq = $_[0];#"TGTACTGCCCTTTCTCCATCCCCCAAATCTTTTGGAGTTTTTAACACTATAAATTGAGATACAGATGGAGATTTCTTGAGGCAGGGAGAGGAGGTCAGTGGCGGAGCTTGGAGATATCAGTAGCAGTGGAAGGGGCATGCAGAGGAGATTATATATGTTGATATGCTTCCTATGCTTCCTCTCTCCTCTGCCTGCCCCATCCACTCCTGCTGTTATCCCCTTCACGCGTCATACTGCGGATTAATCCCGTGTCCTCCTATATTTTTTTTCCAG";
   
  my $short_seq = $_[1];#"TGGAAGGGGCATGCAGAGGAGA";
  #print $short_seq;
  my $file = $short_seq;
  open (FH, "< $file") or die "Can't open $file for read: $!";
  my @lines;
  while (<FH>) {
    push (@lines, $_);
  }
  my $allcount = scalar @lines;
  print("All count:",$allcount);

  #chomp(@lines); #remove a newline from the end of every element in the array
  my $j=0;
  # foreach my $n (@lines){
  #   #print $n;
  #   ++$j;
  # }
  ### get match positions of short sequences and
  # construct postscript macro for RNAplot
  # in the format 
  # "start_pos end_pos 8 color omark" (8 is the thickness of the line)
  # start the macro by putting a circle around the
  # 5' end
  my $rnaplot_ps_macro = '1 cmark ' ;
  
  my $i = 0 ;
  foreach my $short_seq (@lines){
    #print $short_seq;
    my @positions = get_pos(\$long_seq,\$short_seq) ;
    print @positions;
    foreach (@positions) {
     my @shuffled = shuffle(@rgb_colors);
     # print "Shufled:",@shuffled;
     my $rgb_colors_val = $shuffled[0];
     # print "\nShuffled color:",$rgb_colors_val."\n";
     $rnaplot_ps_macro .= $_ .' 8 ' . convert_rgb_to_ps_format($rgb_colors_val) . ' omark ';
    }
    ++$i;
    # last if $i > 25;
  }
  print"\nTotal sRNAs:",$allcount," Matched:",$bcount,"\n";
  
  print "\n\nTranscript:",$_[2],"\nReady for RNAFold and RNAplot...";
  print "\nMacro to Plot: ",$rnaplot_ps_macro;   

### run RNAfold
  my $RNAfold_out_file = $_[2].'_RNAfold_out';#.'/data/RNAfold_out';
  print "\n\nRunning RNAFold for... ",$_[2];
  my $cmd = "echo '$long_seq' | $RNAfold_cmd > $RNAfold_out_file" ; # run RNAfold
  #print $cmd;
  #print "HERE.....\n";
  my $out = `$cmd` ;
  
  
  ### run RNAplot
  # RNAplot always generates its output in a file named after the input sequence
  # or (if not fasta) "rna.ps". Therefore we have to generate another directory
  my $RNAplot_out_file = 'rna.ps' ; # default RNAplot file name
  my $RNAplot_ps_final = $_[2].'_rna.ps_final' ; # file after addition of label for pdf 

  # make temp dir
  my $RNAplot_dir = $working_dir.'/'.$_[2].'_RNAplot_out';#.'/data/RNAplot_dir';
  
  #print $RNAplot_dir;
  mkdir $RNAplot_dir;
  move($RNAfold_out_file,$RNAplot_dir);
  #copy($RNAfold_out_file,$RNAplot_dir);
  chmod 0777, $RNAplot_dir ; 
  chdir $RNAplot_dir ; 

  # RNAplot command
  print "\nRunning RNAplot for... ",$_[2],"\n";
  $cmd = "$RNAplot_cmd --pre \"$rnaplot_ps_macro\" < $RNAfold_out_file " ; 
  #print $cmd;
  $out = `$cmd` ;
  
  $rnaplot_ps_macro = '' ; # reset to empty for the next macro coomand

  ### modify and convert images
  # We add a label to the postscript file, then convert it to pdf
  # and also generate the jpg file to print to the result page
  
  # add label to postscript file
  # open the RNAplot outputfile in its temp dir and another file for writing the modified ps
  my $label ="Secondary structure for_" .$_[2] ;

  open PS , '<', $RNAplot_dir .'/'. $RNAplot_out_file;
  open PS_FINAL, '>' , $RNAplot_dir . '/' . $RNAplot_ps_final;

  # add label to ps file
  # this is a bit of a hack:
  # the ps code is not inserted into the ideal position in the file
  # but by inserting just before "%%BeginProlog" we can print the label before the coordinate system is translated for plotting the RNA
  # post-script stack used:
  # Helvetica findfont -> change font
  # 12 scalefont -> set font size
  # 100 10 moveto -> set cursor position (bottom of page)
  # setfont
  # $label show -> print the label
  while (<PS>) {
    s/%%BeginProlog/\/Helvetica findfont\n12 scalefont\n100 10 moveto\nsetfont\n\($label\) show\n%%BeginProlog/ ;
    print PS_FINAL;
  }
  close PS ;
  close PS_FINAL ;

  #unlink $working_dir.'/'.$RNAplot_out_file; # remove rna.ps file


  # convert to PDF and write to results file
  my $pdf_file = $_[2].'Structure_plot.pdf';

  $cmd = $ps2pdf_cmd . ' '. $RNAplot_dir . '/' . $RNAplot_ps_final . ' ' . $pdf_file ;
  #print $cmd,"\n";
  chmod 0666, $pdf_file ; 
  my $output = `$cmd` ;
  if ($output) {
    print "There was a problem running ps2pdf" ;
  }

 #chdir $working_dir;  

# get match positions of short sequence on long sequence
sub get_pos{
  my ($long_seq_ref, $short_seq_ref) = @_ ;
  my $seq;
  $seq = $$short_seq_ref;
  #print $seq;
  # my $seq= $$short_seq_ref;
  chomp $seq;
  $seq=~ s/\R//g;
  # print $seq,"\n";

  print "\nMatching this sRNAs:",$seq,"\n";
  my @pos = () ;
  my $long_seq_len = length($$long_seq_ref) ;
  my $short_seq_len = length($seq) ;
  #print $short_seq_len;
  for (my $i = 1 ; $i <= $long_seq_len - $short_seq_len + 1; ++$i){#$long_seq_len - $short_seq_len + 1 ; ++$i) {
    my $sub=substr($$long_seq_ref, $i - 1, $short_seq_len);
    if (substr($$long_seq_ref, $i - 1, $short_seq_len) eq $seq) { # we have a match
    my $sub=substr($$long_seq_ref, $i - 1, $short_seq_len);
    print "Matched:",$seq,"\n";
    ++$bcount; 
    my $start = $i ;
    my $stop = $i + $short_seq_len - 1 ;
    push @pos, $start . ' ' . $stop ;
    }
  }
  
  return @pos ;
}

# Convert space delimited RGB into notation for RNAplot macros
# # where 255 = 1
 sub convert_rgb_to_ps_format {
   my $input_rgb = shift ; # in format R G B
   my ($r,$g,$b)=($input_rgb=~/(\d+)/g) ;
   $r = $r / 255 ;
   $g = $g / 255 ;
   $b = $b / 255 ;
   return "$r $g $b" ;
 }


} # end of subroutine run

#Functional Call

my $long_SEQ=$ARGV[0];
my $short_SEQ=$ARGV[1];
my $miR=$ARGV[2];

#print $long_SEQ;
#print $short_SEQ;
#print $miR;
run($long_SEQ,$short_SEQ,$miR);

