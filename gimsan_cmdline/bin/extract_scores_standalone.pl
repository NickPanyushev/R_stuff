#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd qw(abs_path);

my $R_source = "$ENV{GIMSAN}/bin/conf_pval_only.R";

my $DEBUG = 0;

main();

sub main 
{
	if(@ARGV == 0) {
		print "Usage: ./program.pl --ref-result=<filename> --nullset-dir=<directory> --nullscores-out=<filename> --R-out=<filename>\n";
		print "\n";
		print "This program extracts scores for the statistical significance analysis of a particular motif-finder and writes to a *.R file.\n";
		print "This is made as a standalone version of extract_scores.pl for GIMSAN\n";
		print "\n";
		print "Input: --ref-result and --nullset-dir\n";
		print "Output: --nullscores-out and --R-out \n";
		print "\n";
		print <<"End";
Inputs:
--ref-result is the filename of the results from the reference input FASTA file
--nullset-dir is the directory that contains the results from the nullset. The filenames must be in the form null-00001.stdout, null-00002.stdout, etc.

Outputs:
--nullscores-out is the file where the scores from the nullset is written to.
--R-out is the R-script that gets generated
End
		exit;
	}

	my $max_rank = 1;
	my $finder_result_fn = ''; 
	my $nullset_dir = ''; 
	my $Rfilename = '';
	my $scores_fn = '';

	foreach my $a (@ARGV) {
		#if($a =~ /--max-rank=(\d+)/) {
		#	$max_rank = $1;
		#}
		if($a =~ /--ref-result=(\S+)/) {
			$finder_result_fn = $1;
		}
		elsif($a =~ /--nullset-dir=(\S+)/) {
			$nullset_dir = $1;
		}
		elsif($a =~ /--R-out=(\S+)/) {
			$Rfilename = $1;
		}
		elsif($a =~ /--nullscores-out=(\S+)/) {
			$scores_fn = $1;
		}
		
	}
	
	my $ref_score = extract_score_by_rank($finder_result_fn, 1);

	my @null_scores = ();
	my $count = 0;

	#while(-f (my $fn = sprintf("$nullset_dir/gibbsmarkov01-%05d.stdout", $count+1))){
	while(-f (my $fn = sprintf("$nullset_dir/null-%05d.stdout", $count+1))){
		my $score = extract_score_by_rank($fn, 1);
		push(@null_scores, $score);
		$count++;
	}

	if($count != scalar(@null_scores) ) {
		die "Error: Inconsistent number of scores found.";
	}
	if( $count == 0) {
		print STDERR ("Nullset-dir: $nullset_dir\n");
		die "Error: zero null scores found.";
	}

	print "Score of reference motif: $ref_score\n";
	printf ("Size of null set: %d\n", $count);
	
	open(my $scores_ostream, ">$scores_fn") or die "Cannot open $scores_fn for write: $!"; 
	print $scores_ostream join("\n", @null_scores);
	close($scores_ostream);

	#R-code generation
	open(my $Rfh, ">$Rfilename") or die "Cannot open $Rfilename for write: $!"; 
	print $Rfh ("#===================================================================\n\n");
	printf $Rfh ("# Reference filename: %s\n", basename($finder_result_fn));
	print $Rfh ("#\n");
	$R_source = abs_path($R_source);
	print $Rfh ("source(\"$R_source\")" . "\n");
	print $Rfh ("library(MASS)\n");
	print $Rfh ("\n");
	$scores_fn =~ s/[\r\n]//g;
	$scores_fn = abs_path($scores_fn);
	printf $Rfh ("sample<-scan(\"%s\")\n", $scores_fn);
	printf $Rfh ("getConfPvalLat(%s, sample, conf=0.1, mins=7, maxs=200)\n", $ref_score);
	close($Rfh);
}

sub extract_score_by_rank
{
	my ($fn, $rank) = @_;

	open(IN, "< $fn") or die "Cannot open $fn for read.";
	my $score = undef;
	my $start_record;
	if($rank == 1) {
		$start_record = 1; #since there are no "RANK XX" tag
	}
	else {
		$start_record = 0; 
	}
	foreach my $ln (<IN>) {
		if($ln =~ /RANK\s+$rank\s+/) {
			$start_record = 1;
		}
		if($start_record) {
			if($ln =~ /Log Markovian-CLR generated from sites without pseudocount:\s*(\-?[\d\.]+)/) {
				$score = $1;
				last;
			}
		}
	}
	if($rank == 1 && !defined($score)) {
		die "Error: Rank-1 MCLR score not found in $fn\n";
	}
	close(IN);
	return $score;
}



