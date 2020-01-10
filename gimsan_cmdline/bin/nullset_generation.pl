#!/usr/bin/perl -w
use strict;
use File::Basename;

use lib "$ENV{GIMSAN}/bin/";
use Random;
#use Math::Random;

my $DEBUG0 = 1;
my @int2nt = ('A', 'C', 'G', 'T');
my $rand_obj;

main();

sub myrand
{
	#return Math::Random::random_uniform();
	return $rand_obj->random_uniform();
}

sub main 
{
	my $binsize = 5; # 5 percent
	my $wndsize = 100; # in base-pairs
	my $markov_order = 3;

	if(scalar(@ARGV) == 0) {
		print "usage: <program> --iters=<int> --input=<input.fsa>  --output-dir=<output folder> [--genome=<genome.fsa>] [--binsize=<int>] [--wndsize=<int>]\n";
		print "\n";
		print "Example: ./nullset_generation.pl --iters=10 --input=kluyveri.fsa --output-dir=testout_dir --genome=genome.fa\n";
		print "\n";
		print "Stand-alone version of sample_from_GC_content_range. Construct null datasets based on GC-content percentages\n";
		print "This is particularly created for the GIMSAN application.\n";
		print "\n";
		print "Default window size: $wndsize bp\n";
		print "Default GC-content bin size: $binsize%\n";
		print "If --genome is not specified, third-order Markov model is used to generate the nullset.";
		print "\n";
		exit(1);
	}

	my $iters = -1;
	my $genome_fn = '';
	my $input_fn = '';
	my $output_dir = '';
	my $type = '';

	my $coord_fn = '';
	foreach my $a (@ARGV) {
		if( $a =~ /^\--iters=(\d+)/) {
			$iters= $1;
		}
		elsif( $a =~ /^\--genome=(.+)/) {
			$genome_fn= $1;
		}
		elsif( $a =~ /^\--input=(.+)/) {
			$input_fn= $1;
		}
		elsif( $a =~ /^\--output-dir=(.+)/) {
			$output_dir= $1;
		}

		elsif( $a =~ /^\--binsize=(\d+)/) {
			$binsize = $1;
		}
		elsif( $a =~ /^\--wndsize=(\d+)/) {
			$wndsize= $1;
		}
		elsif( $a =~ /^\--coord=(\S+)/) {
			$coord_fn = $1;
		}
		else {
			die "Unrecognized parameter: $a\n";
		}
	}
	$rand_obj = Random->Random::new(10000000);
	chomp($genome_fn);
	chomp($input_fn);
	chomp($output_dir);

	if($genome_fn ne '') {
		$type = 'genome';
	}
	else {
		$type = 'markov';
	}

	if($iters <= 0) {
		die "Number of iterations must be specified.";
	}

	if(!-d $output_dir) {
		die "Directory $output_dir does not exist.";
	}

	#input file may contain non-ACGT
	my $inputfa = parseFa($input_fn);

	if($type eq 'genome') {
		nullset_generation_genome($iters, $inputfa, $genome_fn, $output_dir, $binsize, $wndsize);
	}
	elsif($type eq 'markov') {
		nullset_generation_markov($iters, $inputfa, $output_dir, $markov_order);
	}
	else {
		die "Invalid type.";
	}
}


sub nullset_generation_markov
{
	my ($iters, $inputfa, $output_dir, $markov_order) = @_;

	my $seq = '';
	foreach my $f (@$inputfa) {
		$seq .= $f->{seq};
	}

	my $transcount = estimate_markov_model($seq);
	my $model = normalize_model($transcount, $markov_order);
	#print_transprob(*STDERR, $model, $markov_order);
	my $generate_seqlens = precompute_generation_seqlens($inputfa);

	for(my $i = 1; $i <= $iters; $i++) {
		my %is_used_pos = (); 
		print "Generating null $i of $iters ... \n";

		my $out_fasta = sampleFromMarkovModel($inputfa, $generate_seqlens, $model, $markov_order, $i); 
		
		my $fn = sprintf("$output_dir/null-%05d.fa", $i);
		open( my $outfh, "> $fn") || die "Cannot open $fn for write: $!";
		writeFa($outfh, $out_fasta);
		close($outfh);
	}
}


sub nullset_generation_genome
{
	my ($iters, $inputfa, $genome_fn, $output_dir, $binsize, $wndsize) = @_;

	my $genome_raw= open_genome_file($genome_fn); #genome_raw only has ACGT
	my $genome = binGenomeBackground($genome_raw, $binsize, $wndsize);

	my $with_replace = 0;

	if($with_replace) {
		print "Sampling with replacement within each sequence-set\n";
	}
	else {
		print "Sampling without replacement within each sequence-set\n";
	}

	for(my $i = 1; $i <= $iters; $i++) {
		#use the same hash over all iterations to ensure all windows sampled over ALL iterations are non-overlapping
		#Here, we use a different hash, thus sequence-set A may overlap with sequence-set B, but windows in sequence-set A is non-overlapping.
		my %is_used_pos = (); 
		print "Generating null $i of $iters ... \n";

		my $out_fasta = sampleFromGenome($inputfa, $genome, $i, 0, \%is_used_pos); #sample without replacement
		
		my $fn = sprintf("$output_dir/null-%05d.fa", $i);
		open( my $outfh, "> $fn") || die "Cannot open $fn for write: $!";
		writeFa($outfh, $out_fasta);
		close($outfh);
	}

}
###########################################################################
# Markov sampling
###########################################################################

sub estimate_markov_model
{
	my ($seq) = @_;
	$seq = uc($seq);
	$seq =~ s/[^ACGTacgt]//g;

	my $pseudocount = 1;

	my @mat0 = ();
	my @mat1 = ();
	my @mat2 = ();
	my @mat3 = ();
	my @mat4 = ();
	my @mat5 = ();
	for(my $a = 0; $a < 4; $a++) {
		$mat0[$a] = $pseudocount;
		for(my $b = 0; $b < 4; $b++) {
			$mat1[$a][$b] = $pseudocount;
			for(my $c = 0; $c < 4; $c++) {
				$mat2[$a][$b][$c] = $pseudocount;
				for(my $d = 0; $d < 4; $d++) {
					$mat3[$a][$b][$c][$d] = $pseudocount;
					for(my $e = 0; $e < 4; $e++) {
						$mat4[$a][$b][$c][$d][$e] = $pseudocount;
						for(my $f = 0; $f < 4; $f++) {
							$mat5[$a][$b][$c][$d][$e][$f] = $pseudocount;
						}
					}
				}
			}
		}
	}

	$seq =~ tr/ACGT/0123/;
	my @arr = split(//, $seq);

	for(my $i = 0; $i < scalar(@arr); $i++) {
		$mat0[$arr[$i]]++;
	}
	for(my $i = 1; $i < scalar(@arr); $i++) {
		$mat1[$arr[$i-1]][$arr[$i]]++;
	}
	for(my $i = 2; $i < scalar(@arr); $i++) {
		$mat2[$arr[$i-2]][$arr[$i-1]][$arr[$i]]++;
	}
	for(my $i = 3; $i < scalar(@arr); $i++) {
		$mat3[$arr[$i-3]][$arr[$i-2]][$arr[$i-1]][$arr[$i]]++;
	}
	for(my $i = 4; $i < scalar(@arr); $i++) {
		$mat4[$arr[$i-4]][$arr[$i-3]][$arr[$i-2]][$arr[$i-1]][$arr[$i]]++;
	}
	for(my $i = 5; $i < scalar(@arr); $i++) {
		$mat5[$arr[$i-5]][$arr[$i-4]][$arr[$i-3]][$arr[$i-2]][$arr[$i-1]][$arr[$i]]++;
	}
	my %model = (
		mat0 => \@mat0,
		mat1 => \@mat1,
		mat2 => \@mat2,
		mat3 => \@mat3,
		mat4 => \@mat4,
		mat5 => \@mat5,
	);
	return \%model;
}

sub normalize_model
{
	my ($transcount, $markov_order) = @_;

	my @mat0 = ();
	my @mat1 = ();
	my @mat2 = ();
	my @mat3 = ();
	my @mat4 = ();
	my @mat5 = ();

	#make a copy
	my $offset = 1e-50;
	for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
		$mat0[$a] = $transcount->{mat0}->[$a] + $offset;
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			$mat1[$a][$b] = $transcount->{mat1}->[$a][$b] + $offset ;
			for(my $c = 0; $c < 4 && $markov_order >= 2; $c++) {
				$mat2[$a][$b][$c] = $transcount->{mat2}->[$a][$b][$c] + $offset;
				for(my $d = 0; $d < 4 && $markov_order >= 3; $d++) {
					$mat3[$a][$b][$c][$d] = $transcount->{mat3}->[$a][$b][$c][$d] + $offset;
					for(my $e = 0; $e < 4 && $markov_order >= 4; $e++) {
						$mat4[$a][$b][$c][$d][$e] = $transcount->{mat4}->[$a][$b][$c][$d][$e] + $offset;
						for(my $f = 0; $f < 4 && $markov_order >= 5; $f++) {
							$mat5[$a][$b][$c][$d][$e][$f] = $transcount->{mat5}->[$a][$b][$c][$d][$e][$f] + $offset;
						}
					}
				}
			}
		}
	}

	#normalize
	if(1) {
		my $sum = 0.0;
		for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
			$sum += $mat0[$a];
		}
		for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
			$mat0[$a] /= $sum;
		}
	}
	for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
		my $sum = 0.0;
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			$sum += $mat1[$a][$b];
		}
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			$mat1[$a][$b] /= $sum;
		}
	}
	for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			my $sum = 0.0;
			for(my $c = 0; $c < 4 && $markov_order >= 2; $c++) {
				$sum += $mat2[$a][$b][$c];
			}
			for(my $c = 0; $c < 4 && $markov_order >= 2; $c++) {
				$mat2[$a][$b][$c] /= $sum;
			}
		}
	}
	for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			for(my $c = 0; $c < 4 && $markov_order >= 2; $c++) {
				my $sum = 0.0;
				for(my $d = 0; $d < 4 && $markov_order >= 3; $d++) {
					$sum += $mat3[$a][$b][$c][$d];
				}
				for(my $d = 0; $d < 4 && $markov_order >= 3; $d++) {
					$mat3[$a][$b][$c][$d] /= $sum;
				}
			}
		}
	}
	for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			for(my $c = 0; $c < 4 && $markov_order >= 2; $c++) {
				for(my $d = 0; $d < 4 && $markov_order >= 3; $d++) {
					my $sum = 0.0;
					for(my $e = 0; $e < 4 && $markov_order >= 4; $e++) {
						$sum += $mat4[$a][$b][$c][$d][$e];
					}
					for(my $e = 0; $e < 4 && $markov_order >= 4; $e++) {
						$mat4[$a][$b][$c][$d][$e] /= $sum;
					}
				}
			}
		}
	}

	for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			for(my $c = 0; $c < 4 && $markov_order >= 2; $c++) {
				for(my $d = 0; $d < 4 && $markov_order >= 3; $d++) {
					for(my $e = 0; $e < 4 && $markov_order >= 4; $e++) {
						my $sum = 0.0;
						for(my $f = 0; $f < 4 && $markov_order >= 5; $f++) {
							$sum += $mat5[$a][$b][$c][$d][$e][$f];
						}
						for(my $f = 0; $f < 4 && $markov_order >= 5; $f++) {
							$mat5[$a][$b][$c][$d][$e][$f] /= $sum;
						}
					}
				}
			}
		}
	}
	my %model = (
		mat0 => \@mat0,
		mat1 => \@mat1,
		mat2 => \@mat2,
		mat3 => \@mat3,
		mat4 => \@mat4,
		mat5 => \@mat5,
	);
	return \%model;
}
sub print_transprob
{
	#my ($ofh, $transcount, $markov_order) = @_;
	#my $model = normalize_model($transcount, $markov_order);
	my ($ofh, $model, $markov_order) = @_;

	my $mat0 = $model->{mat0};
	my $mat1 = $model->{mat1};
	my $mat2 = $model->{mat2};
	my $mat3 = $model->{mat3};
	my $mat4 = $model->{mat4};
	my $mat5 = $model->{mat5};

	for(my $a = 0; $a < 4 && $markov_order >= 0; $a++) {
		printf $ofh ("%s %6.4lf ", $int2nt[$a], $mat0->[$a]);
		for(my $b = 0; $b < 4 && $markov_order >= 1; $b++) {
			printf $ofh ("%s%s %6.4lf ", $int2nt[$b], $int2nt[$a], $mat1->[$b][$a]);
			for(my $c = 0; $c < 4 && $markov_order >= 2; $c++) {
				printf $ofh ("%s%s%s %6.4lf ", $int2nt[$c], $int2nt[$b], $int2nt[$a], $mat2->[$c][$b][$a]);
				for(my $d = 0; $d < 4 && $markov_order >= 3; $d++) {
					printf $ofh ("%s%s%s%s %6.5lf ", $int2nt[$d], $int2nt[$c], $int2nt[$b], $int2nt[$a], $mat3->[$d][$c][$b][$a]);
					for(my $e = 0; $e < 4 && $markov_order >= 4; $e++) {
						printf $ofh ("%s%s%s%s%s %6.5lf ", $int2nt[$e], $int2nt[$d], $int2nt[$c], $int2nt[$b], $int2nt[$a], $mat4->[$e][$d][$c][$b][$a]);
						for(my $f = 0; $f < 4 && $markov_order >= 5; $f++) {
							printf $ofh ("%s%s%s%s%s%s %6.5lf ", $int2nt[$f], $int2nt[$e], $int2nt[$d], $int2nt[$c], $int2nt[$b], $int2nt[$a], $mat5->[$f][$e][$d][$c][$b][$a]);
						}
					}
				}
			}
		}
		print $ofh "\n";
	}
}




sub sampleFromMarkovModel
{
	my ($inputfsa, $sample_seqlens, $model, $markov_order, $iter_num) = @_; 

	my $total_len = 0;
	foreach my $s (@$sample_seqlens) {
		$total_len += $s;
	}

	if($markov_order != 3) {
		die "Only implemented for Markov-order 3.";
	}

	my $seqall = sampleSeqFromMarkovModel($model, $total_len);
	$seqall =~ tr/0123/ACGT/;
	my $start_index = 0;
	my @out_fasta = ();
	for(my $i = 0; $i < scalar(@$sample_seqlens); $i++) {
		my $len = $sample_seqlens->[$i];
		my $seq = substr($seqall, $start_index, $len);

		my $newseq = insertCorrespondingGaps($inputfsa->[$i]->{seq}, $seq);
		my $newtag = sprintf("id=%04d-%04d tag=[%s] markov-order=$markov_order", $iter_num, $i, $inputfsa->[$i]->{tag});

		my %struct = (
			seq => $newseq,
			tag => $newtag,
		);
		push(@out_fasta, \%struct);

		$start_index += $len;
	}
	return \@out_fasta;
}

sub sampleNucleotide 
{
	my $prob = shift;

	my $sum = 0.0;
	my $rand = myrand();
	for(my $j = 0; $j < 4; $j++) {
		$sum += $prob->[$j];
		if($sum > $rand) {
			return $j;
		}
	}
	return 3;
}

sub sampleSeqFromMarkovModel
{
	my ($model, $seqlen) = @_;

	my @seq = ();
	$seq[0] = sampleNucleotide( $model->{mat0} );
	$seq[1] = sampleNucleotide( $model->{mat1}->[$seq[0]] );
	$seq[2] = sampleNucleotide( $model->{mat2}->[$seq[0]][$seq[1]] );
	for(my $i = 3; $i < $seqlen; $i++) {
		my $prob = ( $model->{mat3}->[$seq[$i-3]][$seq[$i-2]][$seq[$i-1]] );
		my $sum = 0.0;
		my $rand = myrand();
		my $nt = 3;
		for(my $j = 0; $j < 4; $j++) {
			$sum += $prob->[$j];
			if($sum > $rand) {
				$nt = $j;
				last;
			}
		}
		$seq[$i] = $nt;
	}
	return (join('', @seq));;
}


sub precompute_generation_seqlens
{
	my ($inputfsa) = @_;

	my @seqlen = ();

	for(my $i = 0; $i < scalar(@$inputfsa);$i++) {
		my $seq = $inputfsa->[$i]->{seq};
		$seq =~ s/[^ACGTacgt]//g;
		$seqlen[$i] = length($seq);
	}
	return \@seqlen;
}



###########################################################################
# GC-content
###########################################################################

sub compute_GC_content
{
	my $seq = shift;

	#$seq =~ s/\s//g;
	$seq =~ s/[^ACGTacgt]//g; #modified so that GC + AT content is 100%
	(my $gc_seq = $seq) =~ s/[^GCgc]//g; #remove all non-GC characters

	if(length($seq) == 0) {
		return 0.0;
	}
	else {
		return length($gc_seq) / length($seq) * 100; #convert to percentage
	}
}

sub revcompl
{
	my $seq = shift;
	$seq = reverse($seq);
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return $seq;
}

###########################################################################
# Sampling chunks
###########################################################################
sub sampleFromGenome
{
	my $inputfa = shift;
	my $genome = shift;
	my $iter_num = shift;
	my $with_replace = shift; #0 for without, 1 for with
	my $usedPos = shift;

	if($with_replace != 1 && $with_replace != 0) {
		die "Error: the variable with_replace has value $with_replace.\n";
	}
	
	my @fasta = (); #output
	my $i = 0;
	foreach my $f (@$inputfa) {
		my $tag = $f->{tag};
		my $seq = $f->{seq};

		(my $clean_seq = $seq) =~ s/[^ACGTacgt]//g; #filter out all invalid characters when sampling from genome

		#sample and insert appropriate gaps according to the original seq
		my $sampled = sampleFromSeq($clean_seq, $genome, $usedPos, $with_replace);
		my $newseq;
		if(length($sampled) == length($seq)) {
			$newseq = $sampled;
		}
		else {
			#print STDERR "Inserting corresponding gaps...\n";
			$newseq = insertCorrespondingGaps($seq, $sampled);
		}

		my $old_GC_content = compute_GC_content($clean_seq);
		my $new_GC_content = compute_GC_content($sampled);

		my $genomeName = $genome->{genomeName};
		my $newtag = sprintf("id=%04d-%04d tag=[$tag] original_GC=%.1f new_GC=%.1f binsize=%d wndsize=%d sampled_genome=$genomeName", 
			$iter_num, $i, $old_GC_content, $new_GC_content, $genome->{binsize}, $genome->{wndsize}
		);

		#printf ("Original seqlen: %d, new seqlen: %d \n", length($seq) ,length($newseq));
		
		my %struct = (
			seq => $newseq,
			tag => $newtag,
		);
		push(@fasta, \%struct);
		$i++;
		
	}
	return \@fasta;
}

# Given a clean sequence, it samples from a genomic background (null distribution) based on a specified
# window size of GC-content.
sub sampleFromSeq
{
	my $seq = shift;
	my $genome = shift;
	my $usedPos = shift;
	my $with_replace = shift;

	my $newseq = '';
	my $i = 0; 
	while($i < length($seq)) {
		# Take cares of the right-end if smaller than window size
		my $wndsize = ( length($seq)-$i < $genome->{wndsize} ? length($seq)-$i : $genome->{wndsize});

		my $gc_content = compute_GC_content(substr($seq, $i, $wndsize));
		my $bin_ind = getBinInd($genome, $gc_content); 
		my $numbins = scalar(@{$genome->{bins}});

		my $current_bin_attempts= 0;
		my $reget_bin_attempts = 0;

		my $pos;
		do {
			if(++$current_bin_attempts > 10) {
				print STDERR "Too many attempts in sampling from bin $bin_ind: $current_bin_attempts. ";
				#move toward larger bin
				$bin_ind = regetLargerBinInd($genome, $bin_ind);
				$current_bin_attempts = 0;
				print STDERR "Attempting to use bin $bin_ind...\n";

				if(++$reget_bin_attempts > $numbins * 100) {
					#to avoid non-halting
					die "Error: Too many attempts to get new bin: $reget_bin_attempts\n";
				}
			}
			#sample a position
			my $bin = $genome->{bins}->[$bin_ind];
			$pos = $bin->[ int(myrand()*scalar(@$bin)) ];
		}while(!$with_replace && isUsedPos($pos, $usedPos, $wndsize));
		
		if(!defined($pos)) {
			die "pos is not defined";
		}

		markUsedPos($pos, $usedPos, $wndsize);

		$newseq .= substr($genome->{seq}, $pos, $wndsize);
		$i += $wndsize;
	}

	if(length($newseq) != length($seq)) {
		die "Sequence lengths do not match";
	}
	return $newseq;
}

sub isUsedPos 
{
	my $pos = shift;
	my $usedPos = shift;
	my $wndsize = shift;

	#Hack to make the check a bit faster since "used" occurs mostly on the ends
	if(exists( $usedPos->{$pos + $wndsize - 1} )) {
		return 1;
	}

	for(my $i = 0; $i < $wndsize; $i++) {
		if(exists( $usedPos->{$i+$pos} )) {
			return 1;
		}
	}
	return 0;
}

sub markUsedPos
{
	my $pos = shift;
	my $usedPos = shift;
	my $wndsize = shift;

	for(my $i = 0; $i < $wndsize; $i++) {
		$usedPos->{$i + $pos} = 1;
	}
}

# Insert gaps (e.g. N, X, -, etc) into an ungapped sequence according to a gapped sequence
sub insertCorrespondingGaps
{
	my $gapped_seq = shift;
	my $ungapped_seq = shift;

	(my $tmp_seq = $gapped_seq) =~ s/[^ACGTacgt]//g;
	if($DEBUG0){
		if(length($tmp_seq) != length($ungapped_seq)) {
			die "The number of non-gap characters in the two sequences are inconsistent.";
		}
	}

	my @gapped = split(//, $gapped_seq);
	my @ungapped = split(//, $ungapped_seq);
	my @newseq = ();

	my $i = 0;
	for(my $g = 0; $g < scalar(@gapped); $g++) {
		#$g is the index for gapped, and $i is for ungapped
		if($gapped[$g] =~ /[ACGTacgt]/) { #if it is a non-gap character
			$newseq[$g] = $ungapped[$i];	
			$i++;
		}
		else {
			$newseq[$g] = $gapped[$g];
		}
	}

	if($DEBUG0) {
		if( scalar(@newseq)!= scalar(@gapped)) {
			die "Error: new gapped sequence has different length with the original gapped sequence";
		}
		if( $i != scalar(@ungapped) ) {
			die "Error: not all the characters of ungapped sequence was used";
		}
	}


	my $str = join('', @newseq);
	return $str;
}
###########################################################################
# Hash bins
###########################################################################
sub binGenomeBackground
{
	my $genome_raw = shift;
	my $binsize = shift;
	my $wndsize = shift;

	my $seq = $genome_raw->{seq};

	my $bins = createHashBins($binsize);

	my %genome = (
		genomeFn => $genome_raw->{filename},
		genomeName => $genome_raw->{name},
		seq => $seq,
		bins => $bins,
		binsize => $binsize, #percentage of GC-content
		wndsize => $wndsize, #number of bp 
	);

	print "Start binning... ";
	my $numbins = scalar(@$bins);
	my $startpos = int( myrand() * $wndsize); #allow more variations between sets
	for(my $i = $startpos; $i < length($seq) - $wndsize + 1; $i+=int($wndsize/2)) {
		my $str = substr($seq, $i, $wndsize);
		my $content = compute_GC_content($str); 

		# $content is not exactly the key because it needs to be binned. $i is the value
		my $bin_ind = int($content / $binsize);
		push(@{$bins->[$bin_ind]}, $i);
	}
	print "Done!\n";

	printHashBins(\%genome);

	return \%genome;
}

sub createHashBins
{
	my $binsize = shift;
	
	my $numbins = int( 100 / $binsize) + 1;

	my @bins = ();
	for(my $i = 0; $i < $numbins; $i++) {
		my @lst = ();
		$bins[$i] = \@lst;
	}
	return \@bins;
	
}

sub printHashBins
{
	my $hb = shift;

	printf ("Bin size: %d\n", $hb->{binsize});
	printf ("Number of bins: %d\n", scalar(@{$hb->{bins}}));
	for(my $i = 0; $i < scalar(@{$hb->{bins}}); $i++) {
		printf("Bin %2d [%3.1f%%, %3.1f%%): %d\n", $i, $i * $hb->{binsize}, ($i+1) * $hb->{binsize}, scalar(@{$hb->{bins}->[$i]}));
	}
	#print "\n";
}


#Given a bin index, return a bin index plus/minus 1 (whichever one is a larger bin than the input).
sub regetLargerBinInd
{
	my $genome = shift;
	my $bin_ind = shift;

	my $attempt = 0;
	do{
		if( $bin_ind + 1 >= scalar(@{$genome->{bins}}) ) {
			$bin_ind--; # bin_ind+1 is out of range
		}
		elsif( $bin_ind - 1 < 0 ) {
			$bin_ind++; # bin_ind-1 is out of range
		}
		else {
			my $left = scalar(@{$genome->{bins}->[$bin_ind-1]});
			my $right = scalar(@{$genome->{bins}->[$bin_ind+1]});

			#add pseudocounts
			$left += 100;
			$right += 100;

			my $left_prob = $left / ($left + $right);

			if(myrand() < $left_prob) { 
				$bin_ind--;
			}
			else {
				$bin_ind++;
			}
		}
		if(($attempt++) > 500) {
			die "Fail to find a non-empty bin in getBinInd()";
		}
	} while( scalar(@{$genome->{bins}->[$bin_ind]}) < 100 ); #better lower bound 2/28/2008
	return $bin_ind;
}

# Get the best possible bin index (the bin should not be empty)
# Heuristic: the number of (possibly overlapping) positions in the bin should be >= a constant.
#            If not, then iteratively chooses a larger bin.
sub getBinInd
{
	my $genome = shift;
	my $gc_content = shift;

	my $bins = $genome->{bins};
	my $numbins = scalar(@$bins);

	#Get the exact bin index
	my $bin_ind = int($gc_content / $genome->{binsize});
	if($bin_ind < 0 || $bin_ind >= $numbins) {
		die "Index is incorrect \n";
	}

	if( scalar(@{$genome->{bins}->[$bin_ind]}) < 100 ){
		$bin_ind = regetLargerBinInd($genome, $bin_ind);
	}
	return $bin_ind;
}

###########################################################################
# FASTA
###########################################################################


sub open_genome_file
{
	my ($fn) = @_;

	open(IN, "< $fn") or die "Cannot open $fn for read\n";

	my @new_in = ();
	foreach my $s (<IN>) {
		chomp $s;
		$s =~ s/[\r\n]//g;
		if($s !~ /^\>/) { #if not tag
			push @new_in, $s;
		}
	}
	close(IN);

	my $in_join = join("", @new_in);
	$in_join =~ s/[^ACGTacgt]//g; #filter out all invalid characters

	my $name = File::Basename::basename($fn);

	my %hash = (
		seq => $in_join,
		name => $name,
		filename => $fn,
	);

	return \%hash;
}

sub parseFa
{
	my $fname = shift;
	open(FILE, "< $fname") or die "Couldn't open file $fname: $!\n";
	my @istream = <FILE>;
	close(FILE);

	my $title = undef;
	my $seq = undef;
	my @lst = ();
	foreach my $s (@istream) {
		chomp $s;
		$s =~ s/[\r\n]//g;

		if(length($s) == 0) {
			next;
		}
		elsif($s =~ /^\>(.+)/) {
			my $newTitle = $1;
			#if($seq ne '') {
			#	push(@lst, parseFa_createStruct($seq, $title));
			#	$seq = '';
			#}

			#allow for empty sequences
			if(defined($seq) && defined($title)) {
				push(@lst, parseFa_createStruct($seq, $title));
			}

			$title = $newTitle;
			$seq = '';
		}
		else {
			$s =~ s/\s//g;
			$seq = $seq . $s;
		}

	}
	#allow for empty sequences
	if(defined($seq) && defined($title)) {
		push(@lst, parseFa_createStruct($seq, $title));
	}
	#if($seq =~ /[ACGTacgt]/) {
	#	push(@lst, parseFa_createStruct($seq, $title));
	#}
	printf ("$fname has %d sequences.\n", scalar(@lst));
	return \@lst;
}

sub parseFa_createStruct
{
	my $seq = shift;
	my $title = shift;

	my %hash = (
		tag => $title,
		seq => $seq,
	);
	return \%hash;
}

sub writeFa
{
	my $fh = shift;
	my $fasta = shift;

	foreach my $f (@$fasta) {
		printf $fh (">%s\n", $f->{tag});
		printf $fh ("%s\n", $f->{seq});
	}
}


