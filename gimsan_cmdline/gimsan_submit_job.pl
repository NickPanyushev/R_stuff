#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX qw(ceil floor);
use subs 'abs_path';
use Cwd;
use File::Copy;
use File::Basename;
use File::Compare;
use Time::localtime;

# GIMSAN evaluates the statistical significance of a particular motif-finder. This script qsub 
# jobs sampled from a specified null distribution. 

my $DEBUG = 0;
if(@ARGV < 1) {
	print_help();
	exit;
}

if(!defined($ENV{GIMSAN})) {
	print "Error: GIMSAN environment variable is undefined.\n";
	print "\n";
	gimsan_env_help();
	exit;
}
my $NULLSET_GENERATOR = abs_path("$ENV{GIMSAN}/bin/nullset_generation.pl");
my $GM_EXE = abs_path("$ENV{GIMSAN}/gibbsmarkov/gibbsmarkov.out");
my $CWD = getcwd();
main();

#Overriding Cwd::abs_path
sub abs_path
{
	my $str = shift;
	$str = resolve_tilde($str);
	if(-e $str) {
		return Cwd::abs_path($str);
	}
	else {
		return $str;
	}
}



sub main 
{
	if(!-f $NULLSET_GENERATOR) {
		print "Null set generator cannot be found at:\n $NULLSET_GENERATOR\n\n";
		print "Please check GIMSAN environment variable.\n";
		print "\n";
		gimsan_env_help();
		exit(1);
	}
	if(!-f $GM_EXE) {
		print "Motif-finder GibbsMarkov cannot be found at:\n $GM_EXE\n\n";
		print "Please check GIMSAN environment variable and verify that GibbsMarkov has been compiled correctly.\n";
		print "\n";
		gimsan_env_help();
		exit(1);
	}

	my $spec = parse_params(\@ARGV); 	

	if($spec->{is_batch_mode}) {

		main_batch_mode($spec);
	}
	else {
		main_iter_mode($spec);
	}
}

#############################################################################################
# Parameters 
#############################################################################################
sub print_help
{ 
	print "DESCRIPTION\n";
	print "  Submits a GIMSAN job on your current machine or on a cluster via qsub\n";
	print "\n";
	print "USAGE\n";
	print "  perl gimsan_submit_job.pl --f=<FASTA> --w=<int-list> [OPTIONS]\n";
	print "\n";
	print "EXAMPLES\n";
	print "  Single machine job (dual processor) example:\n";
	print "       ./gimsan_submit_job.pl --f=input.fsa --w=8,14,20 --bg=cer_intergenic.fsa\n";
	print "\n";
	print "  QSUB job (allocating at most 3 machines) example: \n";
	print "       ./gimsan_submit_job.pl --mach=3 --f=input.fsa --w=8,14,20 --bg=cer_intergenic.fsa\n";
	print "\n";
	print "REQUIRED\n";
	print "  --f=<FASTA>         input FASTA file for motif-finding\n";
	print "  --w=<int-list>      motif width (e.g. --w=8,14,20)\n";
	print "\n";
	print "OPTIONAL\n";
	print "  --bg=<FASTA>        file (FASTA format) for background model estimation [default: file specified in --f]\n";
	print "  --nullset=<int>     size of null set [default: 20]\n";
	print "\n";
	print "  --proc=<int>        number of processes to run simultaneously on a single machine [default: 2]\n";
	print "  --mach=<int>        number of machines to allocate using QSUB [default: 1]\n";
	print "  --q=<str-list>      a list of queue names (e.g. --q=general.q,long.q)\n";
	print "  --N=<str>           name of experiment [default: file specified in --f]\n";
	print "  --outdir=<str>      output directory [default: current working directory]\n";
	print "\n";
	print "  --L=<int>           rapid convergence rate [default: 200]\n";
	print "  --ds=<yes|no>       consider double strand [default: yes]\n";
	print "  --markov=<int>      order of Markov background (0 to 5) [default: 5]\n";
	print "\n";
	print "  Occurrence per sequence (choose one of the following: --zoops or --oops) [default: --zoops=0.2]\n";
	print "  --zoops=<flt>       weight of ZOOPS (Zero or One Occurrence Per Sequence) mode\n";
	print "  --oops              OOPS (One Occurrence Per Sequence) mode\n";
	print "\n";
	print "  Stopping criterion (choose one of the following: --cput or --t) [default: --cput=300]\n";
	print "  --cput=<int>        cpu-time of each finder execution\n";
	print "  --t=<int>           number of cycles (random restarts)\n";
	print "\n";
	print "URL\n";
	print "  Visit GIMSAN homepage at http://www.cs.cornell.edu/~ppn3/gimsan/\n";
}

sub parse_params
{
	my ($ARGV) = @_;

	my %spec_hash = ();
	my $spec = \%spec_hash;

	#defaults
	$spec->{qsub_jobname} = undef;
	$spec->{iterind} = undef;
	$spec->{iters_slicer} = undef;
	$spec->{fasta_fn} = undef;
	$spec->{bgfile_fn} = undef;
	$spec->{queue_names} = undef;
	$spec->{destdir} = undef;
	$spec->{tmpdir} = undef;
	$spec->{nullset_dir} = undef;
	$spec->{iters_slicer} = undef;
	$spec->{iter_mode_param} = undef;
	$spec->{expname} = undef;
	$spec->{is_batch_mode} = 1;
	$spec->{submit} = 1;
	$spec->{proc_per_machine} = 2;
	$spec->{nullset_size} = 20;
	$spec->{machines_alloc} = 1; #number of machines to allocate
	$spec->{per_seq_mode} = 'zoops';
	$spec->{ds} = 1; #double-strand mode
	$spec->{zoops_weight} = 0.2;
	$spec->{stop_crit} = 'cput';
	$spec->{stop_crit_int} = 300;
	$spec->{markov_order} = 5;
	$spec->{rapid_conv} = 200;

	my @widths = ();
	my $outdir = $CWD;

	#parsing ARGV
	foreach my $a (@$ARGV) {
		$a =~ s/[\r\n]//g; #filter-out newline chars

		if($a =~ /^--iter=(\d+)$/) { #internal use only
			$spec->{is_batch_mode} = 0;
			$spec->{iterind} = $1;
		}
		elsif($a =~ /^--f=(\S+)$/)  {
			$spec->{fasta_fn} = abs_path($1);
		}
		elsif($a =~ /^--bg=(\S+)$/) {
			$spec->{bgfile_fn} = abs_path($1);
		}
		elsif($a =~ /^--w=([\d\,]+)$/) {
			my $str = $1;
			@widths = split(/,/, $str);
		}
		elsif($a =~ /^--nullset=(\d+)$/) {
			$spec->{nullset_size} = $1;
		}
		elsif($a =~ /^--oops$/) {
			$spec->{per_seq_mode} = 'oops';
		}
		elsif($a =~ /^--zoops=(\d*\.\d+)$/) {
			$spec->{per_seq_mode} = 'zoops';
			$spec->{zoops_weight} = $1;
		}
		elsif($a =~ /^--cput=(\d+)$/) {
			$spec->{stop_crit} = 'cput';
			$spec->{stop_crit_int} = $1;
		}
		elsif($a =~ /^--t=(\d+)$/) {
			$spec->{stop_crit} = 't';
			$spec->{stop_crit_int} = $1;
		}
		elsif($a =~ /^--L=(\d+)$/) {
			$spec->{rapid_conv} = $1;
		}
		elsif($a =~ /^--ds=(\S+)$/) {
			if($1 eq 'yes') {
				$spec->{ds} = 1;
			}
			elsif($1 eq 'no') {
				$spec->{ds} = 0;
			}
			else {
				die "--ds must be either 'yes' or 'no'";
			}
		}
		elsif($a =~ /^--markov=(\d)$/) {
			$spec->{markov_order} = $1;
		}
		elsif($a =~ /^--N=(\S+)$/) {
			$spec->{expname} = $1;
		}
		elsif($a =~ /^--outdir=(\S+)$/) {
			$outdir = abs_path($1);
		}
		elsif($a =~ /^--proc=(\d+)$/) {
			$spec->{proc_per_machine} = $1; 
		}
		elsif($a =~ /^--mach=(\d+)$/) {
			$spec->{machines_alloc} = $1; #needed to compute 'iters_slicer'
		}
		elsif($a =~ /^--q=(\S+)$/) {
			$spec->{queue_names} = $1;
		}
		
		else {
			print "ERROR: Invalid parameter $a in gimsan_submit_job\n\n";
			print_help();
			exit(1);
		}
	}

	#bound-check
	if(! (-d $outdir)) {
		print STDERR "The directory $outdir does not exist!\n";
		exit;
	}
	if(!defined($spec->{fasta_fn}) || scalar(@widths) == 0) {
		print "ERROR: --f and --w must be specified\n";
		exit(1);
	}
	if(!check_fasta($spec->{fasta_fn})) {
		print STDERR "Invalid FASTA file " . $spec->{fasta_fn} . "\n";
		exit;
	}
	if(!check_fasta($spec->{bgfile_fn})) {
		print STDERR "Invalid FASTA file " . $spec->{bgfile_fn} . "\n";
		exit;
	}
	foreach my $w( @widths) {
		if($w < 5) { # if there is non-digit
			print STDERR "ERROR: Motif width must be at least 5\n";
			exit;
		}
	}
	if($spec->{nullset_size} < 0) {
		print STDERR "ERROR: Null set size must be at least 0.\n";
		exit;
	}
	foreach my $a (@$ARGV) {
		my $has_seen = 0;
		if($a =~ /--oops/ || $a=~ /--zoops/) {
			if($has_seen) {
				print STDERR "ERROR: --oops and --zoops can not be used together\n";
				exit;
			}
			$has_seen = 1;
		}
	}
	foreach my $a (@$ARGV) {
		my $has_seen = 0;
		if($a =~ /--cput/ || $a=~ /--t/) {
			if($has_seen) {
				print STDERR "ERROR: --cput and --t can not be used together\n";
				exit;
			}
			$has_seen = 1;
		}
	}
	if($spec->{stop_crit_int} < 1) {
		print STDERR "ERROR: specified --cput or --t must be at least 1.\n";
		exit;
	}
	if($spec->{rapid_conv} < 1) {
		print STDERR "ERROR: rapid convergence rate --L must be at least 1.\n";
		exit;
	}
	if($spec->{markov_order} < 0 || $spec->{markov_order} > 5) {
		print STDERR "ERROR: Markov-order must be between 0 and 5.\n";
		exit;
	}
	if($spec->{proc_per_machine} < 1 ) {
		print STDERR "ERROR: --proc must be at least 1\n";
		exit;
	}
	if($spec->{machines_alloc} < 1 ) {
		print STDERR "ERROR: --mach must be at least 1\n";
		exit;
	}
	
	#actions 
	if(!defined($spec->{expname}) || $spec->{expname} eq '') {
		#my ($DAY, $MONTH, $YEAR) = localtime[3,4,5];
		#print STDERR "year $YEAR";
		my $date = sprintf("%04d%02d%02d", localtime->year() + 1900, localtime->mon()+1, localtime->mday());
		$spec->{expname} = sprintf("%s_%s",$date,basename($spec->{fasta_fn}, ('.fa', '.fsa')));
	}
	if(!defined($spec->{qsub_jobname}) || $spec->{qsub_jobname} eq '') {
		$spec->{qsub_jobname} = 'gimsan_' . $spec->{expname};
	}
	if(defined($spec->{bgfile_fn}) && compare($spec->{fasta_fn},$spec->{bgfile_fn}) == 0) {
		print STDERR "Warning: --f and --bg specified the same file. This is the same as not specifying a --bg file\n";
		#if the user specifies the same file, then statsig should be estimated with Markov model rather than mosaic sampling
		$spec->{bgfile_fn} = undef; 
	}
	if(!defined($spec->{bgfile_fn})) {
		print STDERR "Warning: specifying an intergenic genome file --bg is highly RECOMMENDED\n";
	}
	$spec->{destdir} = sprintf("%s/%s/", $outdir, $spec->{expname});
	$spec->{tmpdir} = sprintf("%s/tmp/", $spec->{destdir});
	$spec->{nullset_dir} = sprintf("%s/nullset/", $spec->{destdir});
	$spec->{statsig_dir} = sprintf("%s/statsig/", $spec->{destdir});
	$spec->{iters_slicer} = ceil( ($spec->{nullset_size}+1) / $spec->{machines_alloc} );
	if($spec->{is_batch_mode}) {
		$spec->{iter_mode_param} = join(' ', @$ARGV);
		if($spec->{iter_mode_param} =~ /--iter=/) {
			die "ERROR: --iter found when under batch-mode";
		}
		if($spec->{iter_mode_param} !~ /--N=/) {
			$spec->{iter_mode_param} .= sprintf(" --N=%s ", $spec->{expname});
		}
	}
	else {
		$spec->{iter_mode_param} = undef;
	}

	#generate finder parameters
	my @finders = ();
	my $ds_str = ($spec->{ds} ? '-ds' : '');
	my $ops_str = '';
	if($spec->{per_seq_mode} eq 'oops') {
		$ops_str = '';
	}
	elsif($spec->{per_seq_mode} eq 'zoops') {
		$ops_str = sprintf("-zoops %f", $spec->{zoops_weight});
	}
	else {
		die "ERROR: invalid per-sequence-mode - it must be OOPS or ZOOPS";
	}
	my $runtime_str = sprintf("-%s %d -L %d", $spec->{stop_crit}, $spec->{stop_crit_int}, $spec->{rapid_conv});
	my $markov_str = sprintf("-markov %d", $spec->{markov_order});
	my $bgfile_str = '';
	if(defined($spec->{bgfile_fn})) {
		$bgfile_str = sprintf("-bfile %s", $spec->{bgfile_fn});
	}
	my $param_base = sprintf("-gibbsamp -best_clr $runtime_str $ds_str $ops_str $bgfile_str $markov_str -p 0.05 -em 0 ");
	foreach my $w (@widths) {
		if($w =~ /^\d+$/) {
			$w = int($w);
			my $name = sprintf("width%03d", $w);
			my %hash = (
				name => $name,
				path => $GM_EXE,
				param => "-l $w $param_base",
			);
			push(@finders, \%hash);
		}
		else {
			die "Invalid width $w";
		}
	}
	$spec->{finders} = \@finders;

	return $spec;
}


sub check_fasta
{
	my ($fn) = @_;
	if(defined($fn)) {
		open(IN, "<$fn"), or die "Cannot open file $fn: $!";
		my $is_fasta = 0;

		foreach my $ln (<IN>) {
			if($ln =~ /^(\S)/) {
				if($1 ne '>') {
					print "ERROR: '>' is not the first character of FASTA file\n";
					return 0;
				}
				last;
			}
		}
	}
	close(IN);
	return 1;
}

sub print_params
{
	my ($fptr, $spec) = @_;

	print $fptr ("Iteration slicer: " . $spec->{iters_slicer} . "\n");

}

sub print_result_param_table
{
	my ($spec) = @_;

	my $fn = sprintf("%s/param_table.txt", $spec->{statsig_dir});

	my $per_seq_str = ''; 
	if($spec->{per_seq_mode} eq 'zoops') {
		$per_seq_str .= sprintf("ZOOPS (prior-weight = %s)", $spec->{zoops_weight});
	}
	elsif($spec->{per_seq_mode} eq 'oops') {
		$per_seq_str .= 'OOPS';
	}
	my $genomic_file_str = '';
	if(defined($spec->{bgfile_fn})) {
		$genomic_file_str .= basename($spec->{bgfile_fn});
	}
	else {
		$genomic_file_str .= 'input FASTA file';
	}
	my @lst = (
		'Input file', basename($spec->{fasta_fn}),
		'Size of null set', sprintf("%d", $spec->{nullset_size}),
		'Occurrences per sequence', $per_seq_str,
		'Process run time limit', sprintf("-%s %d", $spec->{stop_crit}, $spec->{stop_crit_int}),
		'Convergence rate', sprintf("%d", $spec->{rapid_conv}),
		'Double-strand', ($spec->{ds} ? 'yes' : 'no'),
		'Order of Markov background', sprintf("%d", $spec->{markov_order}),
		'Genomic file', $genomic_file_str,
	);

	open(my $ostream, ">$fn") or die "Cannot open $fn for write: $!";
	if(scalar(@lst) % 2 != 0) {
		die "Error: list should be even";
	}
	for(my $i = 0; $i < scalar(@lst);  $i+=2) {
		printf $ostream ("\t<tr><td><b>%s</b></td> <td>%s</td></tr>\n", $lst[$i], $lst[$i+1]);
	}
	close($ostream);
}

sub write_config_readme
{
	my ($ARGV, $spec) = @_;
	my $fn = sprintf("%s/%s.config", $spec->{destdir}, $spec->{expname});
	open(my $ostream, ">$fn") or die "Cannot open $fn for write: $!";
	print $ostream join(' ', ($0, @$ARGV)) . "\n";
	close($ostream);
}

#"spec" struct
# expname - name of experiment
# qsub_jobname
# queue_names
# fasta_fn - just the filename without the directory path
# bgfile_fn 
# destdir - this should be an absolute path
# finders - list of finders
# nullset_size
# iters_slicer - number of (nullset) iters per job
# proc_per_machine - number of processes to run simultaneously on a single machine
# iter_mode_param
# tmpdir - a temporary directory within destdir
# nullset_dir
# is_batch_mode
# iterind - iteration index when under iteration mode
# statsig_dir
# $spec->{per_seq_mode} = 'zoops';
# $spec->{ds} = 1; #double-strand mode
# $spec->{zoops_weight} = 0.2;
# $spec->{stop_crit} = 'cput';
# $spec->{stop_crit_int} = 300;
# $spec->{markov_order} = 5;
# $spec->{rapid_conv} = 200;

#"finders" struct
# name - name of finder
# path
# param - parameters

#############################################################################################

sub main_batch_mode 
{
	my ($spec) = @_;

	printf "Preparing the job " . $spec->{expname} . "...\n";

	#create directory if not exist
	if(! -d $spec->{destdir}) {
		mkdir($spec->{destdir}) or die "Cannot make directory " . $spec->{destdir};
	}
	if(! -d $spec->{tmpdir}) {
		mkdir($spec->{tmpdir}) or die "Cannot make directory " . $spec->{tmpdir};
	}
	if(! -d $spec->{nullset_dir}) {
		mkdir($spec->{nullset_dir}) or die "Cannot make directory " . $spec->{nullset_dir};
	}
	if(! -d $spec->{statsig_dir}) {
		mkdir($spec->{statsig_dir}) or die "Cannot make directory " . $spec->{statsig_dir};
	}

	#print and write config
	write_config_readme(\@ARGV, $spec);
	print_params(*STDOUT, $spec);
	print_result_param_table($spec); #for gimsan_result.pl to use

	#get fileroot of fasta_fn and use the new path
	my $new_fasta_fn = sprintf("%s/%s", $spec->{destdir}, basename($spec->{fasta_fn}));
	copy($spec->{fasta_fn}, $new_fasta_fn) or die "Cannot copy file: $!"; #copy into directory
	$spec->{fasta_fn} = $new_fasta_fn;

	#create finders directories
	create_finders_dir($spec->{finders}, $spec->{destdir} );

	#create batch scripts
	my $batch_dir = sprintf("%s/batch_scripts/", $spec->{destdir});
	if(-d $batch_dir) {
		die "Directory $batch_dir already exists!\n";
	}
	else {
		mkdir($batch_dir) or die "Cannot make directory $batch_dir";
	}
	my $script_lst = generate_batchscripts($spec);
	my $batch_files = write_batchscripts($batch_dir, $script_lst);

	#generate nullset
	my $nullset_bgstr = '';
	if(defined($spec->{bgfile_fn})) {
		$nullset_bgstr = sprintf("--genome=%s", $spec->{bgfile_fn});
	}
	if($spec->{nullset_size} > 0) {
		system(sprintf("perl $NULLSET_GENERATOR --iters=%d --wndsize=50 --input=%s --output-dir=%s/nullset/ $nullset_bgstr ", 
				$spec->{nullset_size}, $spec->{fasta_fn}, $spec->{destdir}));
	}

	#submit job
	if($spec->{submit}) {
		print "Submitting jobs...\n";
		if($spec->{machines_alloc} == 1) {
			run_simple_batchscript($batch_files, $spec);
		}
		else {
			qsub_batchscripts($batch_files, $spec);
		}
	}
}

sub main_iter_mode 
{
	my ($spec) = @_;
	my $iterind = $spec->{iterind};

	if($iterind == -1) {
		die "No iteration index specified for non-batch mode\n";
	}

	print "Iterind: $iterind, Number of motif-finders: " . scalar( @{$spec->{finders}} ) . "\n";
	run_finders($spec, $iterind);
}

################################################################################
# Iteration mode (non-batch) 
################################################################################
sub resolve_fasta_fn_with_iterind 
{
	my ($spec, $iterind) = @_;

	my $fn = '';
	if($iterind > 0) {
		$fn = sprintf("%s/nullset/null-%05d.fa", $spec->{destdir}, $iterind);
	}
	else {
		#original input FASTA 
		#$fn = sprintf("%s/%s", $spec->{fasta_fn}, $spec->{destdir});; #fasta should be in absolute path
		$fn = $spec->{fasta_fn}; #fasta should be in absolute path
	}

	if(! -f $fn) {
		die "Error: the resolved file for motif-finder to run does not exist!\nResolved file: $fn\n";
	}
	
	return $fn;
}

# run the motif-finders
sub run_finders
{
	my ($spec, $iterind) = @_;
	foreach my $f(@{$spec->{finders}}) {
		my $stdout;
		my $stderr;
		if($iterind == 0) {
			$stdout = sprintf("%s/%s/motif-%05d.stdout", $spec->{destdir}, $f->{name}, $iterind);
			$stderr = sprintf("%s/%s/motif-%05d.stderr", $spec->{destdir}, $f->{name}, $iterind);
		}
		else {
			$stdout = sprintf("%s/%s/null-%05d.stdout", $spec->{destdir}, $f->{name}, $iterind);
			$stderr = sprintf("%s/%s/null-%05d.stderr", $spec->{destdir}, $f->{name}, $iterind);
		}
		my $fasta_fn = resolve_fasta_fn_with_iterind($spec, $iterind);
		my $path = $f->{path};
		my $param = $f->{param};
		my $cmd = "$path $fasta_fn $param 1>$stdout 2>$stderr";
		print "CMD: $cmd\n";
		system($cmd);
	}
}


################################################################################
# Generate/write/qsub batchscripts
################################################################################


sub generate_batchscripts
{
	my ($spec) = @_;

	my @script_lst = ();
	for(my $i = 0; $i <= $spec->{nullset_size}; $i+=$spec->{iters_slicer}) { #the 0-th iteration is the original FASTA input file
		my $last_iter = $i + $spec->{iters_slicer} - 1;
		if($last_iter > $spec->{nullset_size}) {
			$last_iter = $spec->{nullset_size};
		}
		my @iters = ( $i .. $last_iter );
		my $script =  generate_batchscripts_slicer($spec, \@iters);
		push(@script_lst, $script);
	}
	return \@script_lst;
}


sub create_finders_dir
{
	my ($finders, $destdir) = @_;

	foreach my $f(@$finders) {

		my $fdir = sprintf("$destdir/%s", $f->{name});
		if(! -d $fdir) {
			mkdir($fdir) or die "Cannot create directory $fdir"; ;
		}
	}
}

sub generate_batchscripts_pid_run
{
	my ($spec, $iter) = @_;
	my $cmd = '';

	my $tmpdir = $spec->{tmpdir};
	$cmd = sprintf("%s --iter=$iter %s  1>>$tmpdir/tmp.stdout 2>>$tmpdir/tmp.stderr &\n", abs_path($0), $spec->{iter_mode_param});
	$cmd .= "pid$iter=\$!\n";
}

sub generate_batchscripts_slicer
{
	my ($spec, $iters) = @_;
	my $script = '#!/bin/bash' . "\n";
	$script .= '#$ -S /bin/bash' . "\n";
	$script .= "\n";
	$script .= "export GIMSAN=$ENV{GIMSAN};\n";
	$script .= "\n";
	$script .= "cd $CWD\n";

	for(my $i = 0; $i < $spec->{proc_per_machine} && $i < scalar(@$iters); $i++) {
		$script .= generate_batchscripts_pid_run($spec, $iters->[$i]);
	}
	for(my $i = $spec->{proc_per_machine}; $i < scalar(@$iters); $i++) {
		$script .= sprintf("wait \$pid%d;\n", $iters->[ $i - $spec->{proc_per_machine} ]);
		$script .= generate_batchscripts_pid_run($spec, $iters->[$i]);
	}
	for(my $i = 0; $i < $spec->{proc_per_machine} && $i < scalar(@$iters); $i++) {
		$script .= sprintf("wait \$pid%d;\n", $iters->[scalar(@$iters)-1-$i]);
	}
	return $script;
}


sub write_batchscripts 
{
	my ($batch_dir, $script_lst) = @_;

	my @fn_lst = ();
	for(my $i = 0; $i < scalar(@$script_lst); $i++) {
		my $outfn = sprintf("$batch_dir/batch_script_%05d.sge", $i+1);
		open(my $fh, ">$outfn") or die "Cannot open $outfn for write.\n";
		print $fh ($script_lst->[$i]);
		close($fh);
		push (@fn_lst, $outfn);
	}
	return \@fn_lst;
}

sub run_simple_batchscript
{
	my ($batchfiles, $spec) = @_;
	
	my $tmpdir = $spec->{tmpdir};
	my $qsub_jobname = $spec->{qsub_jobname};
	if(scalar(@$batchfiles) != 1) {
		die "Error: invalid number of batch-files in run_simple_batchscript()";
	}
	my $file = $batchfiles->[0];

	my $cmd = "bash $file 1>$tmpdir/tmp.stdout 2>$tmpdir/tmp.stderr";
	print "$cmd\n";
	system($cmd);
}

sub qsub_batchscripts
{
	my ($batchfiles, $spec) = @_;
	my $queue_str = '';
	if(defined($spec->{queue_names})) {
		if($spec->{queue_names} =~ /[A-Za-z]/) {
			$queue_str = '-q ' . $spec->{queue_names};
		}
		else {
			die "Invalid queue names " . $spec->{queue_names};
		}
	}	
	my $qsub_jobname = $spec->{qsub_jobname};
	my $tmpdir = $spec->{tmpdir};

	foreach my $f (@$batchfiles) {
		#system("qsub -N $qsub_jobname -q $queue_names -pe make 1 -e ~/sge_tmp -o ~/sge_tmp $f");
		system("qsub -N $qsub_jobname $queue_str -e $tmpdir -o $tmpdir $f");
		sleep(0.1);
	}
}


###################################################################################

sub gimsan_env_help
{
	print "For example, you can initialize the environment variable in bash by:\n";
	print "\n";
	print "    export GIMSAN=~/gimsan-cmdline/\n";
	print "\n";
}


sub resolve_tilde 
{
	my $file = shift;
	$file =~ s{ ^ ~ ( [^/]* ) }
    	          { $1
        	            ? (getpwnam($1))[7]
            	        : ( $ENV{HOME} || $ENV{LOGDIR}
                	         || (getpwuid($>))[7]
                    	   )
	}ex;
	return $file;
}
