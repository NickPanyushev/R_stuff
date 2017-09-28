#!/usr/bin/perl -w
use strict;
use warnings;
use subs 'abs_path';
use Cwd;
use File::Copy;
use File::Basename;

my $DEBUG = 0;
if(@ARGV < 1) {
	print_help();
}

if(!defined($ENV{R_BIN})) {
	print STDERR "Error: R_BIN environment variable is undefined.\n";
	print STDERR "\n";
	print STDERR "Please install the statistical computing software R and then set the R_BIN environment variable.\n";
	print STDERR "For example, you can initialize the environment variable in bash by:\n";
	print STDERR "\n";
	print STDERR "    export R_BIN=/usr/local/R-2.6.1/bin/\n";
	print STDERR "\n";
	exit;
}

if(!defined($ENV{GIMSAN})) {
	print "Error: GIMSAN environment variable is undefined.\n";
	print "\n";
	gimsan_env_help();
	exit;
}

my $STAT_EVAL = abs_path("$ENV{GIMSAN}/bin/extract_scores_standalone.pl");
my $WEBLOGO_EXE = abs_path("$ENV{GIMSAN}/weblogo/seqlogo");
my $R_EXE = abs_path("$ENV{R_BIN}/R");
my $COL_DEPEND_EXE = abs_path("$ENV{GIMSAN}/column_dependency_app/column_dependency.out");
my $CWD = getcwd();
main();



sub main 
{

	my $expdir = abs_path($ARGV[0]);

	if(! -f $STAT_EVAL) {
		print "Statistical evaluation engine cannot be found at:\n $STAT_EVAL\n\n";
		print "Please check GIMSAN environment variable.\n";
		print "\n";
		gimsan_env_help();
		exit(1);
	}
	if(! -f $WEBLOGO_EXE) {
		print "WebLogo generator cannot be found at:\n $WEBLOGO_EXE\n\n";
		print "Please check GIMSAN environment variable.\n";
		print "\n";
		gimsan_env_help();
		exit(1);
	}
	if(! -f $R_EXE) {
		print "R cannot be found at:\n $R_EXE\n\n";
		print "Please check R_BIN environment variable.\n";
	}
	if(! -f $COL_DEPEND_EXE) {
		print "Column-dependency software cannot be found at:\n $COL_DEPEND_EXE\n";
		print "It may not have compiled correctly.\n";
		exit(1);
	}

	#check for batch_scripts, etc
	if(! -d "$expdir/batch_scripts") {
		print "ERROR: invalid experiment directory - does not contain batch_scripts directory\n";
		print "$expdir/batch_scripts does not exist!\n";
		exit(1);
	}

	my $res = construct_result_struct($expdir);

	#print html
	open(my $ostream, '>'. $res->{html_fn}) or die "Cannot open " . $res->{html_fn} . " for write: $!";
	print_html_header($ostream, $res->{expname});
	print_header_config($ostream, $res);
	foreach my $f (sort {$a->{width} <=> $b->{width}} @{$res->{finders}}) {
		print_finder_result($ostream, $res, $f);
	}
	print_html_footer($ostream);
	close($ostream);
}


sub print_header_config
{
	my ($ostream, $res) = @_;

	open(IN, "<" . $res->{param_table_fn}) or die "Cannot open " . $res->{param_table_fn} . " for read: $!";
	printf $ostream ("GIMSAN output from job <b>%s</b> <br /><br />\n", $res->{expname});
	printf $ostream ("<b>Run parameters (<a href=\"%s\">command-line</a>):</b><br>\n", $res->{expname} . '.config');
	print $ostream "<table border=\"2\" cellspacing=\"1\" cellpadding=\"10\">\n";
	foreach my $ln (<IN>) {
		$ln =~ s/[\r\n]//g;
		print $ostream "$ln\n";
	}
	print $ostream "</table>\n";
	print $ostream "<br><hr>\n";
	close(IN);
}

sub print_html_header
{
	my ($ostream, $expname) = @_;
	my $txt = <<"End";
<html>
<body>
<title>$expname</title>
<style type="text/css">
body { margin: 20px;
	  color: black;
	  background: white;
	  font-family:arial,sans-serif;
	  font-size: 10pt;}
ul li { padding: 5;}
a {color: #0033FF; text-decoration: underline; font-weight: none;}
a:visited {color: #8B0000; text-decoration: underline; font-weight: none;}
a:hover {color: #0033FF; text-decoration: underline; font-weight: none;}
a:visited:hover {color: #8B0000; text-decoration: underline; font-weight: none;}
table { font-size: 10pt; }
</style>
<center>
End
	print $ostream $txt;
}

sub print_html_footer
{
	my $ostream = shift;
	print $ostream "</center>\n";
	print $ostream "</body>\n";
	print $ostream "</html>\n";
}

sub print_help
{
	print "DESCRIPTION\n";
	print "  Generates a html output of a GIMSAN job and performs motif significance evaluation\n";
	print "\n";
	print "USAGE\n";
	print "  perl gimsan_result.pl <experiment directory>\n";
	print "\n";
	exit(1);
}

sub construct_result_struct
{
	my ($expdir) = @_;

	my $expname = basename($expdir);
	$expname =~ s/[\/\:\\]//g; #remove : and / and \ characters

	#construct finders
	my $finders = search_finders($expdir);

	my %hash = (
		expname => $expname,
		expdir => $expdir,
		statsig_dir => "$expdir/statsig/",
		param_table_fn => "$expdir/statsig/param_table.txt",
		html_fn => "$expdir/output.html",
		finders => $finders,
	);
	return \%hash;
}

sub search_finders
{
	my ($expdir) = @_;

	my @finders = ();
	opendir(DIR, $expdir) || die "Cannot opendir $expdir: $!";
	while(my $f = readdir(DIR)) {
		if(-d "$expdir/$f" && $f =~ /^width(\d+)/) { #directory that begins with "width"
			(my $name = $f) =~ s/\W//g; #directory should be alphanumeric
			my $motif_res_fn = "$expdir/$name/motif-00000.stdout";
			my $motif_res_info = parse_finder_output($motif_res_fn);
			my $count = 0;
			while(-f sprintf("$expdir/$name/null-%05d.stdout", $count+1)) {
				$count++;
			}
			my %hash = (
				name => $name,
				width => int($1),
				out_dir => "$expdir/$name/",
				motif_res_fn => $motif_res_fn,
				nullset_size => $count, 
				numsites => $motif_res_info->{numsites},
			);
			push(@finders, \%hash);
		}
	}
	closedir DIR;

	if(scalar(@finders) == 0) {
		print "ERROR: there are no results found in $expdir.\n";
		exit(1);
	}

	my $nullset_size = $finders[0]->{nullset_size};
	foreach my $f(@finders) {
		print STDERR $f->{name} . " has null set size of " . $f->{nullset_size} . "\n";
	}
	foreach my $f(@finders) {
		if($nullset_size != $f->{nullset_size}) {
			print STDERR "WARNING: inconsistent null-set size. GIMSAN job did not complete.\n";
		}
	}
	
	return \@finders;
}

sub print_finder_result
{
	my ($fh, $res, $finder) = @_;

	my $rscript_stdout = '';
	if($finder->{nullset_size} > 0) {
		#Create R-script
		my $nullscores_out = sprintf("%s/scores.%s", $res->{statsig_dir}, $finder->{name});
		my $r_out = sprintf("%s/Rscript_%s.R", $res->{statsig_dir}, $finder->{name});
		my $cmd = sprintf("perl $STAT_EVAL --ref-result=%s --nullset-dir=%s --nullscores-out=$nullscores_out --R-out=$r_out",
			$finder->{motif_res_fn}, $finder->{out_dir});
		print STDERR "$cmd\n";
		system($cmd);

		#Run R-script
		$rscript_stdout = `$R_EXE -f $r_out 2>&1`;
	}
	
	#Print span and number of sequences
	print $fh "<br />\n";
	printf $fh ("<b>span: %d</b>, logo constructed from %d sequences<br />\n", $finder->{width}, $finder->{numsites});
	
	#Print R-script output
	foreach my $ln (split(/\n/, $rscript_stdout)) {
		$ln =~ s/[\r\n]//g;
		if($ln =~ /MLE/) {
			print $fh "$ln<br />\n";
		}
	}
	print $fh "<br />\n";

	#k-mers file
	my $kmers_fn = sprintf("%s/%s.kmers", $res->{statsig_dir}, $finder->{name});
	write_kmers_file($finder, $kmers_fn);

	#Generate and print WebLogo
	generate_weblogo($finder, $kmers_fn, sprintf("%s/pict.%s.png", $res->{expdir}, $finder->{width}));
	printf $fh ("<IMG height=\"200\" alt=\"WebLogo\" src=\"pict.%d.png\" /> <br />\n", $finder->{width});

	#Generate and print column-dependency
	my $num_coldep_pairs = run_coldepend($finder, $kmers_fn, sprintf("%s/coldep.%d", $res->{expdir}, $finder->{width}));
	printf $fh ("<a href=\"coldep.%d\">Column pairs with statistically significant dependency (%d pairs) </a><br />\n", $finder->{width}, $num_coldep_pairs);
	
	#print link for details
	printf $fh ("<a href=\"%s/%s\">Motif finder detailed output</a><br /><br />\n", $finder->{name}, 'motif-00000.stdout');
	print $fh "<hr>\n";
}

sub run_coldepend
{
	my ($finder, $kmers_fn, $dest_fn) = @_;
	my $stderr_fn = "$dest_fn.stderr";
	my $cmd = ("$COL_DEPEND_EXE -fsa $kmers_fn 1>$dest_fn 2>$stderr_fn");
	print STDERR "$cmd\n";
	system("$cmd");

	#extract number of pairs
	open(IN, "<$dest_fn") or die "File $dest_fn not found";
	my $num_significant_pairs = undef;
	foreach my $ln (<IN>) {
		if($ln =~ /for statistically significant pairs \((\d+) pairs\)/) {
			$num_significant_pairs = $1;
		}
	}
	close(IN);
	if(!defined($num_significant_pairs)) {
		die "Error: column-dependency did not complete";
	}
	return $num_significant_pairs;
}

sub generate_weblogo
{
	my ($finder, $kmers_fn, $dest_fn) = @_;
	my $logo_opt = "-F png -caYnMb";
	my $logo_width = (20/45) * $finder->{width} + 2;
	my $cmd = "$WEBLOGO_EXE -f $kmers_fn $logo_opt -w $logo_width > $dest_fn";
	print STDERR "$cmd\n";
	system($cmd);
}

sub write_kmers_file
{
	my ($finder, $dest_fn) = @_;

	my $in_fn = $finder->{motif_res_fn};
	my $numsites = undef;
	my @kmers = ();

	open(IN, "<$in_fn") or die "Cannot open $in_fn for read: $!";
	foreach my $ln (<IN>) {
		$ln =~ s/[\r\n]//g;
		if($ln =~ /^\s*\d+\s+[acgtnx]*\s+[\-\d]+ ([ACGTNX]+)\s+[\-\d]+\s+[acgtnx]*\s+\(/) {
			push (@kmers, $1);
		}
	}
	close(IN);

	if(scalar(@kmers) != $finder->{numsites}) {
		die "Internal ERROR: the number of sites does not match number of k-mers";
	}

	#write to file
	open(my $ostream, "> $dest_fn") or die "Cannot open $dest_fn for write: $!";
	foreach my $k (@kmers) {
		if(length($k) != $finder->{width}) {
			print STDERR "Warning: $k does not the correct width of " . $finder->{width} . "\n";
		}
		print $ostream ">FASTA header\n";
		print $ostream "$k\n";
	}
	close($ostream);
}
sub parse_finder_output
{
	my ($fn) = @_;
	open(IN, "<$fn") or die "Cannot open $fn for read:$!\n";
	my %hash = (
		numsites => undef,
	);
	foreach my $ln (<IN>) {
		$ln =~ s/[\r\n]//g;
		if($ln =~ /Number of predicted sites:\s*(\d+)/) {
			$hash{numsites} = int($1);
		}
	}
	return \%hash;
}
	
sub gimsan_env_help 
{
	print "For example, you can initialize the environment variable in bash by:\n";
	print "\n";
	print "    export GIMSAN=~/gimsan-cmdline/\n";
	print "\n";
}

sub resolve_tilde
{
        #taken from Perl Cookbook Recipe 7.3 "Expanding Tildes in Filenames"
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
