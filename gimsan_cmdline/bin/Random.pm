#!/usr/bin/perl -w
use strict;

package Random;

my $RAND_GENERATOR = "$ENV{GIMSAN}/gibbsmarkov/rand.out";
my $SEED_MAX = 4294967295; #unsigned int

sub new 
{
	my $capacity;
	if(scalar(@_) == 2) {
		my $class = shift;
		$capacity = shift;
	}
	else {
		die "incorrect usage";
	}
	
	my $set = gen_randset($capacity);
	my $self = {
		capacity => $capacity,
		randset => $set,
	};
	bless $self;
	return $self;
}

sub random_uniform
{
	my $self = shift;
	my $set = $self->{randset};
	if(scalar(@$set) == 0) {
		push(@$set, @{gen_randset($self->{capacity})});
	}
	return pop(@$set);
}

sub gen_randset
{
	my ($size) = @_;
	my $seed = rand() * $SEED_MAX;
	
	if(!-f $RAND_GENERATOR) {
		print STDERR "Fatal ERROR: cannot find random number generator at:\n$RAND_GENERATOR\n";
		print STDERR "\n";
		print STDERR "Check whether GIMSAN environment variable is set correctly and that\n";
		print STDERR "the random-number generator has been compiled correctly.\n";
		print STDERR "\n";
	}

	my $out = `$RAND_GENERATOR -n $size -s $seed`;
	my @set = split(/\s+/, $out);
	if(scalar(@set) != $size) {
		print STDERR "Error in gen_randset(): found " . scalar(@set) . " when size=$size\n";
		print STDERR "$out\n";
		die;
	}
	return \@set;

}

1;
__END__

