#!/usr/bin/env perl

use strict;
use warnings;
use Math::Random;
use POSIX;

use Getopt::Long;
Getopt::Long::Configure("bundling");

my $fasta_file;
my $error_file;
my $outfile;

GetOptions (
	"f=s"	=> \$fasta_file,
	"e=s"	=> \$error_file,
	"o=s"	=> \$outfile
);

my @mean_qualities;
my @sd_qualities;
my @error;
open(ERROR, $error_file) or die("cannot open the error file\n");
while(<ERROR>)
{
	chomp;
	my $line = $_;
	next if ($line =~ /^base/);
	my @elements = split(/\t/, $line);
	push(@mean_qualities, $elements[1]);
	push(@sd_qualities, $elements[2]);
	push(@error, $elements[3]);	
}
close(ERROR);

my $out_handle;
if ($outfile)
{
	open($out_handle, ">" . $outfile) or die "cannot write to the outfile: $outfile\n";
}
if(! $outfile)
{
	$out_handle = *STDOUT;
}

## If you want to control the seed, here's how you do that:
# my $random_seed = time;
#random_set_seed_from_phrase($random_seed);



open(FASTA, $fasta_file) or die("cannot read the FASTA file\n");
while(<FASTA>)
{
	my $id = $_;
	my $seq = <FASTA>;
	
	chomp $id;
	chomp $seq;
	
	my @seq = split(//, $seq);
	my $qual;
	my @qual;
	
	if (length(@seq) > length(@mean_qualities) + 1)
	{
		print "Error: sequence length is longer than the error profile\n";
		exit 1;
	}
	
	if (length(@seq) == length(@mean_qualities) + 1)
	{
		$mean_qualities[$#mean_qualities + 1] = $mean_qualities[$#mean_qualities];
		$sd_qualities[$#sd_qualities + 1] = $sd_qualities[$#sd_qualities];
		$error[$#error + 1] = $error[$#error];
	}
		
	my @mutations;
	
	for (my $pos = 0; $pos <= $#seq; $pos++)
	{
		my $mean_qual = $mean_qualities[$pos];
		my $sd_qual = $sd_qualities[$pos];
		my $error_rate = $error[$pos];
		
		my $rand_qual = random_normal(1, $mean_qual, $sd_qual);
		if ($rand_qual > 40) {$rand_qual = 40};
		if ($rand_qual < 0) {$rand_qual = 0};
		
		$qual[$pos] = chr($rand_qual + 33);
		
		my $rand_error = random_uniform();
		if ($rand_error <= $error_rate)
		{
			push(@mutations, $pos + 1);
			my $current_base = $seq[$pos];
			my $other_bases = "ATCG";
			$other_bases =~ s/$current_base//;
			my $rand_base = random_uniform_integer(1, 0, 2);
			$seq[$pos] = substr($other_bases, $rand_base, 1)
		}
	}

	$id   =~ s/^>/@/;	
	$seq  = join("", @seq);
	$qual = join("", @qual);

	print $out_handle "$id\n";
	print $out_handle "$seq\n";
	if (@mutations)
	{
		print $out_handle "+ mutated bases: " . join(" ", @mutations) . "\n";
	} else {
		print $out_handle "+\n";
	}
	print $out_handle "$qual\n"	
}
close(FASTA);
close($out_handle);


