#!/usr/bin/perl

# This script accepts the output of bedtools intersect using the -wao option
# and applies multi-read correction. It does this by dividing the base overlap
# by the read frequency. For example if a 100 bp read aligns to 5 positions
# its contribution should be 20bp at each position. Whereas a read that aligns
# uniquely will contribute a full 100 bp.
#
# Follow up by running "bedtools merge -scores sum -i <output of this file>"

use strict;
use warnings;

my %read_frequency;

if (!@ARGV)
{
	print "Usage:\t$0 <intersected bed file with -wao>\n";
	exit;
}

# Read through file once to count the frequency of each read
open(INFILE, $ARGV[0]) or die;
while(<INFILE>)
{
	chomp;
	my $line = $_;
	my @data = split(/\t/, $line);
	next if ($read_frequency{$data[9]} eq ".");
	$read_frequency{$data[9]}++;	
}
close(INFILE);

# read through file again and adjust the contribution of each read
open(INFILE, $ARGV[0]) or die;
while(<INFILE>)
{
	chomp;
	my $line = $_;
	my @data = split("\t", $line);
	my ($chr_ref, $pos1_ref, $pos2_ref, $domain, $score, $strand_ref, $chr_read, $pos1_read, $pos2_read, $read_id, $quality, $strand_read, $base_overlap) = split("\t", $line);
	$base_overlap = $base_overlap + 1;
	my $adjusted_depth;
	if ($read_frequency{$read_id})
	{
		$adjusted_depth = $base_overlap / $read_frequency{$read_id};
	} else {
		$adjusted_depth = $base_overlap;
	}
	my $result = join("\t", $domain, $pos1_ref, $pos2_ref, $chr_ref, $adjusted_depth);
	print "$result\n";
}
close(INFILE);


