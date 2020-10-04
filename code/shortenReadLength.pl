#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
Getopt::Long::Configure("bundling");

my $read_length;
my $out_file;

GetOptions( "l=i" => \$read_length,
            "o=s" => \$out_file
          );

if (!@ARGV || @ARGV > 1 || (!$read_length))
{
        print "$0:  shorten read length of a FASTQ file\n";
        #print "------------------------------------------------------------------------------------------\n";
        print "\n";
        print "Usage:\t$0 [options] <fastq file> \n";
        print "\n";
        print "Options:\n";
        print "    -l <int>   read length                           \n";
        print "    -o <file>  file to write results to (default: print to stdout)\n";
        print "\n";
        exit;
}

open(INFILE, $ARGV[0]) or die "can't open file: $ARGV[0]\n";

if ($out_file)
{
	open(OUTFILE, ">".$out_file) or die "can't write to file: $out_file\n";
}


while (<INFILE>)
{
	my $id       = $_;	
	my $sequence = <INFILE>;
	my $comment  = <INFILE>;
	my $quality  = <INFILE>;
	
	if ($id !~ /^\@/ | $comment !~ /^\+/ | $sequence !~ /^[ACTGUNactgun]+$/)
	{
		warn "Error: malformed FASTQ file\n";
		exit;
	}
	
	chomp $sequence;
	chomp $quality;
	
	$sequence = shortenLine($sequence, $read_length);
	$quality  = shortenLine($quality, $read_length);
	
	if ($out_file)
	{
		print OUTFILE $id;
		print OUTFILE $sequence . "\n";
		print OUTFILE $comment;
		print OUTFILE $quality . "\n";
	} 
	else
	{
		print $id;
		print $sequence . "\n";
		print $comment;
		print $quality . "\n";
	}
}


close(INFILE);

if ($out_file)
{
	close(OUTFILE);
}


sub shortenLine
{
	my ($line, $read_length) = @_;
	
	chomp $line;
	if (length($line) <= $read_length)
	{
		warn "Error: your desired read length: $read_length is greater than or equal to the actual read length: " . length($line)  . "\n";
		exit;
	}
	$line = substr($line, 0, $read_length);
	
	return($line);
}

