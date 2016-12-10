#!/usr/bin/perl

use strict;
use warnings;
use Math::Random;
use Getopt::Long;
Getopt::Long::Configure("bundling");


my $read_length     = 100;
my $number_of_reads = 1000000;
my $random_seed     = time;
my $out_file;
my $out_file_2;
my $paired = 0;
my $fragment_length = 350;
my $fragment_length_sd = 50;
my $coverage;

GetOptions (    "n=i"        => \$number_of_reads,
                "c=f"        => \$coverage,
                "o=s"        => \$out_file,
		"l=i"        => \$read_length,
		"insert=i"   => \$fragment_length,
		"insertsd=i" => \$fragment_length_sd,
                "s=s"        => \$random_seed,
		"paired"     => \$paired,
           );  

if (!@ARGV || @ARGV > 1 || ($paired & ! $out_file))
{
	print "$0:  generate random reads from regions specified in .bed file. Output reads in bed format\n";
	print "------------------------------------------------------------------------------------------\n";
	print "Usage:\t$0 [options] <bed file> \n";
	print "\n";
	print "Options:\n";
	print "    -n <int>   number of reads                       default: " . $number_of_reads . "\n";
	print "    -c <num>   coverage or desired read depth        default: none\n";
	print "    -l <float> read length                           default: " . $read_length. "bp\n";
	print "    --insert   insert size (length of fragment)      default: " . $fragment_length . "bp\n";
	print "    --insertsd insert size standard deviation        default: " . $fragment_length_sd . "bp\n";
	print "    -s <int>   seed for the random number generator  default: time\n";
	print "    -o <file>  file prefix to write results to       default: print to stdout\n";
	print "    --paired   generate paired end reads             default: single end reads\n";
	print "               creates two files with _1 and _2 \n";
	print "               extension for each pair\n";
	print "               Note: you must provide a filename for paired output\n";
	print "\n";
	exit;
}

my $bed_reference = $ARGV[0];
my @chromosome_name;
my @chromosome_size;
my @chromosome_offset;
my @genome_position;
my $total_size = 0;

open(BED_REF, $bed_reference) or die "cannot open $bed_reference";
while(<BED_REF>)
{
	my $line = $_;
	chomp $line;
	my @elements = split(/\t/, $line);
	my $chr_length = $elements[2] - $elements[1] + 1;
	push(@chromosome_name, $elements[0]);
	push(@chromosome_size, $chr_length);
	push(@chromosome_offset, $elements[1]);
	push(@genome_position, 1 + $total_size);
	$total_size += $chr_length;
}
close(BED_REF);

if ($coverage)
{
	$number_of_reads = ($total_size * $coverage) / $fragment_length;
#	$number_of_reads = $number_of_reads / 2 if ($paired);
	$number_of_reads = sprintf("%.0f", $number_of_reads);
}

if ($paired)
{
	# this check is necessary when sampling from short regions
	my $are_any_chromosomes_longer_than_the_fragment = 0;
	for (my $i = 0; $i <= $#chromosome_size; $i++)
	{
		if ($chromosome_size[$i] > $fragment_length)
		{
			$are_any_chromosomes_longer_than_the_fragment = 1;
		} else {
			print STDERR "$chromosome_name[$i] is shorter than the fragment length of $fragment_length\n"
		}
	}
	if (! $are_any_chromosomes_longer_than_the_fragment)
	{
		print STDERR "None of the .bed regions are longer than $fragment_length\n";
		print STDERR "Either shorten the fragment length or sample from a different sequence file\n";
		exit;
	}
}

if (system("bedtools --version >&/dev/null"))
{
	print STDERR "bedtools is not present in your \$PATH. Please add or install\n";
	exit;
}

 
if ($out_file)
{
	if ($out_file !~ /\.bed$/)
	{
	 	$out_file .= ".bed";
	}
	if ($paired)
	{
		$out_file =~ s/\.bed/_1.bed/g;
		$out_file_2 = $out_file;
		$out_file_2 =~ s/_1\.bed/_2.bed/;
		open(OUTFILE2, ">", $out_file_2) or die "can't write to file: $out_file_2\n";
	}
	open(OUTFILE, ">", $out_file) or die "can't write to file: $out_file\n";
} else 
{
	*OUTFILE = *STDOUT;
}

random_set_seed_from_phrase($random_seed);
# srand($random_seed);

for (my $i = 1; $i <= $number_of_reads; $i++)
{
	# my $random_position = sprintf("%.0f", rand($total_size));
	my $random_position = sprintf("%.0f", random_uniform(1, 1, $total_size));
	my $fragment_length = sprintf("%.0f", random_normal(1, $fragment_length, $fragment_length_sd));
	for (my $j = 0; $j <= $#genome_position; $j++)
	{
		next if ($random_position < $genome_position[$j] || $random_position > $genome_position[$j] + $chromosome_size[$j]);
		my $read_mid_point  = $random_position - $genome_position[$j] + $chromosome_offset[$j];
		my $first_position  = $read_mid_point - sprintf("%.0f", $read_length / 2);
		my $second_position = $first_position + $read_length - 1;
		if ($second_position > $chromosome_offset[$j] + $chromosome_size[$j] || 
			$first_position < $chromosome_offset[$j])
		{
			$i--;
			last;
		}
		my $strand = rand();
		if (sprintf("%.0f", $strand) >= 0.5)
		{
			$strand = "+";
		} else {
			$strand = "-";
		}
		

		if ($paired && $out_file_2)
		{
			my $read2_strand;
			my $read2_first_position;
			my $read2_second_position;
			my $mate_distance = $fragment_length - $read_length;
			# Note that the second read is on the opposite strand and points in the other direction
			if ($strand eq "-")
			{
				$read2_strand = "+";
				$read2_first_position  = $first_position - $mate_distance;
				$read2_second_position = $read2_first_position + $read_length;
			} else {
				$read2_strand = "-";
				$read2_first_position  = $first_position + $mate_distance;
				$read2_second_position = $read2_first_position + $read_length;
			}
			
			if ($read_length >= abs($fragment_length))
			{
				$read2_first_position  = $first_position;
				$read2_second_position = $second_position;
			}
			if ($read2_first_position  < $chromosome_offset[$j] || 
			    $read2_second_position < $chromosome_offset[$j] || 
			    $read2_first_position  > $chromosome_offset[$j] + $chromosome_size[$j] || 
			    $read2_second_position > $chromosome_offset[$j] + $chromosome_size[$j]
			   )
			{
				$i--;
				last;
			}
			# Note: The read_ids for both pairs must match
			# optionally they can end with #0/1 for first read and #0/2 for second read
			print OUTFILE2 "$chromosome_name[$j]\t$read2_first_position\t$read2_second_position\t$i\:$chromosome_name[$j]\:$first_position\-$second_position:$strand#0/2\t255\t$read2_strand\n";
		}
		print OUTFILE "$chromosome_name[$j]\t$first_position\t$second_position\t$i:$chromosome_name[$j]\:$first_position\-$second_position:$strand#0/1\t255\t$strand\n";
	}
}

if ($out_file) 
{ 
	close(OUTFILE);
	#cleanup_file($out_file);
}

if ($out_file_2) 
{ 
	close(OUTFILE2); 
	#cleanup_file($out_file_2);
}

sub cleanup_file
{
	my $outfile = shift;
	my $sorted_bed = $outfile;
        my $fasta_result = $outfile;
        $sorted_bed =~ s/\.bed/_sorted.bed/;
        $fasta_result =~ s/\.bed/.fa/;
        system("bedtools getfasta -s -name -fi $bed_reference -bed $outfile -fo $fasta_result");
        system("rm $outfile");
}


sub reverse_complement 
{
	my $dna = shift;
	
	# reverse the DNA sequence
	my $revcomp = reverse($dna);
	
	# complement the reversed DNA sequence using full IUPAC bases
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	return $revcomp;
}

