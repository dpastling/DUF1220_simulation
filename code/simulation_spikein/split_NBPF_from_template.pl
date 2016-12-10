#!/usr/bin/env perl

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;
Getopt::Long::Configure("bundling");

my $numberOfJobs = 10;
my $bed_ref_file;
my $template_prefix;
my $output_folder;

GetOptions (
	        "i=s"        => \$bed_ref_file,
                "t=s"        => \$template_prefix,
                "o=s"        => \$output_folder,
		"n=i"	     => \$numberOfJobs
           );  

if (! $bed_ref_file || ! $template_prefix || ! $output_folder)
{
	print "$0:  split simulated reads by NBPF genes\n";
	print "--------------------------------------------------\n";
	print "Usage:\t$0 -i <bed_reference.bed> -t <simulated_reads.bed> -o <output folder> -n number_of_jobs\n";
	print "\n";
	exit;
}

my @bed_ref;
open(BED_REF, $bed_ref_file) or die "cannot open the bed reference: $bed_ref_file\n";
while(<BED_REF>)
{
	chomp;
	push @bed_ref, [ split(/\t/, $_) ];
}
close(BED_REF);

my $sample_prefix = $template_prefix;
$sample_prefix =~ s/^.+?\/template_([^\/]+?)/$1/; 

# since we plan to append to each file, we need to remove any existing files beforehand 
for (my $i = 0; $i <= $#bed_ref; $i++)
{
	my $domain = $bed_ref[$i][3];
	my $gene = $domain;
	$gene =~ s/^(NBPF[^_]+)_.+?$/$1/;
	my @old_files = glob("$output_folder/$gene/$domain\_$sample_prefix*");
	foreach my $file (@old_files)
	{
		system("rm $file");
	}
	if (! -e "$output_folder/$gene") 
	{
		mkdir "$output_folder/$gene";
	}
}


open(TEMPLATE1, "$template_prefix\_1.bed") or die "cannot open TEMPLATE1: $template_prefix\_1.bed\n";
open(TEMPLATE2, "$template_prefix\_2.bed") or die "cannot open TEMPLATE2: $template_prefix\_2.bed\n";
my $pm = new Parallel::ForkManager($numberOfJobs);
while(<TEMPLATE1>)
{
	chomp;
	my @coords_1 = split(/\t/, $_);

	my $line2 = <TEMPLATE2>;
	chomp $line2;
	my @coords_2 = split(/\t/, $line2);

	my $pid = $pm->start;
	if ($pid != 0)
	{
		## for the parent process
		next;
	}

	for (my $i = 0; $i <= $#bed_ref; $i++)
	{
		my $hit = 0;
		next if ($bed_ref[$i][0] ne $coords_1[0]);
		$hit = 1 if ($bed_ref[$i][1] >= $coords_1[1] && $bed_ref[$i][1] <= $coords_1[2]); 
		$hit = 1 if ($bed_ref[$i][2] >= $coords_1[1] && $bed_ref[$i][1] <= $coords_1[2]);
		$hit = 1 if ($bed_ref[$i][1] >= $coords_2[1] && $bed_ref[$i][1] <= $coords_2[2]); 
		$hit = 1 if ($bed_ref[$i][2] >= $coords_2[1] && $bed_ref[$i][1] <= $coords_2[2]);
		next if (! $hit);
		my $domain = $bed_ref[$i][3];
		my $gene = $domain;
		$gene =~ s/^(NBPF[^_]+)_.+?$/$1/;
		my $suffix = $template_prefix;
		my $file1 = "$output_folder/$gene/$domain\_$sample_prefix\_1.bed";
		my $file2 = "$output_folder/$gene/$domain\_$sample_prefix\_2.bed";
		open(OUT1, ">>", $file1) or die("cannot write results\n");
		print OUT1 join("\t", @coords_1) . "\n";
		close(OUT1);
		open(OUT2, ">>", $file2) or die("cannot write results\n");
		print OUT2 join("\t", @coords_2) . "\n";
		close(OUT2);
		last;
	}
	$pm->finish;
}
$pm->wait_all_children;
close(TEMPLATE1);
close(TEMPLATE2);

