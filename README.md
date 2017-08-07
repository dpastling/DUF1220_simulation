
# WGS Simulation of DUF1220 Regions

Simulation of copy number changes within the DUF1220 domains and testing several alignment and summarization strategies. The details will be published in the
forthcoming paper:

> Astling, DP, Heft IE, Jones, KL, Sikela, JM. "High resolution measurement of
> DUF1220 domain copy number from whole genome sequence data"

## Overview

The types of simulations provided by this code

1. A spike-in study where reads from an individual DUF1220 domain are simulated and aligned back to the genome. The purpose of this is to assess the read alignment ambiguity for each domain.

2. A baseline study to simulate the coverage one might expect for a WGS experiment. Rather than simulating reads for the whole genome, reads are simulated from all DUF1220 domains. The purpose of this study is to assess the accuracy of each alignment strategy to account for all 271 DUF1220 copies.

3. Test the sequencing parameters such as read length, single versus paired-end reads, coverage, quality scores, and sequencing errors on the accuracy of each method. Each of these parameters can be adjusted.

Note: these simulations could be adapted for other segmental duplications and genomic regions.


## Dependencies

See (plethora)[https://github.com/dpastling/plethora] for installing bedtools, samtools, and bowtie2.

- bowtie2 version 2.2.9
- bedtools version 2.17.0
- samtools version: 0.1.19-44428cd
- Perl module: Math::Random
- Perl module: Math::Complex

If you are interesting in comparing the different alignment strategies described in the paper, you will also need to install mrsFastUltra and bowtie1.

- (mrsFast-Ultra)[http://sfu-compbio.github.io/mrsfast/] version 3.3.11
- (Bowtie1)[http://bowtie-bio.sourceforge.net/index.shtml] version 1.1.2

After installing mrsFastUltra and Bowtie1, you will need to build separate genome indicies for each. Follow the instructions included for each program for building these genome references.


## Getting Started




