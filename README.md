
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

See [plethora](https://github.com/dpastling/plethora) for installing bedtools, samtools, and bowtie2.

- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) version 2.2.9
- [Bedtools](http://bedtools.readthedocs.io/en/latest/) version 2.17.0
- [Samtools](http://samtools.sourceforge.net) version: 0.1.19-44428cd
- Perl module: Math::Random
- Perl module: Math::Complex

If you are interesting in comparing the different alignment strategies described in the paper, you will also need to install mrsFast-Ultra and Bowtie1.

- [mrsFast-Ultra](http://sfu-compbio.github.io/mrsfast/) version 3.3.11
- [Bowtie1](http://bowtie-bio.sourceforge.net/index.shtml) version 1.1.2

After installing mrsFastUltra and Bowtie1, you will need to build separate genome indicies for each. Follow the instructions included for each program for building these genome references.


## Getting Started



## Running the spike-in simulation

#### 1. Simulate reads

```bash
bsub < code/simulation_spikein/1_simulate_readlengths.sh
```

#### 2. Trim and filter the reads

```bash
bsub < code/simulation_spikein/2_trim_spikein.sh
```

#### 3. Align reads to the genome reference

Here you can choose one of the alignment strategies described in the paper:

- `code/simulation_spikein/bowtie2_spikein.sh` is the alignment strategy used in the [plethora](https://github.com/dpastling/plethora) pipeline, which finds the best alignment for each read.
- `code/simulation_spikein/bowtie_multiread_spikein.sh` is tailored for the multi-read correction strategy. It uses Bowtie1 to align reads to the genome, first looking for the best possible alignment. But if more than one valid alignment is found, it reports all tied alignments.
- `code/simulation_spikein/mrsfast_spikein.sh` uses mrsFast-Ultra to find all possible alignments for each read as described in the paper.
- `code/simulation_spikein/bowtie_total_spikein.sh` like the script above, but using Bowtie1 instead to find all possible alignments.

#### 4. Calculate coverage for each DUF1220 domain

```bash
bsub < code/simulation_spikein/make_bed_spikein.sh
```


## Running the baseline simulation

#### 1. Simulate reads

```bash
bsub < code/simulation_baseline/1_simulate_replicates.sh
```

#### 2. Trim and filter the reads

```bash
bsub < code/simulation_baseline/2_trim_replicates.sh
```

#### 3. Align reads to the genome reference

Here you can choose one of the alignment strategies described in the paper:

- `code/simulation_baseline/bowtie2_replicates.sh`
- `code/simulation_baseline/bowtie_multiread_replicates.sh`
- `code/simulation_baseline/mrsfast_replicates.sh`
- `code/simulation_baseline/bowtie_total_replicates.sh`

#### 4. Calculate coverage for each DUF1220 domain

```bash
bsub < code/simulation_baseline/make_bed_replicates.sh
```


## Analyzing the results in R

The code in the `code/analysis/` folder are analyzing the result files in R.


- `spikein_analysis.R` processes the result files for the spikein study
- `replicate_analysis.R` processes the result files for the baseline study

Behind the scenes, the following analysis scripts are run:

- `load_data.R` gathers all of the result files into a single data frame
- `normalize.R` normalizes the data depending which type of analysis was run


## Reproducing the figures in the paper

The code in the `code/figures/` folder was used to generate the figures for the paper.

- `heatmap_plots.R` was used to generate the heatmap plots for Figures 2 and 5.
- `combined_rmse_plot.R` was used to generate Figure 3.
- `combined_rmse_plot_by_gene.R` was used to generate Figure 4.




