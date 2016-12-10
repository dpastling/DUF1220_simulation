
options(stringsAsFactors = FALSE)

library(tidyr)
library(dplyr)

stopifnot(exists("X"))
required.columns <- c("read.depth", "read.length", "file", "measured.clade", "pairing")
stopifnot(all(required.columns %in% colnames(X)))
stopifnot(exists("data.folder"))

if (! exists("method") & grepl("(mrsfast)|(total)", data.folder))
{
	method <- "total"
} else if (! exists("method")) {
	method <- "normal"
}

#####################################################################
# Normalize Data
#####################################################################
# for single end reads, our theoretical coverage is lower because the number
# of simulated reads was based on the insert size. instead divide by the adjusted
# coverage.
X <- X %>% 
	mutate(insert.size = 350) %>%
	mutate(insert.size = ifelse(read.length == 300, 700, insert.size)) %>%
	mutate(insert.size = ifelse(read.length == 600, 1400, insert.size))

X <- mutate(X, read.depth = read.depth / haploid.coverage)
if (grepl("(mrsfast)|(total_trim)", data.folder))
{
	X <- X %>%
        	mutate(read.depth = ifelse(pairing == "single", read.depth * (insert.size / 36), read.depth))
} else {
	X <- X %>%
		mutate(read.depth = ifelse(pairing == "single", read.depth * (insert.size / read.length), read.depth))
}

if (method == "total")
{
	# normalize align all strategy by the number of clade domains
	# note this does not effect the percent off.target error rate
#	outlier.genes <- c("NBPF13P", "NBPF17P", "NBPF3", "NBPF2P", "NBPF7", "NBPF21P")
#	X <- filter(X, ! measured.gene %in% outlier.genes)
#	mrsfast.groupings <- c(
#		"NBPF(4|(5P)|6)_",
#		"NBPF(3|(17P)|(2P))_",
#		"NBPF13P",
#		"NBPF12P"
#	)
#	X <- mutate(X, measured.group.id = "majority_")
#	for (group in mrsfast.groupings)
#	{
#		X <- mutate(X, measured.group.id = ifelse(grepl(group, measured.domain), group, measured.group.id))
#	}
#	X <- mutate(X, measured.group.id = paste0(measured.group.id, measured.clade))
#	X <- mutate(X, measured.group.id = measured.clade)
#	X <- mutate(X, measured.group.id = gsub("HLS[123]", "HLS", measured.group.id))
	X <- X %>%
		group_by(file, measured.clade) %>%
#		group_by(file, measured.group.id) %>%
		mutate(read.depth = read.depth / n()) %>%
#		mutate(read.depth = read.depth / sum(coverage)) %>%
		ungroup()
}



