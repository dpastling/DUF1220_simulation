
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
	X <- X %>%
		group_by(file, measured.clade) %>%
		mutate(read.depth = read.depth / n()) %>%
		ungroup()
}



