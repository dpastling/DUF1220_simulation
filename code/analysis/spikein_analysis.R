
options(stringsAsFactors = FALSE)

library(tidyr)
library(dplyr)

#data.folder <- "results/spikein/mrsfast"
#data.folder <- "results/spikein/multiread"
#data.folder <- "results/spikein/total"
data.folder <- "results/spikein/normal"

diploid.coverage <- 30
haploid.coverage <- diploid.coverage / 2

method <- "normal"

#####################################################################
# Load Data
#####################################################################
source("code/analysis/load_data.R")

# ignore the trimmed files
X <- filter(X, ! grepl("trim", file))

#####################################################################
# Gather Metadata
#####################################################################
X <- X %>%
	mutate(temp = file) %>% 
	separate(temp, c("spiked.gene", "spiked.clade", "spiked.clade.index", "read.length", "coverage", "pairing"), 
		sep = "_", extra = "drop", fill = "right") %>%
	mutate(read.length = gsub("bp", "", read.length)) %>%
	mutate(read.length = as.numeric(read.length)) %>%
	mutate(coverage = gsub("x", "", coverage)) %>%
	mutate(coverage = as.numeric(coverage)) %>%
	mutate(spiked.domain = paste(spiked.gene, spiked.clade, spiked.clade.index, sep = "_"))

# Categorize the spiked domain
X <- mutate(X, spiked.gene.id  = paste(spiked.gene, spiked.clade, sep = "_"))
X <- mutate(X, spiked.group.id = spiked.gene.id)
for (group in groupings)
{
	X <- mutate(X, spiked.group.id = gsub(group, group, spiked.group.id))
}


#####################################################################
# Normalize Data
#####################################################################
source("code/analysis/normalize.R")


#####################################################################
# Calculate Error Rate
#####################################################################
rmse <- function(p, m) sqrt(mean((m - p)^2))

error.rate <- X %>%
	group_by(spiked.gene, spiked.clade, spiked.clade.index, read.length, pairing, coverage) %>%
	summarise(
		domain.on.target  = sum(read.depth[measured.domain == spiked.domain]) / sum(read.depth), 
		domain.off.target = sum(read.depth[measured.domain != spiked.domain]) / sum(read.depth),
		gene.on.target    = sum(read.depth[measured.gene.id == spiked.gene.id]) / sum(read.depth), 
		gene.off.target  = sum(read.depth[measured.gene.id != spiked.gene.id]) / sum(read.depth),
		group.on.target  = sum(read.depth[measured.group.id == spiked.group.id]) / sum(read.depth), 
		group.off.target = sum(read.depth[measured.group.id != spiked.group.id]) / sum(read.depth),
		clade.on.target  = sum(read.depth[measured.clade == spiked.clade]) / sum(read.depth), 
		clade.off.target = sum(read.depth[measured.clade != spiked.clade]) / sum(read.depth)
	) %>%
	group_by(read.length, pairing, spiked.clade) %>%
	summarise(
		domain.on.target  = mean(domain.on.target), 
		domain.off.target = mean(domain.off.target),
		gene.on.target    = mean(gene.on.target), 
		gene.off.target   = mean(gene.off.target),
		group.on.target   = mean(group.on.target), 
		group.off.target  = mean(group.off.target),
		clade.on.target   = mean(clade.on.target),
		clade.off.target  = mean(clade.off.target)  
	)


error.rate.predicted <- X %>%
	group_by(spiked.gene, spiked.clade, spiked.clade.index, read.length, pairing, coverage) %>%
	summarise(
		domain.on.target  = sum(read.depth[measured.domain == spiked.domain]) / sum(coverage), 
		domain.off.target = sum(read.depth[measured.domain != spiked.domain]) / sum(coverage),
		gene.on.target    = sum(read.depth[measured.gene.id == spiked.gene.id]) / sum(coverage), 
		gene.off.target  = sum(read.depth[measured.gene.id != spiked.gene.id]) / sum(coverage),
		group.on.target  = sum(read.depth[measured.group.id == spiked.group.id]) / sum(coverage), 
		group.off.target = sum(read.depth[measured.group.id != spiked.group.id]) / sum(coverage),
		clade.on.target  = sum(read.depth[measured.clade == spiked.clade]) / sum(coverage), 
		clade.off.target = sum(read.depth[measured.clade != spiked.clade]) / sum(coverage)
	) %>%
	group_by(read.length, pairing, spiked.clade) %>%
	summarise(
		domain.on.target  = mean(domain.on.target), 
		domain.off.target = mean(domain.off.target),
		gene.on.target    = mean(gene.on.target), 
		gene.off.target   = mean(gene.off.target),
		group.on.target   = mean(group.on.target), 
		group.off.target  = mean(group.off.target),
		clade.on.target   = mean(clade.on.target),
		clade.off.target  = mean(clade.off.target)  
	)


# First summarise the read.depth by each category
# we will be recycling the X.rmse object
X.rmse <- X %>%
group_by(spiked.gene, spiked.clade, spiked.clade.index, read.length, pairing, coverage) %>%
summarise(
	measured.domain  = sum(read.depth[measured.domain == spiked.domain]), 
	measured.gene    = sum(read.depth[measured.gene.id == spiked.gene.id]), 
	measured.group   = sum(read.depth[measured.group.id == spiked.group.id]), 
	measured.clade   = sum(read.depth[measured.clade == spiked.clade]),
	total.coverage   = sum(coverage)
)

# Calculate the RMSE by gene. This table will be used for the gene barplot
X.gene.rmse <- X.rmse %>%
	group_by(read.length, pairing, spiked.clade, spiked.gene) %>%
	summarise(
		domain.rmse = rmse(measured.domain, coverage),
		gene.rmse   = rmse(measured.gene, coverage),
		group.rmse  = rmse(measured.group, coverage),
		clade.rmse  = rmse(measured.clade, coverage)
	)

# Calculate the RMSE by clase. This will be used for the line plots
X.rmse <- X.rmse %>%
	group_by(read.length, pairing, spiked.clade) %>%
	summarise(
		domain.rmse = rmse(measured.domain, coverage),
		gene.rmse   = rmse(measured.gene, coverage),
		group.rmse  = rmse(measured.group, coverage),
		clade.rmse  = rmse(measured.clade, coverage)
	)



#####################################################################
# Write Data
#####################################################################
write.table(error.rate, file = paste(data.folder, "error_rate.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)
write.table(error.rate.predicted, file = paste(data.folder, "error_rate_predicted.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)
write.table(X.rmse, file = paste(data.folder, "rmse.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)
write.table(X.gene.rmse, file = paste(data.folder, "gene_rmse.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)


