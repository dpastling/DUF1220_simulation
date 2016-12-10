options(stringsAsFactors = FALSE)

library(tidyr)
library(dplyr)

data.folder <- "results/replicates/mrsfast"

diploid.coverage <- 30
haploid.coverage <- diploid.coverage / 2


#####################################################################
# Load Data
#####################################################################
source("code/load_data.R")

# let's ignore the 36bp reads for now and consider those reads
# that have been trimmed down from 100bp to 36
X <- filter(X, ! grepl("_36bp_", file))
X <- mutate(X, file = gsub("_trim_", "_36bp_", file))

#####################################################################
# Gather Metadata
#####################################################################
X <- mutate(X, file = gsub("template_", "", file))
X <- X %>% 
	mutate(temp = file) %>%
	separate(temp, c("read.length", "coverage", "replicate", "pairing"), 
		sep = "_", extra = "drop", fill = "right") %>%
	mutate(read.length = gsub("bp", "", read.length)) %>%
	mutate(read.length = as.numeric(read.length)) %>%
	mutate(coverage = gsub("x", "", coverage)) %>%
	mutate(coverage = as.numeric(coverage))


#####################################################################
# Normalize Data
#####################################################################
source("code/normalize.R")


#####################################################################
# Calculate Error Rate
#####################################################################
rmse <- function(p, m) sqrt(mean((m - p)^2))

error_for_group <- function(df, stat = "rmse", summarise.by.clade = FALSE)
{
        df <- group_by(df, read.length, coverage, pairing, replicate, add = TRUE)
	if (summarise.by.clade) df <- group_by(df, measured.clade, add = TRUE)
        df <- summarise(df, read.depth = sum(read.depth), expected = sum(coverage))
        df <- group_by(df, read.length, coverage, pairing, add = FALSE)
	if (summarise.by.clade) df <- group_by(df, measured.clade, add = TRUE)
	if (stat == "rmse")
	{
        	df <- summarise(df, rmse = rmse(read.depth, expected))
	} else if (stat == "sd")
	{
		df <- summarise(df, sd = sd(read.depth))
	}
	return(df)
}

domain.rmse <- group_by(X, file, measured.domain) %>% error_for_group() %>% rename(domain.rmse = rmse)
gene.rmse   <- group_by(X, file, measured.gene.id) %>% error_for_group() %>% rename(gene.rmse = rmse)
group.rmse  <- group_by(X, file, measured.group.id) %>% error_for_group() %>% rename(group.rmse = rmse)
clade.rmse  <- group_by(X, file, measured.clade) %>% error_for_group() %>% rename(clade.rmse = rmse)
total.rmse  <- group_by(X, file) %>% error_for_group() %>% rename(total.rmse = rmse)

X.rmse <- left_join(domain.rmse, gene.rmse) %>% 
	  left_join(., group.rmse) %>% 
	  left_join(., clade.rmse) %>% 
	  left_join(., total.rmse)

# split out error rate by clade
domain.rmse <- group_by(X, file, measured.domain) %>% error_for_group(summarise.by.clade = TRUE) %>% rename(domain.rmse = rmse)
gene.rmse   <- group_by(X, file, measured.gene.id) %>% error_for_group(summarise.by.clade = TRUE) %>% rename(gene.rmse = rmse)
group.rmse  <- group_by(X, file, measured.group.id) %>% error_for_group(summarise.by.clade = TRUE) %>% rename(group.rmse = rmse)
clade.rmse  <- group_by(X, file) %>% error_for_group(summarise.by.clade = TRUE) %>% rename(clade.rmse = rmse)
total.rmse  <- group_by(X, file) %>% error_for_group() %>% rename(total.rmse = rmse)

X.rmse.clade <- left_join(domain.rmse, gene.rmse) %>% 
		left_join(., group.rmse) %>% 
		left_join(., clade.rmse) %>% 
		left_join(., total.rmse)


domain.sd <- group_by(X, file, measured.domain) %>% error_for_group(stat = "sd") %>% rename(domain.sd = sd)
gene.sd   <- group_by(X, file, measured.gene.id) %>% error_for_group(stat = "sd") %>% rename(gene.sd = sd)
group.sd  <- group_by(X, file, measured.group.id) %>% error_for_group(stat = "sd") %>% rename(group.sd = sd)
clade.sd  <- group_by(X, file, measured.clade) %>% error_for_group(stat = "sd") %>% rename(clade.sd = sd)
total.sd  <- group_by(X, file) %>% error_for_group(stat = "sd") %>% rename(total.sd = sd)

X.sd <- left_join(domain.sd, gene.sd) %>% 
	  left_join(., group.sd) %>% 
	  left_join(., clade.sd) %>% 
	  left_join(., total.sd)

# split out error rate by clade
domain.sd <- group_by(X, file, measured.domain) %>% error_for_group(summarise.by.clade = TRUE, stat = "sd") %>% rename(domain.sd = sd)
gene.sd   <- group_by(X, file, measured.gene.id) %>% error_for_group(summarise.by.clade = TRUE, stat = "sd") %>% rename(gene.sd = sd)
group.sd  <- group_by(X, file, measured.group.id) %>% error_for_group(summarise.by.clade = TRUE, stat = "sd") %>% rename(group.sd = sd)
clade.sd  <- group_by(X, file) %>% error_for_group(summarise.by.clade = TRUE, stat = "sd") %>% rename(clade.sd = sd)
total.sd  <- group_by(X, file) %>% error_for_group(stat = "sd") %>% rename(total.sd = sd)

X.sd.clade <- left_join(domain.sd, gene.sd) %>% 
		left_join(., group.sd) %>% 
		left_join(., clade.sd) %>% 
		left_join(., total.sd)


#####################################################################
# Write Data
#####################################################################
write.table(X.rmse, file = paste(data.folder, "replicate_rmse.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)
write.table(X.rmse.clade, file = paste(data.folder, "replicate_rmse_by_clade.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)
write.table(X.sd, file = paste(data.folder, "replicate_sd.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)
write.table(X.sd.clade, file = paste(data.folder, "replicate_sd_by_clade.txt", sep = "/"), sep = "\t", row.names = FALSE, quote = F)




