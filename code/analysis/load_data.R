

options(stringsAsFactors = FALSE)

library(tidyr)
library(dplyr)

stopifnot(exists("data.folder"))

if (! exists("file.pattern"))
{
	file.pattern <- "_read_depth.bed"
}

#####################################################################
# Load Data
#####################################################################
files <- list.files(data.folder, pattern = file.pattern, full.names = TRUE, recursive = TRUE)

# Ignore the combined files
files <- grep("combined_", files, invert = TRUE, value = TRUE)

X <- lapply(files, read.delim, header = FALSE, colClasses = c("character", "numeric"))
files <- gsub("^.+?/([^/]+)$", "\\1", files)
files <- gsub(file.pattern, "", files)
names(X) <- files

X <- bind_rows(X, .id = "file")
colnames(X) <- c("file", "domain", "read.depth")


#####################################################################
# Gather Metadata
#####################################################################
X <- X %>% 
	mutate(temp = domain) %>% 
	separate(temp, c("measured.gene", "measured.clade", "measured.clade.index"), sep = "_") %>%
	rename(measured.domain = domain)


#####################################################################
# Setup categories and groupings
#####################################################################
groupings <- c(
	"NBPF(4|(5P)|6)_", 
	"NBPF(1|(1L))_", 
	"NBPF((10)|(14)|(19))_CON", 
	"NBPF((10)|(14)|(19)|(20))_HLS[12]", 
	"NBPF((10)|(14)|(19)|(20)|(15))_HLS3", 
	"NBPF((11)|(12))_", 
	"NBPF((8)|(26))_"
)

# Categorize the measured domains
X <- mutate(X, measured.gene.id  = paste(measured.gene, measured.clade, sep = "_"))
X <- mutate(X, measured.group.id = measured.gene.id)
for (group in groupings)
{
	X <- mutate(X, measured.group.id = gsub(group, group, measured.group.id))
}

