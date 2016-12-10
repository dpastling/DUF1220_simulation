#!/usr/bin/env R

options(stringsAsFactors = FALSE)

library(gplots)
library(tidyr)
library(dplyr)

arrayRed <- function(dataSet, threshold = 20) {
    n <- length(dataSet)
    dataRange <- seq(length = n, from = min(dataSet, na.rm = T), to = max(dataSet, na.rm = T))
    colorRange <- rep(rgb(1, 1, 1), times = n)
    for (i in 1:n) {
        if (is.na(dataRange[i]) | dataRange[i] < 0)
        { 
            colorRange[i] <- rgb(1, 1, 1)
        } else if (dataRange[i] >= threshold) {
            colorRange[i] <- rgb(1, 0, 0)
        } else {
            colorRange[i] <- rgb(1, 1 - (dataRange[i] / threshold ), 1 - (dataRange[i] / threshold ))    
        }
    }
    return(colorRange)
}

arrayBlack <- function(dataSet, threshold = 20) {
    n <- length(dataSet)
    dataRange <- seq(length = n, from = min(dataSet, na.rm = T), to = max(dataSet, na.rm = T))
    colorRange <- rep(rgb(1, 1, 1), times = n)
    for (i in 1:n) {
        if (is.na(dataRange[i]) | dataRange[i] < 0)
        { 
            colorRange[i] <- rgb(1, 1, 1)
        } else if (dataRange[i] >= threshold) {
            colorRange[i] <- rgb(0, 0, 0)
        } else {
            colorRange[i] <- rgb(1 - (dataRange[i] / threshold ), 1 - (dataRange[i] / threshold ), 1 - (dataRange[i] / threshold ))    
        }
    }
    return(colorRange)
}


arrayHue <- function(dataSet, hue = 1, threshold = 20) {
    n <- length(dataSet)
    dataRange <- seq(length = n, from = min(dataSet, na.rm = T), to = max(dataSet, na.rm = T))
    colorRange <- rep(hsv(0, 0, 0), times = n)
    for (i in 1:n) {
        if (is.na(dataRange[i]) | dataRange[i] < 0)
        { 
            colorRange[i] <- hsv(0, 0, 0)
        } else if (dataRange[i] >= threshold) {
            colorRange[i] <- hsv(hue, 1, 1)
        } else {
            colorRange[i] <- hsv(hue, dataRange[i] / threshold, 1)    
        }
    }
    return(colorRange)
}


data.folder <- "results/spikein/mrsfast"
plot.label <- "mrsfast"

#data.folder <- "results/spikein/normal"
#plot.label <- "normal"

file.pattern <- "_read_depth.bed"
diploid.coverage <- 30
haploid.coverage <- diploid.coverage / 2

method <- "normal"

#####################################################################
# Load Data
#####################################################################
source("code/analysis/load_data.R")

X <- filter(X, grepl("100bp_10x", file))

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


####################################################################
# Normalize Data
#####################################################################
source("code/analysis/normalize.R")


####################################################################
# Summarise data for heatmap
#####################################################################
X <- filter(X, spiked.clade == measured.clade)

X <- X %>%
	group_by(file) %>%
	mutate(normalized.read.depth = read.depth / coverage)

if (grepl("mrsfast", data.folder))
{
	X <- X %>%
		group_by(spiked.domain, spiked.gene.id, measured.gene.id) %>%
		summarise(read.depth = mean(normalized.read.depth))
} else {
	X <- X %>%
		group_by(spiked.domain, spiked.gene.id, measured.gene.id) %>%
		summarise(read.depth = sum(normalized.read.depth))
}

X <- X %>%
	group_by(spiked.gene.id, measured.gene.id) %>%
	summarise(read.depth = mean(read.depth))


####################################################################
# Generate Plots for each clade
#####################################################################
for (domain in c(paste0("CON", 1:3), paste0("HLS", 1:3)))
{
    X.clade <- X %>%
	ungroup() %>%
	filter(grepl(domain, spiked.gene.id)) %>%
	mutate(spiked.gene.id   = gsub("^(NBPF[^_]+)_.+$", "\\1", spiked.gene.id)) %>%
	mutate(measured.gene.id = gsub("^(NBPF[^_]+)_.+$", "\\1", measured.gene.id)) %>%
	spread(measured.gene.id, read.depth, fill = 0)

    X.clade <- as.data.frame(X.clade)
    rownames(X.clade) <- X.clade[["spiked.gene.id"]]
    X.clade <- select(X.clade, -spiked.gene.id)
    X.clade <- as.matrix(X.clade)
   
    stopifnot(all(rownames(X.clade) %in% colnames(X.clade))) 
    stopifnot(all(colnames(X.clade) %in% rownames(X.clade))) 

    X.clade <- X.clade[, rownames(X.clade)]

    #roworder <- c("NBPF10", "NBPF14", "NBPF19", "NBPF15", "NBPF25P", "NBPF26", "NBPF4", "NBPF5P", "NBPF6", "NBPF8", "NBPF9", "NBPF11", "NBPF12", "NBPF1", "Unknown", "NBPF20", "NBPF3", "NBPF7", "NBPF13P", "NBPF17P", "NBPF18P", "NBPF21P")
    #roworder <- roworder[roworder %in% colnames(X)]
    #X <- X[roworder, roworder]
    
    dd.col <- as.dendrogram(hclust(dist(X.clade)))
    col.ord <- order.dendrogram(dd.col)

    X.clade <- X.clade[col.ord, col.ord]
    
    # color.range <- c(202, 188, 360, 360, 34, 120)
    # color.range <- color.range / 360
    # names(color.range) <- c(paste0("CON", 1:3), paste0("HLS", 1:3))
    # colorScheme <- arrayHue(X, hue = color.range[domain], threshold = 0.9)
    # if (domain == "CON3")
    # {
    #     colorScheme <- arrayBlack(X, threshold = 0.7)
    # }

    colorScheme = arrayRed(X.clade, threshold = 1)
    
    
    pdf(file = paste0("plots/heatmap_", plot.label, "_", domain, ".pdf"))
    heatmap.2(
	X.clade, 
	Rowv = FALSE, 
	Colv = FALSE, 
	scale = "none", 
	key = FALSE, 
	trace = "none", 
	dendrogram = "none", 
	col = colorScheme, 
	lwid = c(0.25, 5.5), 
	lhei = c(0.25, 5.5), 
	margins = c(7.5, 7.5), 
	cexRow = 1.7, 
	cexCol = 1.7, 
	adjCol = c(NA, 0.6)
	)
    dev.off()
    #title(main = domain)
}


