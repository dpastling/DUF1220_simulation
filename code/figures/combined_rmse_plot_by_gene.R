options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
 
#X <- read.delim("data/simulation_oct2015/X_gene_rmse_clade_groups.txt")
#X <- read.delim("data/simulation_aug2016/X_gene_rmse.txt")
#X <- read.delim("data/simulation_2016-10-15/X_gene_rmse.txt")
X <- read.delim("results/spikein/normal/gene_rmse.txt")
X <- filter(X, read.length == 100, pairing == "paired")

X.sum <- gather(X, Summarization, RMSE, domain.rmse:clade.rmse)
X.sum <- mutate(X.sum, Summarization = gsub(".rmse", "", Summarization))
X.sum <- mutate(X.sum, Summarization = factor(Summarization, levels = c("domain", "gene", "group", "clade")))

gene.names <- unique(X.sum[["spiked.gene"]])
gene.index <- gsub("NBPF", "", gene.names)
gene.index <- gsub("P", "", gene.index)
gene.index <- gsub("L", ".1", gene.index)
gene.index <- as.numeric(gene.index)
gene.names <- gene.names[order(gene.index)]
manual.order <- c("NBPF1", "NBPF1L", "NBPF4", "NBPF5P", "NBPF6", "NBPF10", "NBPF14", "NBPF19", "NBPF20", "NBPF15", "NBPF8", "NBPF26", "NBPF12", "NBPF11", "NBPF2P", "NBPF3", "NBPF7", "NBPF9", "NBPF13P", "NBPF17P", "NBPF21P", "NBPF25P")

X.sum <- mutate(X.sum, spiked.gene = factor(spiked.gene, levels = manual.order))

pdf(file = "plots/RMSE_by_gene.pdf", width = 12, height = 6)
p <- ggplot(X.sum, aes(spiked.gene, RMSE)) + geom_bar(stat = "identity")
p <- p + facet_grid(spiked.clade ~ Summarization) 
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p <- p + xlab("")
print(p)
dev.off()


