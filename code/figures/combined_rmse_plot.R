options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
 
#X <- read.delim("data/simulation_oct2015/X_rmse.txt")
#X <- read.delim("data/simulation_oct2015/X_rmse_clade_groups.txt")
#X <- read.delim("data/simulation_june2016/padded_domains/diagnostic/X_rmse.txt")
#X <- read.delim("data/simulation_aug2016/X_rmse.txt")
#X <- read.delim("data/simulation_2016-10-15/X_rmse.txt")
X <- read.delim("results/spikein/normal/rmse.txt")

X <- filter(X, read.length < 600)

X.sum <- X %>% gather(summarization, RMSE, domain.rmse:clade.rmse)
X.sum <- mutate(X.sum, summarization = gsub(".rmse", "", summarization))
X.sum <- mutate(X.sum, summarization = factor(summarization, levels = c("domain", "gene", "group", "clade")))


pdf(file = "plots/RMSE_by_read_length.pdf", width = 12, height = 6)
p <- ggplot(X.sum, aes(read.length, RMSE, color = pairing)) + geom_point() + geom_line()
p <- p + facet_grid(summarization ~ spiked.clade)
p <- p + xlab("Read Length (in bp)")
p <- p + theme_bw()
print(p)
dev.off()

