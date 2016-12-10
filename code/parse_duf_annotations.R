
# this script imports a table of raw DUF1220 coordinates and prepares a bed file for read depth processing

options(stringsAsFactors = FALSE)

library(dplyr)

# This is the most we will extend the domains
max.extension <- 250

bed <- read.delim("reference/annotation-clade-based-numbering-full-domains-2016-11-29.bed", header = FALSE)
colnames(bed) <- c("chr", "pos1", "pos2", "domain", "gene", "strand", "length")

bed <- bed[order(bed$chr, bed$pos1), ]

bed <- mutate(bed, domain = gsub("AC240274.1", "NBPF1L", domain))

intron.size <- c(bed$pos1, 1e9) - c(0, bed$pos2)

# if more than one chromosome is present, there will be strange negative gaps (ch1:100,000,000 vs chr2:20,000)
chromosomal.breaks <- which(bed$chr[2:nrow(bed)] != bed$chr[1:(nrow(bed) - 1)])
intron.size[chromosomal.breaks + 1] <- 1e5

# Fragments are typically up to 400 bp in length, so we want to create a buffer 
# to eliminate fragments that align near the midpoint between domains
extension <- intron.size - 800
extension <- extension / 2
extension <- floor(extension)
extension[extension < 0] <- 0
extension[extension > max.extension] <- max.extension

bed <- mutate(bed, pos1 = pos1 - extension[1:(length(extension) - 1)] )
bed <- mutate(bed, pos2 = pos2 + extension[2:length(extension)])
bed <- mutate(bed, length = pos2 - pos1)
bed <- mutate(bed, score = 255)
bed <- select(bed, chr, pos1, pos2, domain, score, strand)


write.table(bed, file = "reference/DUF1220_full_domains_hg38_v2.2.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

bed <- filter(bed, grepl("((CON)|(HLS))[1-3]", domain))
bed <- filter(bed, ! grepl("HLS\\d/CON", domain))
bed <- filter(bed, ! grepl("CON\\d/HLS", domain))
bed <- filter(bed, ! grepl("/", domain))
bed <- filter(bed, ! grepl("long", domain, ignore.case = TRUE))
bed <- filter(bed, ! grepl("short", domain, ignore.case = TRUE))

write.table(bed, file = "reference/DUF1220_canonical_domains_hg38_v2.2.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

