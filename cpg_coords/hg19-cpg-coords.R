library(BSgenome)
library(GenomicRanges)
library(tidyverse)
library(data.table)

CHR_NUMS <- as.character(c(1:22, 'X', 'Y'))
CHR_NAMES <- paste0('chr', CHR_NUMS)

# CpG coordinates in reference genome
genome <- BSgenome.Hsapiens.UCSC.hg19
hg19.cpg.coords <- map_dfr(
  CHR_NAMES,
  ~ matchPattern('CG', genome[[.x]]) %>% as.data.frame %>% mutate(chr = .x) %>% select(chr, start, end)
)
hg19.cpg.coords <- GenomicRanges::makeGRangesFromDataFrame(hg19.cpg.coords)

save(hg19.cpg.coords, 'cpg_coords/hg19.cpg.coords.rda')


# `hg19.cpg.coords` is An object of class \code{GRanges} that contains genomic
# information of CpG positions in hg19. This dataset represents a GRanges object
# that contains the collection of CpGs in chr1-chr22, chrX and chrY. Each row is a CpG.
