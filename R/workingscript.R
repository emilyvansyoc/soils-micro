## DESEQ2

if(!require(DESeq2)) {
  
  BiocManager::install("DESeq2")
  require(DESeq2)
}

require(tidyverse)
require(phyloseq)
require(pairwiseAdonis)
require(RColorBrewer)
# load phyloseq objects and metadata
if(!require(repmis)) {
  install.packages("repmis")
  require(repmis)
}
# read in data
source_data("https://github.com/EmilyB17/soils-micro/blob/master/data/phyloseq16S_ITS.RData?raw=true")

set.seed(123)

# convert to deseq2 format

test <- phyloseq_to_deseq2(physeq = psITS18,
                           design = ~ Treatment)

# perform DeSeq testing
de <- DESeq(test, test = "Wald", fitType = "parametric")

# get results
res = results(de, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psITS18)[rownames(sigtab), ], "matrix"))
head(sigtab)

# plot for results
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
