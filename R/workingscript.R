## DESEQ2 at the genus level

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

## ---- Fungal: pre-process ----

## pre-process phyloseq object to replace DNA strings with arbitrary ASV name
dna <- Biostrings::DNAStringSet(taxa_names(psITS18))
names(dna) <- taxa_names(psITS18)
ps <- merge_phyloseq(psITS18, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(psITS18)))

# first agglomerate taxa at the genus level
psgen <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

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

### ---- Fungal: bulk vs rhizo ----

# our "null" comparison has to be the first factor for comparisons
sample_data(psgen)$sample_type <- relevel(sample_data(psgen)$sample_type, ref = "bulk")
str(sample_data(psgen)$sample_type) # bulk is level 1 so will be the "Null" comparison

# convert phyloseq to DESeq2 object
desob <- phyloseq_to_deseq2(physeq = psgen,
                           design = ~ sample_type)

# perform DeSeq testing
de <- DESeq(desob, test = "Wald", fitType = "parametric")

# get results
res = results(de, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psgen)[rownames(sigtab), ], "matrix"))
head(sigtab)

# sort by log2fold change
sigtab <- sigtab[order(sigtab$log2FoldChange), ]
# write table of results
write.table(sigtab, "./data/deseq2-its-bulkvsrhizo.txt", sep = "\t", row.names= FALSE)

## ---- Fungal: treatments ----

# first, remove PRE time so we are only comparing samples with grazing impact
nopre <- prune_samples(sample_data(psgen)$GrazeTime != "PRE", psgen)

# make the NG treatment the NULL level
sample_data(nopre)$Treatment <- relevel(sample_data(nopre)$Treatment, ref = "NO")

# convert to deseq2 format

detrt <- phyloseq_to_deseq2(physeq = nopre,
                           design = ~ Treatment)

# perform DeSeq testing
de <- DESeq(detrt, test = "Wald", fitType = "parametric")

# because we did multiple tests, get each pairwise comparison separately
resultsNames(de)
# store results in a dataframe
resHI <- as.data.frame(DESeq2::results(de, format = "DataFrame", 
                                     name = "Treatment_HI_vs_NO",
                                     cooksCutoff = FALSE))
resLO <- as.data.frame(DESeq2::results(de, format = "DataFrame", 
                                       name = "Treatment_LO_vs_NO",
                                       cooksCutoff = FALSE))
# get only significant results
alpha = 0.05
sigtabHI = resHI[which(resHI$padj < alpha), ]
sigtabHI = cbind(as(sigtabHI, "data.frame"), as(tax_table(nopre)[rownames(sigtabHI), ], "matrix"))
head(sigtabHI)
# sort by log2fold change
sigtabHI <- sigtabHI[order(sigtabHI$log2FoldChange),]
# write table of results
write.table(sigtabHI, "./data/deseq2-its-HIvsNO.txt", sep = "\t", row.names= FALSE)


sigtabLO = resLO[which(resLO$padj < alpha), ]
sigtabLO = cbind(as(sigtabLO, "data.frame"), as(tax_table(nopre)[rownames(sigtabLO), ], "matrix"))
head(sigtabLO)
# sort by log2fold change
sigtabLO <- sigtabLO[order(sigtabLO$log2FoldChange), ]
# write table of results
write.table(sigtabLO, "./data/deseq2-its-LOvsNO.txt", sep = "\t", row.names= FALSE)

## ---- Bacteria: pre-process ----

## pre-process phyloseq object to replace DNA strings with arbitrary ASV name
dna <- Biostrings::DNAStringSet(taxa_names(ps16S18))
names(dna) <- taxa_names(ps16S18)
ps <- merge_phyloseq(ps16S18, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps16S18)))

# first agglomerate taxa at the genus level
psgen <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

### ---- Fungal: bulk vs rhizo ----

# our "null" comparison has to be the first factor for comparisons
sample_data(psgen)$sample_type <- relevel(sample_data(psgen)$sample_type, ref = "bulk")
str(sample_data(psgen)$sample_type) # bulk is level 1 so will be the "Null" comparison

# convert phyloseq to DESeq2 object
desob <- phyloseq_to_deseq2(physeq = psgen,
                            design = ~ sample_type)

# perform DeSeq testing
de <- DESeq(desob, test = "Wald", fitType = "parametric")

# get results
res = results(de, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psgen)[rownames(sigtab), ], "matrix"))
head(sigtab)

# sort by log2fold change
sigtab <- sigtab[order(sigtab$log2FoldChange), ]
# write table of results
write.table(sigtab, "./data/deseq2-16S-bulkvsrhizo.txt", sep = "\t", row.names= FALSE)

## ---- Bacteria: treatments ----

# first, remove PRE time so we are only comparing samples with grazing impact
nopre <- prune_samples(sample_data(psgen)$GrazeTime != "PRE", psgen)

# make the NG treatment the NULL level
sample_data(nopre)$Treatment <- relevel(sample_data(nopre)$Treatment, ref = "NO")

# convert to deseq2 format

detrt <- phyloseq_to_deseq2(physeq = nopre,
                            design = ~ Treatment)

# perform DeSeq testing
de <- DESeq(detrt, test = "Wald", fitType = "parametric")

# because we did multiple tests, get each pairwise comparison separately
resultsNames(de)
# store results in a dataframe
resHI <- as.data.frame(DESeq2::results(de, format = "DataFrame", 
                                       name = "Treatment_HI_vs_NO",
                                       cooksCutoff = FALSE))
resLO <- as.data.frame(DESeq2::results(de, format = "DataFrame", 
                                       name = "Treatment_LO_vs_NO",
                                       cooksCutoff = FALSE))
# get only significant results
alpha = 0.05
sigtabHI = resHI[which(resHI$padj < alpha), ]
sigtabHI = cbind(as(sigtabHI, "data.frame"), as(tax_table(nopre)[rownames(sigtabHI), ], "matrix"))
head(sigtabHI)
# sort by log2fold change
sigtabHI <- sigtabHI[order(sigtabHI$log2FoldChange),]
# write table of results
write.table(sigtabHI, "./data/deseq2-16S-HIvsNO.txt", sep = "\t", row.names= FALSE)


sigtabLO = resLO[which(resLO$padj < alpha), ] # NO SIGNIFICANCE
#sigtabLO = cbind(as(sigtabLO, "data.frame"), as(tax_table(nopre)[rownames(sigtabLO), ], "matrix"))
#head(sigtabLO)
# sort by log2fold change
#sigtabLO <- sigtabLO[order(sigtabLO$log2FoldChange), ]
# write table of results
#write.table(sigtabLO, "./data/deseq2-its-LOvsNO.txt", sep = "\t", row.names= FALSE)

## NO SIGNIFICANT DIFFERENCE BETWEEN LO AND NO IN 16S