## Author: The Winx Club: Andrea Fernandez Veloso, Adrian Perera Bonaño and Emma Serrano Perez.
## Email: fvandrea99@gmail.com, adrianpererads@gmail.com and emma.serrano.perez.99@gmail.com.

## This R script works with arguments passed from a bash script. 
## To run this from the Rstudio environment, the vector "arguments" is not necessary, and every
## file needs to be read individually. 

# Saving arguments in individual variables:

arguments <- commandArgs(trailingOnly = TRUE)
print(arguments)

narrow_peak <- arguments[1]
summits_bed <- arguments[2]
up1 <- as.numeric(arguments[3])
down1 <- as.numeric(arguments[4])
up2 <- as.numeric(arguments[5])
down2 <- as.numeric(arguments[6])
p_value_go <- as.numeric(arguments[7])
p_value_kegg <- as.numeric(arguments[8])
name_sample <- arguments[9]


## Installing packages:

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# General packages.

#if (!require("ChIPseeker")) BiocManager::install("ChIPseeker")
library("ChIPseeker")

#if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
library("clusterProfiler")

#if (!requireNamespace("ggplot2", quietly = TRUE))
#  install.packages("ggplot2")
library("ggplot2")

# Organism specific packages.

#if (!require(txdb)) BiocManager::install(txdb)
library("TxDb.Athaliana.BioMart.plantsmart28")
txdb <- TxDb.Athaliana.BioMart.plantsmart28

#if (!require(annoDb_name)) BiocManager::install(annoDb_name)
library("org.At.tair.db")


## Getting genetic information from organism

genes.organism <- as.data.frame(genes(txdb))
head(genes.organism)

## Defining the experiment universe

my.universe<-genes.organism$gene_id

## Reading narrow peaks

peak <- readPeakFile(peakfile = narrow_peak)
peak

## Reading summits

summit <- readPeakFile(peakfile = summits_bed)
summit 

## Peaks localization

name_covplot <- paste(name_sample, "_covplot.png", sep = "")
title_covplot <- paste(name_sample, "peak coverage")

png(file = name_covplot, width = 1536, height = 801, units = "px")
covplot(peak, weightCol="V5", title = title_covplot, ylab = "ChIP Peaks")
dev.off()

## Defining promoters

Promoters_peak <- getPromoters(TxDb=txdb, upstream=up1, downstream=down1)

Promoters_summit <- getPromoters(TxDb=txdb, upstream=up2, downstream=down2)

## Narrow peaks and summits annotation

name_annopie_narrow <- paste(name_sample, "_annopie_narrow.png", sep = "")
title_annopie_narrow <- paste(name_sample, "narrow peaks location")

peakAnno_peak <- annotatePeak(peak = peak, tssRegion = c(-up1,down1), TxDb = txdb, annoDb = "org.At.tair.db" )
png(file = name_annopie_narrow, width = 768, height = 400, units = "px")
plotAnnoPie(peakAnno_peak)
title(title_annopie_narrow, line = -5, cex = 15)
dev.off()

name_annopie_summit <- paste(name_sample, "_annopie_summit.png", sep = "")
title_annopie_summit <- paste(name_sample, "summits location")

peakAnno_summit <- annotatePeak(peak = summit, tssRegion = c(-up2,down2), TxDb = txdb, annoDb = "org.At.tair.db")
png(file = name_annopie_summit, width = 768, height = 400, units = "px")
plotAnnoPie(peakAnno_summit)
title(title_annopie_summit, line = -5, cex = 15)
dev.off()

df.annotation_peak <- as.data.frame(peakAnno_peak)
head(df.annotation_peak)

df.annotation_summit <- as.data.frame(peakAnno_summit)
head(df.annotation_summit)


## Regulome determination

## First, peaks located in promoter regions are extracted.

promoter.df.annotation_peak <- subset(df.annotation_peak, annotation == "Promoter")
head(promoter.df.annotation_peak)

## Getting names of genes located downstream the previously defined promoters, which is the regulome.

organism.regulome <- promoter.df.annotation_peak$geneId
length(organism.regulome)

name_regulome <- paste(name_sample, "_regulome.txt", sep = "")

write.table(organism.regulome, file = name_regulome, col.names = FALSE, row.names = FALSE)

## GSEA with clusterProfiler

## 1. Ontology terms enrichment

ego <- enrichGO(gene          = organism.regulome,
                universe      = my.universe,
                OrgDb         = org.At.tair.db,
                ont           = "ALL",
                pvalueCutoff  = p_value_go,
                keyType = 'TAIR')
print(ego)
summary(as.data.frame(ego))
head(ego)

name_barplot_go <- paste(name_sample, "_barplot_go.png", sep = "")
title_barplot_go <- paste(name_sample, "bar plot using GO enrichment analysis")
png(file = name_barplot_go, width = 826, height = 622, units = "px")
barplot(ego, showCategory = 30, title = title_barplot_go)
dev.off()

name_dotplot_go <- paste(name_sample, "_dotplot_go.png", sep = "")
title_dotplot_go <- paste(name_sample, "dot plot using GO enrichment analysis")
png(file = name_dotplot_go, width = 826, height = 622, units = "px")
dotplot(ego, showCategory = 30, title = title_dotplot_go)
dev.off()

name_cnetplot_go <- paste(name_sample, "_cnetplot_go.png", sep = "")
title_cnetplot_go <- paste(name_sample, "gene-concept network using GO enrichment analysis")
png(file = name_cnetplot_go, width = 1536, height = 801, units = "px")
cnetplot(ego, colorEdge = TRUE, node_label = "gene", cex_gene = 0.5) + ggtitle(title_cnetplot_go) + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## 2. KEGG Pathways enrichment

ekegg <- enrichKEGG(gene = organism.regulome, 
                           organism = "ath",
                           universe = my.universe, 
                           pvalueCutoff = p_value_kegg)

print(ekegg)
summary(as.data.frame(ekegg))
head(ekegg)

name_barplot_kegg <- paste(name_sample, "_barplot_kegg.png", sep = "")
title_barplot_kegg <- paste(name_sample, "bar plot using KEGG Pathways enrichment analysis")
png(file = name_barplot_kegg, width = 826, height = 622, units = "px")
barplot(ekegg, showCategory = 30, title = title_barplot_kegg)
dev.off()

name_dotplot_kegg <- paste(name_sample, "_dotplot_kegg.png", sep = "")
title_dotplot_kegg <- paste(name_sample, "dot plot using KEGG Pathways enrichment analysis")
png(file = name_dotplot_kegg, width = 826, height = 622, units = "px")
dotplot(ekegg, showCategory = 30, title = title_dotplot_kegg)
dev.off()

name_cnetplot_kegg <- paste(name_sample, "_cnetplot_kegg.png", sep = "")
title_cnetplot_kegg <- paste(name_sample, "gene-concept network using KEGG Pathways enrichment analysis")
png(file = name_cnetplot_kegg, width = 1536, height = 801, units = "px")
cnetplot(ekegg, colorEdge = TRUE, node_label = "gene", cex_gene = 0.5) + ggtitle(title_cnetplot_kegg) + theme(plot.title = element_text(hjust = 0.5))
dev.off()
