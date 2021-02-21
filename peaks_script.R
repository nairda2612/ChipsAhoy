## Author: The Winx Club
## Email: winx@alum.alfea.es
## Date: January 30

##This script works with arguments in order to submit it with R script. 
##If you want to run it in Rstudio, don't use arguments. ARREGLAR!!!

# Arguments:
arguments <- commandArgs(trailingOnly = TRUE)

narrow_peak <- arguments[1] #Peaks file
summits_bed <- arguments[2] #Summits file
up1 <- as.numeric(arguments[3])
down1 <- as.numeric(arguments[4])
up2 <- as.numeric(arguments[5])
down2 <- as.numeric(arguments[6])
p_value_go <- as.numeric(arguments[7])
p_value_kegg <- as.numeric(arguments[8])
title <- arguments[9]

## Installing packages:

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#if (!require("ChIPseeker")) BiocManager::install("ChIPseeker")
library("ChIPseeker")

# paketes de organismos concretos

#if (!require(txdb)) BiocManager::install(txdb)
library("TxDb.Athaliana.BioMart.plantsmart28")
txdb <- TxDb.Athaliana.BioMart.plantsmart28

#if (!require(annoDb_name)) BiocManager::install(annoDb_name)
library("org.At.tair.db")

#if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
library("clusterProfiler")

#BiocManager::install("ggnewscale")
library("ggnewscale")


## Getting genetic information from organism

genes.organism <- as.data.frame(genes(txdb))
head(genes.organism)

## La primera columna es el cromosoma donde se encuentra cada gen, PONER ENGLISH
## y en este caso nos interesa solo el del cromosoma 1, que sera
## nuevtro universo

my.universe<-genes.organism$gene_id

## Reading peaks

peak <- readPeakFile(peakfile = narrow_peak)
peak
#ACORDARSE QUE LOS NARROW IBAN EN COMILLA EN LOS ARGUMENTOS
## Reading summits, maximo de los picos

summit <- readPeakFile(peakfile = summits_bed)
summit 

## Peak localization

covplot(peak, weightCol="V5", title = title)

## Defining promoters

Promoters_peak <- getPromoters(TxDb=txdb, upstream=up1, downstream=down1)

Promoters_summit <- getPromoters(TxDb=txdb, upstream=up2, downstream=down2)

## Peak and summit annotation

peakAnno_peak <- annotatePeak(peak = peak, tssRegion = c(-up1,down1), TxDb = txdb, annoDb = "org.At.tair.db" )
plotAnnoPie(peakAnno_peak, main=title)

peakAnno_summit <- annotatePeak(peak = summit, tssRegion = c(-up2,down2), TxDb = txdb, annoDb = "org.At.tair.db")
plotAnnoPie(peakAnno_summit, main=title)

df.annotation_peak <- as.data.frame(peakAnno_peak)
head(df.annotation_peak)

df.annotation_summit <- as.data.frame(peakAnno_summit)
head(df.annotation_summit)


## En la columna de anotación del df.annotation, te indica con qué
## zona del gen se corresponde la localización de la cumbre. Solo
## interesan aquellas que estén en el promotor porque serán las que
## probablemente se corresponden con zonas de regulación de la 
## transcripción.

##LINEAS PARA PODER HACER EL EXAMEN:
#total.annotations <- unique(df.annotation_summit$annotation)
#excluir <- c("Distal Intergenic","Downstream (1-2kb)","Downstream (2-3kb)","Downstream (<1kb)")
#my.annotation <- setdiff(total.annotations,excluir)


#epi.df.annotation <- subset(df.annotation_summit, annotation %in% my.annotation)
#head(epi.df.annotation)


## FIN LINEA NUEVA

promoter.df.annotation_peak <- subset(df.annotation_peak, annotation == "Promoter")
head(promoter.df.annotation_peak)

## Solo queremos además el nombre del gen que se encuentra aguas abajo
## de esta región que suponemos que será su promotor, por eso hacemos
## lo siguiente.

#organism.regulome <- epi.df.annotation$geneId
#length(organism.regulome)
#write.table(organism.regulome, file = "regulome.txt", sep = "\t",
           # row.names = FALSE)


organism.regulome <- promoter.df.annotation_peak$geneId
length(organism.regulome)
write.table(organism.regulome, file = "regulome.txt", sep = "\t",
           row.names = FALSE)

## Ontology terms with clusterprofiler

## ont = BP de biological process, mirar tb los demas pal trabajo de verdad

ego <- enrichGO(gene          = organism.regulome,
                universe      = my.universe,
                OrgDb         = org.At.tair.db,
                ont           = "ALL",
                pvalueCutoff  = p_value_go,
                keyType = 'TAIR')

head(ego)

dotplot(ego, showCategory = 30, title = title)
barplot(ego, showCategory = 30, title = title)
cnetplot(ego, title = title)


## key type es la base de datos usada para el gene id

## que se supone q me tiene q salir al hacer head(ego?) a el le salen 30
## enriched terms found y a mi 256 xd 
## y los plots a partir de aki me he keao atras

#barplot(ego, showCategory=20)

## le sale que regula respuesta a carriquina pero a mi no. es un
## compuesto de incendio, ayuda a germinacion despues de incendio. 
## germinacion tras incendio es dependiente del reloj circadiano y 
## especificamente de prr5. 

#dotplot(ego)

## A el le gusta el siguiente, que es una red de concepto. Une genes
## a terminos de GO

#cnetplot(ego)

## Tb esta bien porque dice que prr5 esta involucrado en desarrollo de
## organos (a mi no me sale) y en respuestas a falta de agua, sust org,
## y luego de forma aislada a carriquina y frio 

#emapplot(ego)

## Otro tipo de enriquecimiento que no hemos visto es el de las rutas
## KEGG (mirar por nuestra cuenta). 

AT<- search_kegg_organism('Arabidopsis thaliana', by= 'scientific_name')
dim(AT)

pathway.enrich <- enrichKEGG(gene=organism.regulome, 
                             organism="ath",
                             keyType = "kegg", 
                             pvalueCutoff = p_value_kegg, 
                             pAdjustMethod = "BH")

head(pathway.enrich)

dotplot(ego,showCategory = 30 , title = title)
barplot(ego,showCategory = 30 , title = title)
cnetplot(ego , title = title)

