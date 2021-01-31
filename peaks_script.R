## Author: The Winx Club
## Email: winx@alum.alfea.es
## Date: January 30

##This script works with arguments in order to submit it with R script. 
##If you want to run it in Rstudio, don't use arguments. ARREGLAR!!!

# Arguments:
arguments <- commandArgs(trailingOnly = TRUE)

narrow_peak <- arguments[1] #Peaks file
summits_bed <- arguments[2] #Summits file
# posible argumento para up y down stream COMENTAR
up1 <- as.numeric(arguments[3])
down1 <- as.numeric(arguments[4])
up2 <- as.numeric(arguments[5])
down2 <- as.numeric(arguments[6])
maintitle <-(arguments[7])

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

## Getting genetic information from organism

genes.atha <- as.data.frame(genes(txdb))
head(genes.atha)

## La primera columna es el cromosoma donde se encuentra cada gen, PONER ENGLISH
## y en este caso nos interesa solo el del cromosoma 1, que sera
## nuevtro universo

genes.chr1<-subset(genes.atha, seqnames=="1")

my.universe<-genes.chr1$gene_id



## Reading peaks

peak <- readPeakFile(peakfile = narrow_peak)
peak
#ACORDARSE QUE LOS NARROW IBAN EN COMILLA EN LOS ARGUMENTOS
## Reading summits, maximo de los picos

peak2 <- readPeakFile(peakfile = summits_bed)
peak2 
#ACORDARSE QUE LOS SUMMIT IBAN EN COMILLA EN LOS ARGUMENTOS

## En covplot es la quinta columna, normalmente "V5" aunque mi ordenador no,
## por eso he puesto yo X211, pq asi se llama en mi caso

## Peak localization

covplot(peak, weightCol="V5")

## Defining promoters

Promoters_peak <- getPromoters(TxDb=txdb, upstream=up1, downstream=down1)

Promoters_peak2 <- getPromoters(TxDb=txdb, upstream=up2, downstream=down2)

## Peak annotation

#help("annotatePeak")

#asegurar para cual se hace
#peakAnno_peak <- annotatePeak(peak = peak, tssRegion = c(-1000,0), TxDb = txdb, annoDb = "org.At.tair.db" )
#plotAnnoPie(peakAnno)

peakAnno_peak2 <- annotatePeak(peak = peak2, tssRegion = c(-up2,down2), TxDb = txdb, annoDb = "org.At.tair.db")
plotAnnoPie(peakAnno_peak2)
# guardar en el futuro

df.annotation <- as.data.frame(peakAnno_peak2)
head(df.annotation)

## En la columna de anotación del df.annotation, te indica con qué
## zona del gen se corresponde la localización de la cumbre. Solo
## interesan aquellas que estén en el promotor porque serán las que
## probablemente se corresponden con zonas de regulación de la 
## transcripción.

promoter.df.annotation <- subset(df.annotation, annotation == "Promoter")
head(promoter.df.annotation)

## Solo queremos además el nombre del gen que se encuentra aguas abajo
## de esta región que suponemos que será su promotor, por eso hacemos
## lo siguiente.

prr5.regulome <- promoter.df.annotation$geneId
length(prr5.regulome)

## Ontology terms with clusterprofiler


## ont = BP de biological process, mirar tb los demas pal trabajo de verdad

ego <- enrichGO(gene          = prr5.regulome,
                universe      = my.universe,
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                keyType = 'TAIR')
head(ego)

## key type es la base de datos usada para el gene id

## que se supone q me tiene q salir al hacer head(ego?) a el le salen 30
## enriched terms found y a mi 256 xd 
## y los plots a partir de aki me he keao atras
#Abro un canal para generar el archivo png



## Genero la imagen 

barplot(ego, showCategory=20)
##Cierro el canal



## le sale que regula respuesta a carriquina pero a mi no. es un
## compuesto de incendio, ayuda a germinacion despues de incendio. 
## germinacion tras incendio es dependiente del reloj circadiano y 
## especificamente de prr5.

## Genero la imagen 

dotplot(ego)
##Cierro el canal



## A el le gusta el siguiente, que es una red de concepto. Une genes
## a terminos de GO

## Genero la imagen 

cnetplot(ego)

## Tb esta bien porque dice que prr5 esta involucrado en desarrollo de
## organos (a mi no me sale) y en respuestas a falta de agua, sust org,
## y luego de forma aislada a carriquina y frio

## Genero la imagen 

emapplot(ego)
##Cierro el canal




## Otro tipo de enriquecimiento que no hemos visto es el de las rutas
## KEGG (mirar por nuestra cuenta). 
