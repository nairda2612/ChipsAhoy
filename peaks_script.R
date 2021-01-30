## Author: The Winx Club
## Email: winx@alum.alfea.es
## Date: January 30

##This script works with arguments in order to submit it with R script. 
##If you want to run it in Rstudio, don't use arguments. ARREGLAR!!!

# Arguments:
arguments <- commandArgs(trailingOnly = TRUE)

txdb <- arguments[1] #Total genes of AT
narrow_peak <- arguments[2] #Peaks file
summit_peak <- arguments[3] #Summits file
# posible argumento para up y down stream COMENTAR
annoDb_name <- arguments[4] #Anotaci�n de AT
annoDb_param <- arguments[5]
up1 <- as.numeric(arguments[6])
down1 <- as.numeric(arguments[7])
up2 <- as.numeric(arguments[8])
down2 <- as.numeric(arguments[9])
maintitle <- as.numeric(arguments[10])


## Installing packages:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("ChIPseeker")) BiocManager::install("ChIPseeker")
library("ChIPseeker")

# paketes de organismos concretos

if (!require(txdb)) BiocManager::install(txdb)
library(txdb, character.only= TRUE)

if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
library("clusterProfiler")

if (!require(annoDb_name)) BiocManager::install(annoDb_name)
library(annoDb_name, character.only= TRUE)


## Getting information from organism 

####txdb <- TxDb.Athaliana.BioMart.plantsmart28

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

peak2 <- readPeakFile(peakfile = summit_peak)
peak2 
#ACORDARSE QUE LOS SUMMIT IBAN EN COMILLA EN LOS ARGUMENTOS

## En covplot es la quinta columna, normalmente "V5" aunque mi ordenador no,
## por eso he puesto yo X211, pq asi se llama en mi caso

## Peak localization

covplot(peak, weightCol="X211")

## Defining promoters

Promoters_peak <- getPromoters(TxDb=txdb, upstream=up1, downstream=down1)

Promoters_peak2 <- getPromoters(TxDb=txdb, upstream=up2, downstream=down2)

## Peak annotation

#help("annotatePeak")

#asegurar para cual se hace
#peakAnno_peak <- annotatePeak(peak = peak, tssRegion = c(-1000,0), TxDb = txdb, annoDb = "org.At.tair.db" )
#plotAnnoPie(peakAnno)

peakAnno_peak2 <- annotatePeak(peak = peak2, tssRegion = c(-up2,down2), TxDb = txdb, annoDb = annoDb_name)
plotAnnoPie(peakAnno)


df.annotation <- as.data.frame(peakAnno_peak)
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
                OrgDb         = annoDb_param,
                ont           = "BP",
                keyType = 'TAIR')
head(ego)

## key type es la base de datos usada para el gene id

## que se supone q me tiene q salir al hacer head(ego?) a el le salen 30
## enriched terms found y a mi 256 xd 
## y los plots a partir de aki me he keao atras
#Abro un canal para generar el archivo png


png(file=maintitle[1],
    width     = 10,
    height    = 10,
    units     = "in",
    res       = 600
)

## Genero la imagen 

barplot(ego, showCategory=20)
##Cierro el canal
dev.off()



## le sale que regula respuesta a carriquina pero a mi no. es un
## compuesto de incendio, ayuda a germinacion despues de incendio. 
## germinacion tras incendio es dependiente del reloj circadiano y 
## especificamente de prr5.

png(file=maintitle,
    width     = 10,
    height    = 10,
    units     = "in",
    res       = 600
)

## Genero la imagen 

dotplot(ego)
##Cierro el canal
dev.off()


## A el le gusta el siguiente, que es una red de concepto. Une genes
## a terminos de GO

png(file=maintitle,
    width     = 10,
    height    = 10,
    units     = "in",
    res       = 600
)

## Genero la imagen 

cnetplot(ego)
##Cierro el canal
dev.off()




## Tb esta bien porque dice que prr5 esta involucrado en desarrollo de
## organos (a mi no me sale) y en respuestas a falta de agua, sust org,
## y luego de forma aislada a carriquina y frio
png(file=maintitle,
    width     = 10,
    height    = 10,
    units     = "in",
    res       = 600
)

## Genero la imagen 

emapplot(ego)
##Cierro el canal
dev.off()



## Otro tipo de enriquecimiento que no hemos visto es el de las rutas
## KEGG (mirar por nuestra cuenta). 
