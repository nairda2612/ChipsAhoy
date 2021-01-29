
## Installing packages


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
library("ChIPseeker")

# paketes de organismos concretos

BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")
library("TxDb.Athaliana.BioMart.plantsmart28")

BiocManager::install("org.At.tair.db")
library("org.At.tair.db")

BiocManager::install("clusterProfiler")
library("clusterProfiler")

## Getting information from organism 

txdb <- TxDb.Athaliana.BioMart.plantsmart28

## Getting genetic information from organism

genes.atha <- as.data.frame(genes(txdb))
head(genes.atha)


## La primera columna es el cromosoma donde se encuentra cada gen, PONER ENGLISH
## y en este caso nos interesa solo el del cromosoma 1, que sera
## nuevtro universo

genes.chr1<-subset(genes.atha, seqnames=="1")

my.universe<-genes.chr1$gene_id



## Reading peaks

peak <- readPeakFile(peakfile = "1_peaks.narrowPeak")
peak

## Reading cumbres, maximo de los picos

peak2 <- readPeakFile(peakfile = "1_summits.bed")
peak2 

## En covplot es la quinta columna, normalmente "V5" aunque mi ordenador no,
## por eso he puesto yo X211, pq asi se llama en mi caso

## Peak localization

covplot(peak, weightCol="X211")

## Defining promoters

Promotores_peak <- getPromoters(TxDb=txdb, upstream=1000, downstream=0)

Promotores_peak2 <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)

## Peak annotation

help("annotatePeak")

#asegurar para cual se hace
#peakAnno_peak <- annotatePeak(peak = peak, tssRegion = c(-1000,0), TxDb = txdb, annoDb = "org.At.tair.db" )
#plotAnnoPie(peakAnno)

peakAnno_peak2 <- annotatePeak(peak = peak2, tssRegion = c(-1000,1000), TxDb = txdb, annoDb = "org.At.tair.db" )
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
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                keyType = 'TAIR')
head(ego)

## key type es la base de datos usada para el gene id

## que se supone q me tiene q salir al hacer head(ego?) a el le salen 30
## enriched terms found y a mi 256 xd 
## y los plots a partir de aki me he keao atras

barplot(ego, showCategory=20)

## le sale que regula respuesta a carriquina pero a mi no. es un
## compuesto de incendio, ayuda a germinacion despues de incendio. 
## germinacion tras incendio es dependiente del reloj circadiano y 
## especificamente de prr5.

dotplot(ego)

## A el le gusta el siguiente, que es una red de concepto. Une genes
## a terminos de GO
cnetplot(ego)


## Tb esta bien porque dice que prr5 esta involucrado en desarrollo de
## organos (a mi no me sale) y en respuestas a falta de agua, sust org,
## y luego de forma aislada a carriquina y frio
emapplot(ego)

## Otro tipo de enriquecimiento que no hemos visto es el de las rutas
## KEGG (mirar por nuestra cuenta). 
