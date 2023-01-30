rm(list =ls())

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(TBSignatureProfiler)
  library(DESeq2)
  library(Rtsne)
  library(umap)
  library(ggplot2)
  library(limma)
  library(DT)
  library(ComplexHeatmap)
  library(tidyverse)
})

dif <-read.csv("normali_hk.csv", sep = ";", header = T)
dif <- dif[!(dif$Code.Class  %in% c( "Positive", "Negative", 
                                     "Housekeeping")),]

traists <-read.csv("traist.csv",sep = ";", header = T)

sample.names <- names(dif)[-c(1:3)]

row.names(dif)<-dif$Name
dif$Code.Class <- NULL
dif$Name <- NULL
dif$Accession <-NULL

trait.strain <- data.frame(
  row.names = sample.names,
  strain1 = traists
)
trait.strain$strain1.ID <- NULL

#Mudando a variavel para 0 ou 1 verificar se esta correto
desings <- mutate(trait.strain,
                  strain1.RCP = replace(strain1.CASO, strain1.CASO == "2", 0),
                  strain1.DR = replace(strain1.CONTROLE, strain1.CONTROLE == "2", 0))

coldata <- as.matrix(desings)

counts <- dif[,match(rownames(coldata),colnames(dif))]

#aqui colocar o numero de genes que serao analisados 
#cirando o modelo
nano_data_all<- SummarizedExperiment(
                        assays = list(counts = as.matrix(
                        counts[1:758,])),colData = coldata)

nano_data_all <- mkAssay(nano_data_all, log = TRUE, 
                         counts_to_CPM = TRUE)
#fitando o modelo
fit_status <- lmFit(assay(nano_data_all, fc = 2), desings)

#criando a matrix de contraste para comparar os grupos
contrast.matrix_status<- makeContrasts(strain1.CASO - strain1.CONTROLE, 
                                       levels = desings)
fit_status <- contrasts.fit(fit_status,contrast.matrix_status)
fit_status <- eBayes(fit_status)

#tabela com a analise diferencial - ordenada por p- valor 
limmaRes_status <- topTable(fit_status, adjust.method = "fdr",
                            n = Inf, sort.by = "P")
#genes significante e o rowname
adjpcutoff <- 0.05
limmaRes_status.sig <- subset(limmaRes_status, P.Value < adjpcutoff)
sigGenes <-rownames(limmaRes_status.sig)

#criando um viewer com os genes significantes
datatable(limmaRes_status.sig ,options=list(
          scrollX=T,pageLength=10),rownames = T)

#salvando o resultado
write.csv(limmaRes_status.sig, file = "diferencial_genes.csv")

