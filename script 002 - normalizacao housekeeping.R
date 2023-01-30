rm(list =ls())

library(NanoStringNorm)

#carregando os dados sem os housekeeping filtrados
genes_norm<-read.csv("data_filter_gnorm.csv", header = T, sep = ",")

#arquivo que indica os grupos
traists <-read.csv("traist.csv",sep = ";", header = T)

#separando somente o nome das amostras
sample.names <- names(genes_norm)[-c(1:3)]


trait.strain <- data.frame(
  row.names = sample.names,
  strain1 = traists
)

trait.strain$strain1.ID <- NULL

#separando o annotation dos dados de expressãoa
NanoString.mRNA.anno <- genes_norm[,c(1:3)];
NanoString.mRNA.data <- genes_norm[,-c(1:3)];

############ normalizacao por house keeping

data.norm.HK <- NanoStringNorm(
  x = NanoString.mRNA.data,
  anno = NanoString.mRNA.anno,
  CodeCount = "none",
  Background = "none",
  SampleContent = "housekeeping.geo.mean",
  OtherNorm = "none",
  round.values = FALSE,
  take.log = TRUE,
  traits = trait.strain
)

#gerando o QC da normalizacao
pdf('noma_housekeeping.pdf');
Plot.NanoStringNorm(
  x = data.norm.HK,
  label.best.guess = TRUE,
  plot.type = 'all'
);
dev.off()

# salvando os dados normalizados
df_hk <- as.data.frame(data.norm.HK[["normalized.data"]])

write_csv(df_hk, "normali_hk.csv")
