rm(list =ls())


library(nanostringr)
library(NormqPCR)
library(ggplot2)
library(ggpubr)

#ler os RCC's path = nome da pasta que estão os arquivos
Nano <- read_rcc(path = "GENE" )

Nano_gene <- Nano[["raw"]]

sample.names <- names(Nano_gene)[-c(1:3)]

#separando os controles, endogenos e housekeepings
datacontrols <- subset(Nano_gene, Code.Class %in% c("Positive","Negative"))
dataEndogenous <- subset(Nano_gene, Code.Class %in% "Endogenous")
dataHK <- subset(Nano_gene, Code.Class %in% "Housekeeping")

# Aplicando o algoritimo genorm 
hk_expr <- as.matrix(t(dataHK[,-1:-3]))
hk <- dataHK$Name

genorm <- selectHKs(hk_expr, method = "geNorm",
                    Symbols = hk , minNrHK = 2, 
                    log = FALSE)

rank <- as.data.frame(genorm$ranking)
M <- as.data.frame(genorm$meanM)
var <- as.data.frame(genorm$variation)


#separando a minima variacao
min <- rownames(var)[which.min(apply(var,MARGIN=1,min))]
cut1 <- which(rownames(var) == min) 
var$selection <- "unselected"
var[cut1:nrow(var),]$selection <- "selected"

var_names <- factor(rownames(var), levels = unique(rownames(var)))
first.last <- as.vector(var_names)
first <- first.last[1]
n <- nrow(var)
last <- first.last[n]

#plotando os housekeeping
png("dot_genorm.png", height = 300, width = 400)
ggplot(var, aes(var_names, var$`genorm$variation`)) +
  geom_point(
    aes(color = var$selection), 
    position = position_dodge(0.3), size = 3
  )+
  labs(color = "")+
  scale_x_discrete(name = "",breaks=c(first,min,last))+
  scale_y_continuous(name = "Pairwise variation")+
  scale_color_manual(values = c("unselected" = "red", "selected" = "blue"))+
  theme_pubclean()
dev.off()

  # selecionando os housekeepings
cut <- as.numeric(gsub(".*/","", min))
rank$selection <- "unselected"
rank[1:cut,]$selection <- "selected"


# Gerando Planilha com os housekeeping
write.csv(rank, "housekeeping_selecionados.csv")

# Remover os housekeeping nao selecionados
genes_removed <- rank[rank$selection == "unselected",]
rownames(genes_removed) <- genes_removed$`genorm$ranking`
genes_removed <- rownames(genes_removed)

tmp_data <- Nano_gene
data_filter <- tmp_data[!(tmp_data$Name %in% genes_removed),]

##salvar os dados sem os housekeeping nao selecionados
write.csv(data_filter, "data_filter_gnorm.csv", col.names = T, row.names = F)
