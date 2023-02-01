rm(list = ls())
load(file = "GSE3189_exprSet_plot_20181203.Rda")
library(ggplot2)

exprSet <- exprSet_sub

test <- exprSet[1:10,1:10]


index <- which(colnames(exprSet) %in% c("NTSR1","AATK","RARA-AS1","WNT9A","GAP43","AREG","MAST4","IGFBP5","COL5A3","FLT1","ISLR","PLAT","KRT16",
                                        "KCTD11","KCNJ12","STRA6","NKILA","FOSB","SLC25A45","IGFBP3","PFKFB4","CEACAM1","LTBP3","PSPH","ESM1",
                                        "CCNA1","C5orf45","SPRR2D","FOXC1","PINLYP","NRIP1","ITGA1","JAG1","NPTX2","STC1","EGR2","ADAMTS4","ACPP",
                                        "UBD","AP3B2","EGR3","CHAC1","KIAA0040","MN1","ADM2","FERMT1","FCRLA","INTS4P2","SESN2"))
#index <- which(colnames(exprSet) %in% c("FOXC1","EGR2","EGR3"))
#index <- which(colnames(exprSet) %in% c("CXCL3 ","PTGS2 ","IL1B","TNFRSF9","IL11","LCN2","ADAM19","AMACR","TRAF1","MMP9"))

expr10 <- exprSet[c(1:2,index)]
expr10 <- expr10[-1]
expr_gather <- tidyr::gather(expr10,key = Genenames, value = Geneexpr,-c(sample))


library(ggplot2)


# 箱线图,分基因查看癌和癌旁
ggplot(
  expr_gather, 
  aes(x = sample, y = Geneexpr,color=sample),
  title="GSE3189"
) +
  geom_boxplot()+
  facet_wrap(~Genenames,nrow = 1)



##############################
load(file = "GSE3189_exprSet_plot_20181203.Rda")
library(ggplot2)

exprSet <- exprSet_sub
library(tidyr)
library(dplyr)

test <- exprSet[1:10,1:10]
down <-  read.table("7_DHC_down.txt",header = F,stringsAsFactors = F,sep = "\t")
load(file = "diffLab.Rda")
GSE3189 <- diffLab
GSE3189[,7] <- rownames(GSE3189)
GSE3189_7DHC <- inner_join(GSE3189,down,by=c("V7"="V1"))


index <- which(colnames(exprSet) %in% c("CXCL3 ","PTGS2 ","IL1B","TNFRSF9","IL11","LCN2","ADAM19","AMACR","TRAF1","MMP9
                                        "))




expr10 <- exprSet[c(1:2,index)]
expr10 <- expr10[-1]
expr_gather <- tidyr::gather(expr10,key = Genenames, value = Geneexpr,-c(sample))


library(ggplot2)


# 箱线图,分基因查看癌和癌旁
ggplot(
  expr_gather, 
  aes(x = sample, y = Geneexpr,color=sample),
  title="GSE3189"
) +
  geom_boxplot()+
  facet_wrap(~Genenames,nrow = 2)
#########黑素瘤中差异基因中的转录因子，上调的

if(! require("ggstatsplot")) install.packages("ggstatsplot")
library(ggstatsplot)
library(cowplot)
####TFs在黑素瘤表达值，前三名
p1 <- ggbetweenstats(data = exprSet_GSE3189, 
                     x = sample, 
                     y = LEF1,
                     title = "melanoma_GSE3189")
p2 <- ggbetweenstats(data = exprSet_GSE7553, 
                     x = sample, 
                     y = LEF1,
                     title = "melanoma_GSE7553")
p3 <- ggbetweenstats(data = exprSet_GSE46517, 
                     x = sample, 
                     y = LEF1,
                     title = "melanoma_GSE46517")
p4 <- ggbetweenstats(data = exprSet_GSE3189, 
                     x = sample, 
                     y = CDK2,
                     title = "melanoma_GSE3189")
p5 <- ggbetweenstats(data = exprSet_GSE46517, 
                     x = sample, 
                     y = CDK2,
                     title = "melanoma_GSE46517")
p6 <- ggbetweenstats(data = exprSet_GSE3189, 
                     x = sample, 
                     y = MITF,
                     title = "melanoma_GSE3189")
p7 <- ggbetweenstats(data = exprSet_GSE46517, 
                     x = sample, 
                     y = MITF,
                     title = "melanoma_GSE46517")
plot_grid(p1,p2,p3,labels = LETTERS[1:3],nrow = 1)
plot_grid(p4,p5,labels = LETTERS[1:2],nrow = 1)
plot_grid(p6,p7,labels = LETTERS[1:2],nrow = 1)
####其余TFs在黑素瘤表达值

rm(list = ls())
load(file = "GSE3189_exprSet_plot_20181203.Rda")
library(ggplot2)

exprSet <- exprSet_sub

test <- exprSet[1:10,1:10]


index <- which(colnames(exprSet) %in% c("HOXD13","SOX5","CDK9","LEF1","CDK2","POU3F2","ELF1","MITF","LEF1","TBX2",
                                        "FOXD1","E2F6","POU3F4","HOXB7","EN2","SOX4","HIP1","DBF4","USF2",
                                        "POU3F3","BRCA1","MXI1","LHX1","YY1","EMX1","STAT1","PAX5","DAP","HOXA3",
                                        "RFX2","EVX1","HLTF","STAT3","GLI3","NEUROD1","FOXD2","HOXB7","MAZ",
                                        "SOX12","GBX2","PBX2","ATF6","E2F1","HOXB13","DDB1"))


expr10 <- exprSet[c(1:2,index)]
expr10 <- expr10[-1]
expr_gather <- tidyr::gather(expr10,key = Genenames, value = Geneexpr,-c(sample))


library(ggplot2)


# 箱线图,分基因查看癌和癌旁
ggplot(
  expr_gather, 
  aes(x = sample, y = Geneexpr,color=sample),
  title="GSE3189"
) +
  geom_boxplot()+
  facet_wrap(~Genenames,nrow =5)

