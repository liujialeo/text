#############
#############
####GSEA分析来啦

#获得基因列表
rm(list = ls())
load(file = "allDiff_6.Rda")
eg <- rownames(allDiff)
#基因名称转换，返回的是数据框
  eg = bitr(eg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

eg <- dplyr::distinct(eg,SYMBOL,.keep_all=TRUE)
# geneList 三部曲
geneList = allDiff[,1]
names(geneList) = eg$ENTREZID
geneList = sort(geneList, decreasing = TRUE)
geneList[1:5]

# GSEA 采用GO数据集
require(DOSE)
library(clusterProfiler)
library(stringr)
library(ggplot2)
# 需要网络，大家会很拥挤
gg_BP <- gseGO(geneList     = geneList,
            OrgDb        = org.Mm.eg.db,
            ont          = "all",
            nPerm        = 1000,
            minGSSize    = 1,
            maxGSSize    = 500,
            pvalueCutoff =1,
            verbose      = FALSE)
save(gg_BP,file = "gg_BP.Rda")
load(file="gg_BP.Rda")
b1 <- as.data.frame(gg_BP)
dotplot(gg_BP,showCategory=12,split=".sign")+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))

gg_CC <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               ont          = "CC",
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
save(gg_CC,file = "gg_CC.Rda")
load(file="gg_CC.Rda")
b2 <- as.data.frame(gg_CC)
dotplot(gg_CC,showCategory=12,split=".sign")+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))

gg_MF <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               ont          = "MF",
               nPerm        = 1000,
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
save(gg_MF,file = "gg_MF.Rda")
load(file="gg_MF.Rda")
b3 <- as.data.frame(gg_MF)
dotplot(gg_MF,showCategory=12,split=".sign")+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))

## 解决文字太长的问题

dotplot(gg,showCategory=12,split=".sign",title="gseGO_BP")+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
##gse_KEGG富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
save(kk2,file = "kk2.Rda")
head(kk2)
b4 <- as.data.frame(kk2)
dotplot(kk2,showCategory=12,split=".sign")+facet_grid(~.sign)+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
##单个通路作图
#gesaplot
pdf(file = "gseplot.pdf")
gg_BP_df <- as.data.frame(gg_BP)
which("NIK/NF-kappaB signaling"==gg_BP_df$Description)
index <- 374
as.data.frame(gg_BP_df)$Description[index]
gseaplot(gg_BP, geneSetID = as.data.frame(gg_BP_df)$ID[index],title = "GSEA_NIK/NF-kappaB signaling")
dev.off()

gg_BP_df <- as.data.frame(gg_BP)
which("steroid biosynthetic process"==gg_BP_df$Description)
index <- 270
as.data.frame(gg_BP_df)$Description[270]
gseaplot(gg_BP, geneSetID = as.data.frame(gg_BP_df)$ID[index],title = "steroid biosynthetic process")

gg_BP_df <- as.data.frame(gg_BP)
which("cholesterol homeostasis"==gg_BP_df$Description)
index <- 354
as.data.frame(gg_BP_df)$Description[354]
gseaplot(gg_BP, geneSetID = as.data.frame(gg_BP_df)$ID[index],title = "cholesterol homeostasis")

gg_BP_df <- as.data.frame(gg_BP)
which("sterol homeostasis"==gg_BP_df$Description)
index <- 355
as.data.frame(gg_BP_df)$Description[355]
gseaplot(gg_BP, geneSetID = as.data.frame(gg_BP_df)$ID[index],title = "sterol homeostasis")

gg_BP_df <- as.data.frame(gg_BP)
which("steroid metabolic process"==gg_BP_df$Description)
index <- 57
as.data.frame(gg_BP_df)$Description[57]
gseaplot(gg_BP, geneSetID = as.data.frame(gg_BP_df)$ID[index],title = "steroid metabolic process")

kk2_df <- as.data.frame(kk2)
which("Steroid hormone biosynthesis"==kk2_df$Description)
index <- 7
as.data.frame(kk2_df)$Description[7]
gseaplot(kk2, geneSetID = as.data.frame(kk2_df)$ID[index],title = "Steroid hormone biosynthesis")
############
sbp <- gg_BP_df[200,11]
sbp1 <- strsplit(sbp,"/")
sbp2 <- as.data.frame(sbp1)
colnames(sbp2) <- "col"
library(dplyr)
library(tidyr)
steroid_syn <- inner_join(sbp2,eg,by=c("col"="ENTREZID"))
allDiff2 <- allDiff
allDiff2[,7] <-rownames( allDiff2)
steroid_syn_geneexp <- inner_join(steroid_syn,allDiff2,by=c("SYMBOL"="V7")) 
## 其他的KEGG和等数据集，可以自行练习,看clusterProfiler的说明

###提取NFKB基因
sbp <- gg_BP_df[374,11]
sbp1 <- strsplit(sbp,"/")
sbp2 <- as.data.frame(sbp1)
colnames(sbp2) <- "col"
library(dplyr)
library(tidyr)
NFKB <- inner_join(sbp2,eg,by=c("col"="ENTREZID"))
allDiff2 <- allDiff
allDiff2[,7] <-rownames( allDiff2)
NFKB_geneexp <- inner_join(NFKB,allDiff2,by=c("SYMBOL"="V7")) 
colnames(NFKB_geneexp)[1] <-"gene_id" 
write.table(NFKB_geneexp,file="NFKB_gene.xls",sep="\t",quote=F)
