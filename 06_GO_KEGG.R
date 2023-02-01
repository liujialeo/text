#########
#########
# GO分析，KEGG分析
###########################
###########################
rm(list = ls())
if(! require("org.Mm.eg.db")) BiocInstaller::biocLite("org.Mm.eg.db")

library("org.Mm.eg.db")
load(file = "diffLab_6.Rda")
devtools::install_github("GuangchuangYu/clusterProfiler")
library(clusterProfiler)

#获得基因列表
gene <- rownames(diffLab)
#基因名称转换，返回的是数据框
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# 如果网络不好就运行下面这个命令，我已经帮大家储存好了
#load(file = "gene_cluster_20180802.Rda")
head(gene)


#**GO分析**
#细胞组分
ego_CC <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Mm.eg.db,
                   ont = "all",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)



save(ego_CC, file = "go_enrichment.Rda")
load(file = "go_enrichment.Rda")
abc <- as.data.frame(ego_CC)

#**作图
library(ggplot2)
p <- dotplot(ego_CC, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p

#**KEGG分析**
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism ="mouse",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)


save(kk,file = "kegg_enrichment.Rda")
load(file = "kegg_enrichment.Rda")
abc <- as.data.frame(kk)
dotplot(kk)
###将图导出为PPT格式
library(export)
graph2ppt(file="KEGG.pptx", width=7, height=5)


##################################################
#######把上调的差异基因拿出来进行富集分析，步骤同上
rm(list = ls())
library(tidyr)
library(dplyr)
load(file = "allDiff_6.Rda")
load(file = "diffLab_6.Rda")
library(clusterProfiler)
diffLab_up <- subset(allDiff,(logFC> 1 ) & adj.P.Val < 0.05)

gene <- rownames(diffLab_up)
#基因名称转换，返回的是数据框
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(gene)


#**GO分析**

ego_CC <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Mm.eg.db,
                   ont = "all",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

save(ego_CC, file = "go_enrichment_up.Rda")
load(file = "go_enrichment_up.Rda")
abc <- as.data.frame(ego_CC)

#**作图
p <- dotplot(ego_CC, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p


#**KEGG分析**
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism ="mouse",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)


save(kk,file = "kegg_enrichment_up.Rda")
load(file = "kegg_enrichment_up.Rda")
abc <- as.data.frame(kk)
dotplot(kk)


##################################################
#######把下调的差异基因拿出来进行富集分析，步骤同上
rm(list = ls())
library(tidyr)
library(dplyr)
load(file = "allDiff_6.Rda")
load(file = "diffLab_6.Rda")
library(clusterProfiler)
diffLab_down <- subset(allDiff,(logFC< (-1) ) & adj.P.Val < 0.05)

gene <- rownames(diffLab_down)
#基因名称转换，返回的是数据框
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
head(gene)


#**GO分析**

ego_CC <- enrichGO(gene = gene$ENTREZID,
                   OrgDb= org.Mm.eg.db,
                   ont = "all",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.01,
                   readable = TRUE)

save(ego_CC, file = "go_enrichment_down.Rda")
load(file = "go_enrichment_down.Rda")
abc <- as.data.frame(ego_CC)

#**作图
p <- dotplot(ego_CC, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free")
p


#**KEGG分析**
kk <- enrichKEGG(gene = gene$ENTREZID,
                 organism ="mouse",
                 pvalueCutoff = 0.01,
                 qvalueCutoff = 0.01,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)


save(kk,file = "kegg_enrichment_down.Rda")
load(file = "kegg_enrichment_down.Rda")
abc <- as.data.frame(kk)
dotplot(kk)
