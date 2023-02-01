###############
####heatmap热图
#用行名提取数据
rm(list = ls())
load(file = "exprSet.Rdata")
load(file = "diffLab_6.Rda")
if(! require("pheatmap")) install.packages("pheatmap")
#制作一个分组信息
heatdata <- exprSet[rownames(diffLab),]
class <- c(rep("con",5),rep("2d",6),rep("4d",6),rep("6d",6)) 
annotation_col <- data.frame(class)
rownames(annotation_col) <- colnames(heatdata)

#如果注释出界，可以通过调整格子比例和字体修正
pheatmap(heatdata, #热图的数据
         cluster_rows = TRUE,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = F,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("green", "black","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 10, cellheight = 0.6,# 格子比例
         fontsize = 10)


###将胆固醇合成酶和ABCA转运蛋白相关基因从黑素瘤差异基因中挑选出来
rm(list = ls())
load(file = "exprSet_GSE3189.Rdata")
load(file = "diffLab.Rda")
if(! require("pheatmap")) install.packages("pheatmap")
library(tidyr)
library(dplyr)
heatdata1 <- exprSet[rownames(diffLab),]
heatdata1[,53] <- rownames(heatdata1)
heatdata2 <- read.table("cholesterol_biosynthesis.txt")
heatdata3 <- read.table("ABCA_transporters.txt")
heatdata5 <- rbind(heatdata2,heatdata3)
heatdata6 <- inner_join(heatdata1,heatdata5,by=c("V53"="V1"))
heatdata7 <- diffLab
heatdata7[,7] <- rownames(heatdata7)
heatdata8 <- inner_join(heatdata6,heatdata7,by=c("V53"="V7"))

rownames(heatdata8) <- heatdata8[,53]

heatdata8 <- heatdata8[,-(55:59)]
heatdata8 <- heatdata8 %>%
  arrange(desc(abs(logFC)))
rownames(heatdata8) <- heatdata8[,53]
heatdata8 <- heatdata8[,-(53:54)]

class <- c(rep("con",7),rep("tumor",45)) 
annotation_col <- data.frame(class)
rownames(annotation_col) <- colnames(heatdata8)
pheatmap(heatdata8)
pheatmap(heatdata8, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 15, cellheight = 36,# 格子比例
         fontsize = 10)#,display_numbers=T,number_format="%.2f",number_color="red",,fontsize_number=6)

###cholesterol relative genes from prostate cancer 
rm(list = ls())
load(file = "exprSet_GSE3189.Rdata")
load(file = "diffLab.Rda")
if(! require("pheatmap")) install.packages("pheatmap")
library(tidyr)
library(dplyr)
heatdata1 <- exprSet[rownames(diffLab),]
heatdata1[,53] <- rownames(heatdata1)
heatdata2 <- read.table( "cholesterol_relative_genes_from_prostate_cancer.txt")
heatdata3 <- inner_join(heatdata1,heatdata2,by=c("V53"="V1"))
heatdata5 <- diffLab
heatdata5[,7] <- rownames(heatdata5)
heatdata6<- inner_join(heatdata3,heatdata5,by=c("V53"="V7"))

rownames(heatdata6) <- heatdata6[,53]
write.table(heatdata6,file="GSE3189_cholesterol_meatbolic_genes.xls",sep="\t",quote=F)

heatdata6 <- heatdata6[,-(55:59)]
heatdata6 <- heatdata6 %>%
  arrange(desc(abs(logFC)))
rownames(heatdata6) <- heatdata6[,53]
heatdata6 <- heatdata6[,-(53:54)]

class <- c(rep("con",7),rep("tumor",45)) 
annotation_col <- data.frame(class)
rownames(annotation_col) <- colnames(heatdata6)
pheatmap(heatdata6)
pheatmap(heatdata6, #热图的数据
         cluster_rows = T,#行聚类
         cluster_cols = T,#列聚类，可以看出样本之间的区分度
         annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         scale = "row", #以行来标准化，这个功能很不错
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 15, cellheight = 30,# 格子比例
         fontsize = 10)#,display_numbers=T,number_format="%.2f",number_color="red",,fontsize_number=6)
