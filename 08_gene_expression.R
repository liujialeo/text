
getwd()
#看一下我们的工作目录下现在有的文件
dir()

#获取需要处理的数据,这个数据是1208个乳腺癌样本的转录组数据，我用Rdata的形式储存起来
#我们先把他load进来

load(file = "exprSet_GSE3189.Rdata")

test <- exprSet[1:5,1:5]

dim(exprSet )

test <- mRNA_exprSet[1:5,1:10]
exprSet1 <- t(exprSet)
test <- exprSet1[1:10,1:10]
dim(exprSet1)
metadata <- data.frame( rownames(exprSet1)) 
metadata[1:7,2] <- "Normal"
metadata[8:52,2] <- "Tumor"
#制作metadata，不要管这个单词，这一步就是区别肿瘤和正常组

#在右边看一下metadata，可以打开的，没问题，但是我们发现他没有像样的列名
names(metadata) <- c("sample_id","sample")
#class(metadata$sample) 是字符向量
#将sample转化成factor,这个是为了和后面的表达矩阵的匹配，当然我是到了后面才知道的

metadata$sample <- as.factor(metadata$sample) #as系列转换

#感受一下group_by 结合summarise的用法
library(dplyr)
metadata %>% group_by(sample) %>% summarise(n())
#发现有112例normal，1096例肿瘤组织，确实这个包很耐用
table(metadata$sample)




#加入肿瘤和对照的信息，在metadata里面
exprSet_sub <- cbind(metadata,exprSet1) #列合并
test <- exprSet_sub[1:20,1:20]

dim(exprSet_sub)
save(exprSet_sub,file = "GSE3189_exprSet_plot_20181203.Rda")
rm(list=ls())
load(file = "GSE3189_exprSet_plot_20181203.Rda")

cancer_names <- data.frame (exprSet_sub[,1])
cancer_names[,2] <- "melanoma_GSE3189"
colnames(cancer_names) <- c("sample_id","cancer")
exprSet_GSE3189 <- cbind(cancer_names,exprSet_sub)
TEST <- exprSet_GSE3189[1:10,1:10]
exprSet_GSE3189 <- exprSet_GSE3189[,-3]
save(exprSet_GSE3189,file = "GSE3189_exprSet_plot_cancernames_20181206.Rda")


##好啦，作图测试
##前面无论出现什么结果都不要紧，清空，加载，我们开始上路
rm(list = ls())
load(file = "GSE3189_exprSet_plot_20181203.Rda")
library(ggplot2)

exprSet <- exprSet_sub
#看BRCA1基因的是癌和癌旁的表达
#成功后，尝试一下 TP53，ERBB2，ESR1，修改y的值就可以
ggplot(exprSet,aes(x=sample,y=UBE2M))+
  geom_boxplot()


#试试你知道的gene，要注意这里没有非编码的gene

## 美图是个刚需，任何时候都需要，但是不用刻意去学
## 试试下面这个更加简洁的代码，你会有新的认识。
# 需要安装这个包
if(! require("ggstatsplot")) install.packages("ggstatsplot")

## 开始！！
library(ggstatsplot)
ggbetweenstats(data = exprSet, 
               x = sample, 
               y = UBE2M,
               title = "UBE2M in GSE3189")



## 那么问题来了，美图如何保存呢？，比例随便你来调，plot，export。

## 相关性分析如何作图？？
## 柱状图
library(ggstatsplot)
ggscatterstats(data = exprSet, 
               y = FOXA1, 
               x = ESR1,
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#CC79A7", 
               yfill = "#009E73", 
               marginal.type = "histogram",
               title = "Relationship between FOXA1 and ESR1")

## 密度
library(ggstatsplot)
ggscatterstats(data = exprSet, 
               y = FOXA1, 
               x = ESR1,
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#CC79A7", 
               yfill = "#009E73", 
               marginal.type = "density",
               title = "Relationship between FOXA1 and ESR1")
## 箱型图
library(ggstatsplot)
ggscatterstats(data = exprSet, 
               y = FOXA1, 
               x = ESR1,
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#CC79A7", 
               yfill = "#009E73", 
               marginal.type = "boxplot",
               title = "Relationship between FOXA1 and ESR1")

#这些简单的操作结合在一起，最终完成了一个相对复杂的事情，This is your victory!