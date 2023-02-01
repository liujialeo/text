##############
##############
####差异分析
#############
#加载limma包，用于校正和比较差异
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
if(! require("limma")) BiocInstaller::biocLite("limma")


#####该数据集是DSS诱导C57小鼠结肠炎，分别在给药0、2、4、6天进行转录组测序
#####首先比较0和2天的转录组数据差异，提取这两天的转录组数据
rm(list = ls())
library(limma)
load(file = "exprSet.Rdata")
exprSet_2 <- exprSet[,1:11]

#differential差异分析
##构建分组矩阵,这一段有两种方法，但是初学的时候特别容易误解，这一步分完全就是独立的
class <- c(rep("con",5),rep("treat",6)) 
design <- model.matrix(~factor(class))
colnames(design) <- c("con","treat")
design

#线性模型拟合
fit <- lmFit(exprSet_2,design)  ###这一步需要改一下数据集名称
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
save(allDiff,file = "allDiff_2.Rda")  ###改一下名称
load(file = "allDiff_2.Rda")  ###改一下名称
#写入表格
write.table(allDiff,file="limmaTab_2.xls",sep="\t",quote=F) ##改名称

#找出差异两倍以上，pvalue小于0.05，1078个
diffLab <- subset(allDiff,(logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05)
save(diffLab,file = "diffLab_2.Rda") ##改名称
write.table(diffLab,file="diffExp_2.xls",sep="\t",quote=F)  ##改名称

#####首先比较0和4天的转录组数据差异，提取这两天的转录组数据
rm(list = ls())
library(limma)
load(file = "exprSet.Rdata")
exprSet_4 <- exprSet[,c(1:5,12:17)]

#differential差异分析
##构建分组矩阵,这一段有两种方法，但是初学的时候特别容易误解，这一步分完全就是独立的
class <- c(rep("con",5),rep("treat",6)) 
design <- model.matrix(~factor(class))
colnames(design) <- c("con","treat")
design

#线性模型拟合
fit <- lmFit(exprSet_4,design)  ###这一步需要改一下数据集名称
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
save(allDiff,file = "allDiff_4.Rda")  ###改一下名称
load(file = "allDiff_4.Rda")  ###改一下名称
#写入表格
write.table(allDiff,file="limmaTab_4.xls",sep="\t",quote=F) ##改名称

#找出差异两倍以上，pvalue小于0.05，1078个
diffLab <- subset(allDiff,(logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05)
save(diffLab,file = "diffLab_4.Rda") ##改名称
write.table(diffLab,file="diffExp_4.xls",sep="\t",quote=F)  ##改名称

#####首先比较0和6天的转录组数据差异，提取这两天的转录组数据
rm(list = ls())
library(limma)
load(file = "exprSet.Rdata")
exprSet_6 <- exprSet[,c(1:5,18:23)]

#differential差异分析
##构建分组矩阵,这一段有两种方法，但是初学的时候特别容易误解，这一步分完全就是独立的
class <- c(rep("con",5),rep("treat",6)) 
design <- model.matrix(~factor(class))
colnames(design) <- c("con","treat")
design

#线性模型拟合
fit <- lmFit(exprSet_6,design)  ###这一步需要改一下数据集名称
#贝叶斯检验
fit2 <- eBayes(fit)
#输出基因
allDiff=topTable(fit2,adjust='fdr',coef=2,number=Inf) 
save(allDiff,file = "allDiff_6.Rda")  ###改一下名称
load(file = "allDiff_6.Rda")  ###改一下名称
#写入表格
write.table(allDiff,file="limmaTab_6.xls",sep="\t",quote=F) ##改名称

#找出差异两倍以上，pvalue小于0.05，1078个
diffLab <- subset(allDiff,(logFC> 1 | logFC< (-1)) & adj.P.Val < 0.05)
save(diffLab,file = "diffLab_6.Rda") ##改名称
write.table(diffLab,file="diffExp_6.xls",sep="\t",quote=F)  ##改名称

#####以上结果显示随着DSS诱导时间的增加，差异基因数量也逐渐增加，因此，接下来我想分别计算诱导6天的
#####上调差异基因和下调差异基因