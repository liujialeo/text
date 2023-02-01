################################################
### 作者：果子
### 更新时间：2018-10-08
### 联系方式: guoshipeng2008@126.com
### 微信：guo2sky
## 练习GEO数据的处理流程

rm(list = ls())

#更换工作目录，可以在Rstudio中完成
# Session,Set working directory,Choose working directory,选择GEO

getwd()
#安装bioconductor包
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
if(! require("limma")) BiocInstaller::biocLite("limma")
if(! require("GEOquery")) BiocInstaller::biocLite("GEOquery")


##本次处理的GEO数据编号：GSE42872
#载入GEOquery包
library(GEOquery)

#如果网络比较好就往下运行，本次看看就行，不运行
#eSet <- getGEO("GSE42872", GSEMatrix=T, AnnotGPL=FALSE, destdir = "C:/Users/guozi/Desktop/R_for_biotrainee/GEO_practice/")
#exprSet=exprs(eSet[[1]])
#pdata=pData(eSet[[1]])

#如果网络不是很通畅，手工获取上一步matrix的链接，下载到目的文件夹
#在浏览器中打开
#ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42872/matrix/GSE42872_series_matrix.txt.gz
#或者先登录
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
#在GEO accession 中输入GSE42872 也可以找到

#手动读入数据这一步可以得到pData，但是在本次演示中没有特别重要的作用，可以不用运行
#Data <- getGEO(filename="GSE42872_series_matrix.txt.gz")
#pData <- pData(phenoData(Data))


###########
###########
#正式的练习从这里开始
###########
#先解压GSE42872_series_matrix.txt.gz，注意解压到当前目录，再读入
exprSet <- read.table("GSE22307_series_matrix.txt",comment.char="!",stringsAsFactors=F,header=T)
#comment.char="!" 意思是！后面的内容不要读取，可以打开文件看一下?read.table

save(exprSet,file = "exprSet_readGSE.Rdata")
