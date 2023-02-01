rm(list = ls())

load(file = "exprSet_readGSE.Rdata")
#探针基因名转换
##platformMap 中有常见的平台个R注释包的对于关系
platformMap <- data.table::fread("platformMap.txt")

index <- "GPL1261"
paste0(platformMap$bioc_package[grep(index,platformMap$gpl)],".db")


if(!require("mouse430a2.db")) BiocInstaller::biocLite("mouse430a2.db")

#获取探针
probe2symbol_df <- toTable(get("mouse430a2SYMBOL"))

##如果没有安装好"hugene10sttranscriptclusterSYMBOL"，说明你不是很听话
#没有安装就load吧，没有关系的，我们就是要练习
#load(file = "probe2symbol_df.Rdata")

#两个数据框的行名第一个使用相同的名称，便于以后合并
names(exprSet)[1] <- names(probe2symbol_df)[1]

#看一下symbol有没有重复，发现只有18837个，
length(unique(probe2symbol_df$symbol))
# 所以需要去重


#之后要合并数据，看看数据类型
class(exprSet$probe_id)
class(probe2symbol_df$probe_id)
#发现这两个类型不一样，实际上我是后来知道的,报错的时候发现问题解决问题
exprSet$probe_id <- as.character(exprSet$probe_id)

###探针转换以及去重，获得最终的表达矩阵
library(dplyr)
exprSet  <- exprSet %>% 
  inner_join(probe2symbol_df,by="probe_id") #合并探针的信息
exprSet1 <- exprSet[,c(25,2:24)]%>%
  distinct(symbol,.keep_all = T)



exprSet <- exprSet1

rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
range(exprSet)

exprSet <- log2(exprSet)
save(exprSet,file = "exprSet.Rdata")
write.table(exprSet,file="exprSet.xls",sep="\t",quote=F)


