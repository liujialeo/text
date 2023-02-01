###############
###############
##volcano火山图
###############
###############
##用ggplot2
rm(list = ls())
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggrepel")
library(ggplot2)
library(ggrepel)
library(dplyr)
load(file = "allDiff_6.Rda")
data <- allDiff
data$significant <- as.factor(data$P.Value<0.05 & abs(data$logFC) > 1)
data$gene <- rownames(data)

ggplot(data=data, aes(x=logFC, y =-log10(P.Value),color=significant)) +
  geom_point(alpha=0.8, size=1.2)+
  scale_color_manual(values =c("black","red"))+
  labs(title="Volcanoplot", x="log2 (fold change)",y="-log10 (q-value)")+
  theme(plot.title = element_text(hjust = 0.4))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8)+
  #theme(legend.position='none')
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),   
        axis.line = element_line(colour = "black")) +
  #geom_text(data=subset(data, abs(logFC) > 3), aes(label=gene),col="red",alpha = 1)
geom_text_repel(data=subset(data, abs(logFC) > 5), aes(label=gene),col="black",alpha = 0.8)
#ggsave("vocanol.pdf",,width = 7.09, height =5.6,dpi = 300)
