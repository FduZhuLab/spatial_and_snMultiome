#setwd("E:/Rworkplace/01-17")
#remove.packages("rlang")
#BiocManager::install(c("lifecycle"))
#install.packages("BiocManager")
#library(BiocManager)
#.libPaths(c("",""))
#BiocManager::install("ggplot2",force=TRUE)
#BiocManager::install("rlang",force=TRUE)
#BiocManager::install(c("ggplot2","ggpubr","ggsignif"))
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(readr)

data<-read.csv("feature_ExN.csv")
colnames(data)[1]<-"brain_region"
data$brain_region<-factor(data$brain_region)
# my_comparisons <- list(c("Frontal", "Occipital"),
#                        c("Occipital", "Temporal"),
#                        c("Frontal", "Temporal"))


#i=12 #12 #18
i=18
colnames(data)[i]
files<-paste(paste("./again2/",colnames(data)[i],sep = ""),"csv",sep = ".")
ksfiles<-paste(paste("./again2/",colnames(data)[i],sep = ""),"ks.csv",sep = ".")
res.kruskal <- data %>% kruskal_test(RepoR ~ brain_region)%>%
  write_csv(ksfiles)
res.kruskal
stat.test <- data %>%
  dunn_test(RepoR~brain_region,p.adjust.method = "BH")%>%
  add_xy_position(x = "brain_region",
                  fun="max",
                  step.increase=0.2)%>%
  write_csv(files)

stat.test

feature<-colnames(data)[i]
p<-ggplot(data,aes(brain_region,data[,i]))+
  geom_violin(aes(fill=brain_region),
              cex=1.2,
              color="white")+           
  geom_boxplot(aes(fill=brain_region),
               width=0.1,
               lwd=0.3,
               cex=1.2,
               color="white")+
  geom_jitter(shape=21,
              size=0.5,
              fill="#c0c0c0",
              col="#c0c0c0",
              position = position_jitter(width = 0.2))+
  stat_pvalue_manual(stat.test, label = "p.adj.signif",tip.length = 0,
                     size=2)+
  #coord_cartesian(ylim = c(0,0.05))+
  coord_cartesian(ylim = c(-8,0))+
  labs(y="",x="brain area",title = feature)+
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE)
  )+
  theme_classic(base_size = 20)+        
  scale_fill_manual(values = c('#a6212c','#f2911b','#0b9c9a'))+
  theme(
        text = element_text(size = 6),
        axis.text = element_text(color = 'black'),
        plot.title = element_text(hjust = 0.5,size=8),
        plot.subtitle = element_text(size=6),
        axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75,size =6),
        axis.text.y = element_text(size =6),
        axis.title.x=element_text(size=6),
        axis.title.y=element_text(vjust=2, size=6),
        legend.position = 'none')      

p
name<-paste(paste("./again2/",feature,sep = ""),"pdf",sep = ".")
ggsave(p, file=name, width=2.5,height=4)
 
 