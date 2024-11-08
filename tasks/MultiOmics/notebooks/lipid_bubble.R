library(NMF)
library(stringr)
# library(eoffice)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
# library(rJava)
# library(xlsx)
# library(RColorBrewer)
library(circlize)
library(umap)
library(Rtsne)
library(ggrepel)
library(ggplot2)
library(clusterProfiler)
library(KEGGREST)
library(ConsensusClusterPlus)
library(viridis)
library(survminer)
library(survival)
library(gridExtra)
library(factoextra)


##### Prepare Data ######

filepath.metab='~/workstation/MetabSubtype/tasks/Subtype/data/Using/lipid.csv'
filepath.sample='~/workstation/MetabSubtype/tasks/Subtype/data/Using/sample.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
dim(df.metab)
df.sample<-read.csv(filepath.sample, header= TRUE, check.names=F)
length(unique((df.metab[,1])))
df.metab[df.metab=='']=NA
dim(df.metab)

df.sample[df.sample=='']=NA
df.sample[df.sample=='Neg']=NA
# df.sample=df.sample[!is.na(df.sample['Sample Name']),c('Sample Name','主要分型','生存时间分组','MT2结构')]
df.sample=df.sample[!is.na(df.sample['Sample Name']),]
row.names(df.sample)=df.sample[['Sample Name']]
df.sample=df.sample[rownames(df.metab),]
df.sample[match('1520',rownames(df.sample)),'oss']=0
dim(df.metab)
metab_num=ncol(df.metab)
df.raw_metab=df.metab[,1:metab_num]
dim(df.sample)
dim(merge(df.metab,df.sample,by.x="row.names",by.y="Sample Name"))
# Odd Chain
dim(df.raw_metab)
df.raw_metab=drop_odd_chain_cols(df.raw_metab)
dim(df.raw_metab)
metab_num=dim(df.raw_metab)[2]


##### Build DF #####
class_label='kmeans_2_clusters'
df.use_sample=df.cluster_result
df.use=merge(df.raw_metab,df.cluster_result.correction[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

volcano.test.results=pvalue_by_wilcox(df.use,metab_num,class_label)
volcano.test.results$FDR=p.adjust(volcano.test.results$pvalue,method="BH")

type1='1'
type2='2'
volcano.test.results$FC=apply(df.use[,1:metab_num], 2, 
                              function(x) 
                                mean(as.numeric(x[which(df.use[class_label] == type2)]))/
                                mean(as.numeric(x[which(df.use[class_label] == type1)])))
head(volcano.test.results)
# Log2 (FC)
volcano.test.results$log2FC <- log2(volcano.test.results[,'FC'])

volcano.test.results

volcano.test.results$lipid_class=sapply(volcano.test.results$Metabolites,function(col)
{str_extract(col,'([A-Z]*[a-z]*)*')})

#### Draw #####
# pdf("lipid_bubble.pdf")
lipo.bubble=volcano.test.results
lipo.bubble %>%
  filter(lipid_class %in% c('Cer','Hex','LPC','LPE','SM','PC','PE','PG','PI','PS','TAG','DAG','TG')) %>%
  arrange(desc(`FDR`)) %>%
  ggplot(aes(x=lipid_class, y=log2FC, size=-log10(FDR), color=lipid_class)) +
  # geom_point(alpha=0.4) +
  geom_jitter(alpha=0.4,width = 0.4)+
  scale_size(range = c(0.1, 3.5), name="-Log(P)" 
             +
               scale_fill_viridis(discrete=TRUE, guide=TRUE, option="A")
               )+
  theme(plot.title = element_text( size = 6, vjust = 6, hjust = 0.5, face = "plain"), 
        axis.title.x = element_text(size = 6, vjust = -2, hjust = 0.55, face = "plain"), 
        axis.title.y = element_text(size = 10, vjust = 2, face = "plain"),
        axis.text.x = element_text(size = 10, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
        legend.title = element_text(size = 6, face = "plain", colour = "black"),
        legend.text = element_text(size = 10), 
        legend.key.size = unit(0.1, "cm"), 
        legend.spacing = unit(0.1, "cm"))+
  # theme_bw() +
  # theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 0.45)+
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())+
  geom_hline(yintercept = c(-log2(1.2),log2(1.2)),lty=4,col="black",lwd=0.4)
# dev.off()