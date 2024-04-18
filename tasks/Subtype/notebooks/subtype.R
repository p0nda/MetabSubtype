library(NMF)
library(stringr)
# library(eoffice)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(rJava)
library(xlsx)
library(RColorBrewer)
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


getwd()
setwd("D:/repositories/liver-cancer-tasks/Tissue/notebooks/")
source('utils.R')

set.seed(6)

##### Metab #####
# Prepare Data
LABEL_NUM=0
# filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231202/csvs/lipid.csv'
# filepath.metab='D:/repositories/liver-cancer-tasks/Tissue/data/Using/metab.csv'
filepath.metab='D:/repositories/liver-cancer-tasks/Tissue/data/Using/lipid.csv'
filepath.sample='D:/repositories/liver-cancer-tasks/Tissue/data/Using/sample.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
df.sample<-read.csv(filepath.sample, header= TRUE, check.names=F)
df.metab[1:5,1:5]
dim(df.metab)
length(unique((df.metab[,1])))
unique(df.metab[,1])
df.metab[,1]
df.metab[df.metab=='']=NA
df.sample[df.sample=='']=NA
df.sample[df.sample=='Neg']=NA
# df.sample=df.sample[!is.na(df.sample['Sample Name']),c('Sample Name','主要分型','生存时间分组','MT2结构')]
df.sample=df.sample[!is.na(df.sample['Sample Name']),]
df.sample['Sample Name']
row.names(df.sample)=df.sample[['Sample Name']]
df.sample=df.sample[rownames(df.metab),]
df.sample[match('1520',rownames(df.sample)),'oss']=0
dim(df.metab)
metab_num=ncol(df.metab)
df.raw_metab=df.metab[,1:metab_num]
df.raw_metab.log=log2(df.raw_metab)
df.raw_metab[1]
# Normalize Data
df.raw_metab.scaled=df.raw_metab.log
raw_metab_mean=apply(df.raw_metab.log,1,mean)
raw_metab_mean
length(raw_metab_mean)
raw_metab_std=apply(df.raw_metab.log,1,sd)
df.raw_metab.scaled=(df.raw_metab.scaled-raw_metab_mean)/raw_metab_std
df.raw_metab.scaled=as.data.frame(nneg(as.matrix(df.raw_metab.scaled),method='min'))
colnames(df.raw_metab.scaled)
dim(df.raw_metab.scaled)
order(df.sample[['os']])
df.sample[order(df.sample[['os']]),c('os','oss')]
##### Cluster #####
df.raw_metab.scaled
# Consensus 
setwd("D:/repositories/liver-cancer-tasks/Tissue/results/20240406/")
# rcc = ConsensusClusterPlus(data.matrix(df.raw_metab.scaled), 
                           maxK = 8, 
                           reps = 1000, 
                           pItem = 0.8, 
                           pFeature = 1,
                           clusterAlg = "km",
                           plot="pdf")
# Sample Subset
df.use_sample=df.sample
df.use_sample=df.use_sample[,'batch',drop=FALSE]

class_label='batch'
df.use=merge(df.raw_metab.scaled,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

nmf_k=3
kmeans_k=3
df.cluster_result=perform_clustering(df.raw_metab.scaled, metab_num,nmf_k, kmeans_k)
df.cluster_result.lipid.correction=df.cluster_result

dim(df.cluster_result)
dim(df.use)

# Use Saved Matrix
filepath.cluster='D:/repositories/liver-cancer-tasks/Tissue/results/20240406/cluster_7e-2.csv'
df.cluster_result<-read.csv(filepath.cluster, header= TRUE, check.names=F,row.names=1)
df.cluster_result.metab.correction=df.cluster_result
.Random.seed[1:5]
df.sample
# Save Sample Data
df.sample_cluster=merge(df.sample,df.cluster_result.correction,by="row.names")
row.names(df.sample_cluster)=df.sample_cluster[[1]]
df.sample_cluster=df.sample_cluster[,2:ncol(df.sample_cluster)]
df.sample_cluster
# write.csv(df.sample_cluster,"D:/repositories/liver-cancer-tasks/Tissue/results/20240406/sample_cluster_7e-2.csv")
## Test
df=df.cluster_result
reference_labels=df[,1]
label_counts <- table(reference_labels, df[, 2])
table(df[,1], df[, 3])
table(df[,2], df[, 4])
table(df[,1], df[, 2])
label_counts 
##
# df.cluster_result.correction=cluster_label_correction(df.cluster_result,'kmeans_3_clusters',c('kmeans_3_clusters'))
# df.cluster_result.correction=cluster_label_correction(df.cluster_result.correction,'nmf_3_clusters',c('nmf_2_clusters'))
df.cluster_result.correction
table(df.cluster_result.correction[,4])
dim(df.cluster_result.correction)
dim(df.raw_metab.scaled)
# df.use.clustered[,metab_num:ncol(df.use.clustered)]
df.use=merge(df.raw_metab.scaled,df.cluster_result.correction,by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

dim(df.use)

##### Analysis #####
p_cutoff=5e-2
kmeans_pvalue_df=pvalue_by_wilcox(df.use,metab_num,'kmeans_clusters')
nmf_pvalue_df=pvalue_by_wilcox(df.use,metab_num,'nmf_clusters')
draw_single_heatmap(df.use,metab_num,'kmeans_3_clusters',p_cutoff,FALSE )
draw_single_heatmap(df.use,metab_num,'nmf_3_clusters',p_cutoff,FALSE )
draw_single_heatmap(df.use,metab_num,'kmeans_2_clusters',p_cutoff,FALSE )
draw_single_heatmap(df.use,metab_num,'nmf_2_clusters',p_cutoff,FALSE )

##### Visualization #####
# TSNE
# data.matrix(df.raw_metab.scaled)
colnames(df.cluster_result.metab.correction)=lapply(colnames(df.cluster_result.metab.correction), function(col){
   paste0("metab_",col)})
colnames(df.cluster_result.lipid.correction)=lapply(colnames(df.cluster_result.lipid.correction), function(col){
  paste0("lipid_",col)})

df.cluster_result.metab.correction
df.cluster_result.combine=merge(df.cluster_result.lipid.correction[,c(3,4)],df.cluster_result.metab.correction[,c(3,4)],by="row.names")
rownames(df.cluster_result.combine)=df.cluster_result.combine[,1]
df.cluster_result.combine=df.cluster_result.combine[,2:ncol(df.cluster_result.combine)]
df.cluster_result.lipid.correction
df.raw_metab.scaled[match(rownames(df.cluster_result.combine),rownames(df.raw_metab.scaled)),]
X=normalize_input(data.matrix(df.raw_metab.scaled[match(rownames(df.cluster_result.combine),rownames(df.raw_metab.scaled)),]))
tsne_results <- Rtsne(X, dims = 2,perplexity=20)
tsne_coords <- tsne_results$Y

df.draw=data.frame(tsne_coords[, 1], tsne_coords[, 2])
colnames(df.draw)=c('V1','V2')
rownames(df.draw)=rownames(df.cluster_result.combine)
df.draw=merge(df.draw,df.cluster_result.combine,by.x = "row.names",by.y = "row.names")
row.names(df.draw)=df.draw[[1]]
df.draw=df.draw[,2:ncol(df.draw)]
df.draw=merge(df.draw,df.sample[,c('os','oss')],by.x = "row.names",by.y = "row.names")
row.names(df.draw)=df.draw[[1]]
df.draw=df.draw[,2:ncol(df.draw)]

df.draw
color_palette <- brewer.pal(n = 9, name = "Set1")
par(mfrow=c(2,2))
class_labels=colnames(df.cluster_result.combine)
myplots <- vector('list', length(class_labels))
index=1
for (class_label in class_labels){
  tmp_p=ggplot((data = df.draw), aes(x = V1, y = V2, shape = factor(!!sym(class_label)),color=os )) +
    geom_point() +
    # scale_color_manual(values = color_palette(n = 3))+
    scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    labs(title = "t-SNE Visualization with Cluster Colors", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2",fill='cluster') +
    theme_bw()
  myplots[[index]]=local({
    index=index
    ggplot((data = df.draw), aes(x = V1, y = V2, shape = factor(!!sym(class_label)),color=os )) +
      geom_point() +
      scale_color_gradient2(low = "red", mid = "white", high = "blue")+
      # scale_color_manual(values = color_palette(n = 3))+
      # scale_color_brewer(palette="Set1")+
      labs(title = "t-SNE Visualization with Cluster Colors", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2",fill='cluster') +
      theme_bw()
  })
  index=index+1
}
myplots[[1]]+myplots[[2]]+myplots[[3]]+myplots[[4]]

# Without Colored OS
class_labels=colnames(df.cluster_result.combine)
myplots <- vector('list', length(class_labels))
index=1
for (class_label in class_labels){
  tmp_p=ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)))) +
    geom_point() +
    # scale_color_manual(values = color_palette(n = 3))+
    # scale_color_gradient2(low = "red", mid = "white", high = "blue")+
    labs(title = "t-SNE Visualization with Cluster Colors", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2",fill='cluster') +
    theme_bw()
  myplots[[index]]=local({
    index=index
    ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)))) +
      geom_point() +
      # scale_color_manual(values = color_palette(n = 3))+
      scale_color_brewer(palette="Set1")+
      labs(title = "t-SNE Visualization with Cluster Colors", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2",fill='cluster') +
      theme_bw()
  })
  index=index+1
}
myplots[[1]]+myplots[[2]]+myplots[[3]]+myplots[[4]]

# UMAP
X=normalize_input(data.matrix(df.raw_metab.scaled))
umap_results <- umap(X, n_components = 2)
umap_results
umap_coords <- umap_results$layout
umap_coords
class_label='kmeans_clusters'
df.draw=data.frame(umap_coords[, 1], umap_coords[, 2])
colnames(df.draw)=c('V1','V2')
rownames(df.draw)=rownames(df.cluster_result.correction)
df.draw=merge(df.draw,df.cluster_result.correction,by.x = "row.names",by.y = "row.names")
row.names(df.draw)=df.draw[[1]]
df.draw=df.draw[,2:ncol(df.draw)]

df.draw
color_palette <- brewer.pal(n = 3, name = "Set1")
par(mfrow=c(2,2))
class_labels=colnames(df.cluster_result.correction)
myplots <- vector('list', length(class_labels))
index=1
for (class_label in class_labels){
  tmp_p=ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)) )) +
    geom_point() +
    #scale_color_manual(values = color_palette)+
    scale_color_brewer(palette="Set1")+
    labs(title = "Umap Visualization with Cluster Colors", x = "Umap Dimension 1", y = "Umap Dimension 2",fill='cluster') +
    theme_bw()
  myplots[[index]]=local({
    index=index
    ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)) )) +
      geom_point() +
      #scale_color_manual(values = color_palette)+
      scale_color_brewer(palette="Set1")+
      labs(title = "Umap Visualization with Cluster Colors", x = "Umap Dimension 1", y = "Umap Dimension 2",fill='cluster') +
      theme_bw()
  })
  index=index+1
}
myplots[[1]]+myplots[[2]]+myplots[[3]]+myplots[[4]]

##### Survival #####
df.sample
df.use.os=merge(df.cluster_result.lipid.correction,df.sample[,c('os','oss')],by='row.names')
df.use.os=df.use.os[order(df.use.os$os),,drop=FALSE]
# df.use.os=df.use.os[1:(nrow(df.use.os)-1),]
dim(df.use.os)
surv_object=Surv(df.use.os$os,df.use.os$oss)
df.use.os
# K-MEANS
kmeans_fit=survfit(surv_object~kmeans_2_clusters,data=df.use.os)
kmeans_plot=ggsurvplot(kmeans_fit, data = df.use.os,pval = TRUE,risk.table = TRUE)
# NMF
table(df.use.os$nmf_2_clusters)
table(df.use.os$nmf_2_clusters)
nmf_fit=survfit(surv_object~nmf_2_clusters,data=df.use.os)
nmf_plot=ggsurvplot(nmf_fit, data = df.use.os,pval = TRUE,risk.table = TRUE)
# nmf_plot
# Combine
arrange_ggsurvplots(list(kmeans_plot,nmf_plot),ncol=2,nrow=1, data = df.use.os,pval = TRUE)
# write.csv(df.cluster_result,"D:/repositories/liver-cancer-tasks/Tissue/results/20240406/cluster_7e-2.csv")

##### PCA Cluster #####
# Check percentage
# princomp
coor_matrix=cor(df.raw_metab.scaled)
pca.metab.prin=princomp(coor_matrix)
summary(pca.metab.prin)
pca.metab.prin$loadings

# prcomp
pca.metab=prcomp(df.raw_metab.scaled,scale.=TRUE)

df.metab.pca=pca.metab$x
props=apply(df.metab.pca,2,var)/sum(apply(df.metab.pca,2,var))
pca_percentage=cumsum(props)
pca_percentage

using_num=26
df.use=df.metab.pca[,1:using_num]
df.use=as.data.frame(nneg(as.matrix(df.use),method='min'))
nmf_k=3
kmeans_k=3
df.cluster_result=perform_clustering(df.use, using_num,nmf_k, kmeans_k)
df.cluster_result.correction=cluster_label_correction(df.cluster_result,'kmeans_3_clusters',c('kmeans_3_clusters'))
df.cluster_result.correction=cluster_label_correction(df.cluster_result.correction,'nmf_3_clusters',c('nmf_2_clusters'))

# df.use.clustered[,metab_num:ncol(df.use.clustered)]
df.use=merge(df.use,df.cluster_result.correction,by.x = "row.names",by.y = "row.names")
df.use=merge(df.raw_metab.scaled,df.cluster_result.correction,by.x = "row.names",by.y = "row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

dim(df.use.os)

##### Analysis #####
p_cutoff=5e-2
# kmeans_pvalue_df=pvalue_by_wilcox(df.use.clustered,using_num,'kmeans_clusters')
# nmf_pvalue_df=pvalue_by_wilcox(df.use.clustered,using_num,'nmf_clusters')
kmeans_pvalue_df=pvalue_by_wilcox(df.use,metab_num,'kmeans_clusters')
nmf_pvalue_df=pvalue_by_wilcox(df.use,metab_num,'nmf_clusters')
draw_single_heatmap(df.use,using_num,'kmeans_3_clusters',p_cutoff,FALSE )
draw_single_heatmap(df.use,using_num,'nmf_3_clusters',p_cutoff,FALSE )
draw_single_heatmap(df.use,using_num,'kmeans_2_clusters',p_cutoff,FALSE )
draw_single_heatmap(df.use,using_num,'nmf_2_clusters',p_cutoff,FALSE )

# TSNE
# data.matrix(df.raw_metab.scaled)
X=normalize_input(data.matrix(df.raw_metab.scaled))
tsne_results <- Rtsne(X, dims = 2,perplexity=20)
tsne_coords <- tsne_results$Y

dim(df.cluster_result.correction)
df.draw=data.frame(tsne_coords[, 1], tsne_coords[, 2])
colnames(df.draw)=c('V1','V2')
rownames(df.draw)=rownames(df.cluster_result.correction)
df.draw=merge(df.draw,df.cluster_result.correction,by.x = "row.names",by.y = "row.names")
row.names(df.draw)=df.draw[[1]]
df.draw=df.draw[,2:ncol(df.draw)]

df.draw
color_palette <- brewer.pal(n = 9, name = "Set1")
par(mfrow=c(2,2))
class_labels=colnames(df.cluster_result)
myplots <- vector('list', length(class_labels))
index=1
for (class_label in class_labels){
  tmp_p=ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)) )) +
    geom_point() +
    # scale_color_manual(values = color_palette(n = 3))+
    labs(title = "t-SNE Visualization with Cluster Colors", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2",fill='cluster') +
    theme_bw()
  myplots[[index]]=local({
    index=index
    ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)) )) +
      geom_point() +
      # scale_color_manual(values = color_palette(n = 3))+
      scale_color_brewer(palette="Set1")+
      labs(title = "t-SNE Visualization with Cluster Colors", x = "t-SNE Dimension 1", y = "t-SNE Dimension 2",fill='cluster') +
      theme_bw()
  })
  index=index+1
}
myplots[[1]]+myplots[[2]]+myplots[[3]]+myplots[[4]]

# UMAP
X=normalize_input(data.matrix(df.raw_metab.scaled))
umap_results <- umap(X, n_components = 2)
umap_results
umap_coords <- umap_results$layout
umap_coords
class_label='kmeans_clusters'
df.draw=data.frame(umap_coords[, 1], umap_coords[, 2])
colnames(df.draw)=c('V1','V2')
rownames(df.draw)=rownames(df.cluster_result.correction)
df.draw=merge(df.draw,df.cluster_result.correction,by.x = "row.names",by.y = "row.names")
row.names(df.draw)=df.draw[[1]]
df.draw=df.draw[,2:ncol(df.draw)]

df.draw
color_palette <- brewer.pal(n = 3, name = "Set1")
par(mfrow=c(2,2))
class_labels=colnames(df.cluster_result)
myplots <- vector('list', length(class_labels))
index=1
for (class_label in class_labels){
  tmp_p=ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)) )) +
    geom_point() +
    #scale_color_manual(values = color_palette)+
    scale_color_brewer(palette="Set1")+
    labs(title = "Umap Visualization with Cluster Colors", x = "Umap Dimension 1", y = "Umap Dimension 2",fill='cluster') +
    theme_bw()
  myplots[[index]]=local({
    index=index
    ggplot((data = df.draw), aes(x = V1, y = V2, color = factor(!!sym(class_label)) )) +
      geom_point() +
      #scale_color_manual(values = color_palette)+
      scale_color_brewer(palette="Set1")+
      labs(title = "Umap Visualization with Cluster Colors", x = "Umap Dimension 1", y = "Umap Dimension 2",fill='cluster') +
      theme_bw()
  })
  index=index+1
}
myplots[[1]]+myplots[[2]]+myplots[[3]]+myplots[[4]]

##### Survival #####

surv_object=Surv(df.sample$os,df.sample$oss)
# K-MEANS
kmeans_fit=survfit(surv_object~kmeans_3_clusters,data=df.cluster_result)
kmeans_plot=ggsurvplot(kmeans_fit, data = df.cluster_result,pval = TRUE,risk.table = TRUE)
# NMF
nmf_fit=survfit(surv_object~nmf_3_clusters,data=df.cluster_result)
nmf_plot=ggsurvplot(nmf_fit, data = df.cluster_result,pval = TRUE,risk.table = TRUE)
# Combine
arrange_ggsurvplots(list(kmeans_plot,nmf_plot),ncol=2,nrow=1, data = df.use.os,pval = TRUE)

