
# install.packages("survminer") # 安装survminer包
# install.packages("survival")

library(stringr)
# library(eoffice)
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
packageVersion("rlang")    #查看指定R包版本

getwd()
setwd("D:/repositories/liver-cancer/tasks/Tissue/notebooks/")
source('utils.R')
setwd("D:/repositories/liver-cancer/tasks/Tissue/results/20231103/full_figure/")
# install.packages("rlang")

# detach("package:rlang", unload=TRUE)    

# BiocManager::install("ConsensusClusterPlus")

##### Cluster ######
filepath.metab="D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv"
df.metab=read.csv(filepath.metab,, header= TRUE, check.names=T,row.names=1)
df.metab=df.metab[,1:(ncol(df.metab)-4)]
df.metab=t(df.metab)
x <- df.metab
df.metab[104:106,1:5]


# UMAP降维  
umap <- umap(df.metab)
umap <- cbind(umap, k_means_3 = df.result$k_means_3)
# 绘制UMAP可视化
ggplot(umap, aes(V1, V2, color = k_means_3)) + 
  geom_point()

# t-SNE降维
tsne <- Rtsne(df)
tsne <- cbind(tsne$Y, k_means_3 = df$k_means_3)

# 绘制t-SNE可视化
ggplot(tsne$Y, aes(V1, V2, color = k_means_3)) +
  geom_point()

##### Consensus Clustering #####

LABEL_NUM=4
filepath.metab="D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv"
df.metab=read.csv(filepath.metab,, header= TRUE, check.names=T,row.names=1)
df.metab=t(df.metab)
df.metab[nrow(df.metab),]
dim(df.metab)
na_pos <- which(is.na(df.metab), arr.ind = T)
na_pos
# df.metab
colnames(df.metab)
rownames(df.metab[1:(nrow(df.metab)-LABEL_NUM),])
dim(df.metab)
dim(df.metab[1:(nrow(df.metab)-LABEL_NUM),])
rcc = ConsensusClusterPlus(df.metab[1:(nrow(df.metab)-LABEL_NUM),], 
                           maxK = 8, 
                           reps = 100, 
                           pItem = 0.8, 
                           pFeature = 1,
                           clusterAlg = "km",
                           plot="pdf")

# plot(res, ylab="Cluster consistency", xlab="Samples")
data2=data.frame(rcc[[2]]$consensusMatrix,row.names=colnames(df.metab))
colnames(data2)=colnames(df.metab)
data2
anno2=as.data.frame(rcc[[2]]$consensusClass)
colnames(anno2)=c("Cluster")
anno2$Cluster=as.factor(anno2$Cluster)
anno2=arrange(anno2,Cluster)

pheatmap(
    data2[,row.names(anno2)],
    color = colorRampPalette(c("white","blue"))(10),
    clustering_method = "average",
    cluster_rows = T,
    cluster_cols = T,
    annotation_col = anno2,
    annotation_colors = list(
        Cluster=c('1'="#1F78B4",'2'="#A6CEE3")
    ),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    border_color = "white"
)
pheatmap(
    rcc[[3]]$consensusMatrix,
    color = colorRampPalette(c("white","blue"))(50),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation"

)

pheatmap(
    rcc[[4]]$consensusMatrix,
    color = colorRampPalette(c("white","blue"))(50),
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",

)
# 检查聚类数Elbow点
plot(res, type="elbow")

##### Get Significant Metabolites #####
df.name_map=read.csv("D:/repositories/liver-cancer/tasks/Tissue/data/20230830_delete_outliners/name_map_1.csv", header= TRUE, check.names=F)
LABEL_NUM=4
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20230915/csvs/metab.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
typeof(df.metab)
dim(df.metab)
metab_num=dim(df.metab)[2]-LABEL_NUM
# df.metab=cbind(df.metab.raw_data,df.metab[,(dim(df.metab)[2]-LABEL_NUM):dim(df.metab)[2]])
class_label="k_cluster_3"
colnames(df.metab)
which(colnames(df.metab)==class_label)
df.metab[class_label]
df.metab.raw_data=df.metab[,1:metab_num]
df.metab=cbind(df.metab.raw_data,df.metab[class_label])
dim(df.metab)
type1=0
type2=1
type3=2
df.test=rbind(df.metab[which(df.metab$k_cluster_3 == type1),], 
                 df.metab[which(df.metab$k_cluster_3 == type2),],
                 df.metab[which(df.metab$k_cluster_3 == type3),])
# df.test[,-1]                
df.test.results=NA
df.test.results=data.frame(Metabolites = (colnames(df.test)[1:metab_num]))            
df.test.results=df.test.results[]
colnames(df.test)
as.numeric(df.test[,2])
df.test$k_cluster_3
typeof(df.test)
dim(df.test)
typeof(df.test)
typeof(df.test$k_cluster_3)
# kruskal.test(df.test$k_cluster_3~df.test[,2],data=df.test,EXACT=FALSE)
apply(df.test[,1:metab_num],2,function(x) print(typeof(x)))
df.test.results$kruskal=apply(df.test[, 1:metab_num], 2, function(x) unlist(kruskal.test(as.numeric(x), df.test$k_cluster_3,exact=TRUE)[3]))
df.test.results$kruskal_BH=p.adjust(df.test.results$kruskal,method="BH")
df.test.results$cpd_id=df.name_map$KEGG[match(df.test.results$Metabolites,df.name_map$Query)]
df.test.results=df.test.results[!(is.na(df.test.results$cpd_id)),]
dim(df.metab)
df.test.results
df.test.results$Metabolites[df.test.results$kruskal_BH<=5e-2]
df.test.results$is_use=df.test.results$kruskal_BH<=5e-2
df.test.results
# write.table(df.test.results,'D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/test_result_after.csv',sep=',',row.names = FALSE)
listDatabases()

##### Survival #####

# Prepare Data
LABEL_NUM=4
filepath.metab="D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv"
filepath.sample="D:/repositories/liver-cancer/tasks/Tissue/data/SCLC_tissue/第二批/sample_info.csv"
df.metab=read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
df.sample=read.csv(filepath.sample, header= TRUE, check.names=F,row.names=1)
df.metab.survival=cbind(df.metab[,1:(ncol(df.metab)-4)],df.metab['k_cluster_3'])
df.metab.survival=merge(df.metab.survival,df.sample[,c('Sample Name','os','oss')],by.x='row.names',by.y='Sample Name')
colnames(df.metab.survival)[170:ncol(df.metab.survival)]
df.metab.survival[1:4,170:ncol(df.metab.survival)]
# Survival
surv_object=Surv(df.metab.survival$os,df.metab.survival$oss)
fit=survfit(surv_object~k_cluster_3,data=df.metab.survival)
fit
summary(fit)
ggsurvplot(fit, data = df.metab.survival,pval = TRUE)


##### KEGG TEST #####

# 获取hsa通路
hsa_paths <- keggGet("pathway", organism="hsa")
db.hsa=keggList("hsa") ## list all human genes
typeof(db.hsa)
db.hsa[1:10]
keggInfo('hsa')
keggGet('hsa','c00022')
keggFind('hsa','c00022')
pathways <- keggGet("pathway")
## Build Map
db.hsa<-keggList("pathway",'hsa')
db.hsa[1:5]
df.hsa=data.frame(pathway_id=names(db.hsa),pathway_name=db.hsa)
df.hsa$pathway_id=str_replace(df.hsa$pathway_id,'hsa','map')
df.hsa[1:3,]
dim(df.hsa)
listDatabases()
keggList('')
db.c_p=keggLink('compound','pathway')
# db.c_p
df.c_p=data.frame(pathway_id=names(db.c_p),cpd_id=db.c_p)
df.c_p$pathway_id=unlist(df.c_p$pathway_id)
df.c_p$cpd_id=unlist(df.c_p$cpd_id)
df.c_p[1:5,]
df.c_p$pathway_id=lapply(df.c_p$pathway_id,function(x) {
  substr(x, nchar(x)-7, nchar(x))
})
df.c_p=df.c_p[df.c_p$pathway_id %in% df.hsa$pathway_id,]
dim(df.c_p)
df.c_p$pathway_name=df.hsa$pathway_name[match(df.c_p$pathway_id,df.hsa$pathway_id)]
df.c_p$cpd_id=lapply(df.c_p$cpd_id,function(x) {
  substr(x, nchar(x)-5, nchar(x))
})
# pathway_name_dict=c()
# for (i in unique(df.c_p$pathway_id)){
#     tmp_pathway=keggGet(i)
#     tmp_name=tmp_pathway[[1]]$NAME[[1]][1]
#     pathway_name_dict[i]=tmp_name
# }
# typeof(unlist(df.c_p['pathway_id']))
# df.c_p[['pathway_id']][[1]]
# pathway_name_dict[df.c_p[['pathway_id']]]
# df.c_p$pathway_name=pathway_name_dict[unlist(df.c_p$pathway_id)]
# df.c_p[c(1,8,50,89,200),]
# df.c_p$pathway_id=unlist(df.c_p$pathway_id)
df.c_p$cpd_id=unlist(df.c_p$cpd_id)
class(df.c_p$pathway_id)
class(df.c_p$pathway_name)
class(df.c_p$cpd_id)
df.c_p[1:5,]
# write.table(df.c_p,'D:/repositories/liver-cancer/tasks/Tissue/data/kegg_hsa_pathway_map.csv',sep=',',row.names = FALSE)


pathways 
tmp_id="hsa05130"
tmp_pathway=keggGet(tmp_id)
tmp_pathway
tmp_vec=c()
tmp_pathway[[1]]$NAME[[1]]
tmp_vec[tmp_id]=tmp_pathway[[1]]$NAME[1]
tmp_vec
# Connection between map and hsa
hsa_paths <- keggGet("pathway", organism="hsa")
map_paths <- keggList("pathway")
hsa_paths<-keggList("pathway",'hsa')
length(hsa_paths)
map_paths
matches <- map_paths$ENTRY[grepl("^hsa", map_paths$ENTRY)]
map_paths
db.c_p[1:5]
names(db.c_p[1:5])
db.c_p[str_detect(names(db.c_p),'map')]
# 获取 Compounds 数据
compounds <- keggGet("compound")
# 合并 Pathways 和 Compounds 表
combined_df <- merge(pathways, compounds, by="cpd_id")
# 过滤出 hsa 开头的代谢通路
filtered_df <- combined_df[grepl("hsa", combined_df$Pathway_ID),]
# 选择需要的字段
results_df <- filtered_df[, c("Pathway_ID", "Pathway_Name")]

dim(df.name_map)
dim(df.metab)
df.metab
df.name_map
df.name_map
df.new
comp_id=df.name_map$KEGG[12]
comp_id
path_info <- keggGet("C00345")
path_info <- keggGet(comp_id, "pathway")
path_info <- keggLink("pathway", "C14156")
path_info[1]
colnames(merge(x=df.new,y=df.name_map,by.x="cpd_id",by.y="KEGG"))
unique(merge(x=df.new,y=df.name_map,by.x="cpd_id",by.y="KEGG")[1:200,c(1,4,5)])
colnames(df.new)

data(geneList, package='DOSE')
de <- names(geneList)[1:100]
de
yy <- enrichKEGG(de, pvalueCutoff=0.01)
head(yy)

# KEGG for more than 2 classes
# df.test.results
colnames(df.test)
as.numeric(df.test[,2])
df.test$k_cluster_3
typeof(df.test)
dim(df.test)
typeof(df.test)
typeof(df.test$k_cluster_3)
# kruskal.test(df.test$k_cluster_3~df.test[,2],data=df.test,EXACT=FALSE)
apply(df.test[,1:metab_num],2,function(x) print(typeof(x)))
df.test.results$kruskal=apply(df.test[, 1:metab_num], 2, function(x) unlist(kruskal.test(as.numeric(x), df.test$k_cluster_3,exact=TRUE)[3]))
df.test.results$kruskal_BH=p.adjust(df.test.results$kruskal,method="BH")
dim(df.metab)
df.test.results
df.test.results$cpd_id=df.new$cpd_id[match(df.test.results$Metabolites,df.new$query)]
df.test.results$Metabolites[df.test.results$kruskal_BH<=5e-2]
df.test.results$is_use=df.test.results$kruskal_BH<=5e-2
df.test.results
# KEGG enrichment analysis
path1 <- df.path.metab %>% dplyr::select(pathway_id, cpd_id) 
path2 <- df.path.metab %>% dplyr::select(pathway_id, pathway_name)
path1[1:5,]
# df.new=df.new[!is.na(df.new$Metabolism_Name),]
# write.table(df.test.results,'D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/test_result.csv',sep=',',row.names = FALSE)
length(df.test.results$cpd_id[df.test.results$is_use])
df.test.results$cpd_id[df.test.results$is_use]
path2=path2[!is.na(path2$pathway_name),]
path1$cpd_id=toupper(path1$cpd_id)
df.test.results$cpd_id=toupper(df.test.results$cpd_id)
used_cpds=df.test.results$cpd_id[df.test.results$is_use]
used_cpds
df.test.results[df.test.results$is_use,]
length(used_cpds)
# write.table(df.test.results,'D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/test_result.csv',sep=',',row.names = FALSE)
x <- enricher(used_cpds,
              TERM2GENE = path1, 
              TERM2NAME = path2,
              minGSSize = 2, pvalueCutoff = 10, pAdjustMethod = "BH")
x
kegg_table <- na.omit(as.data.frame(x))
kegg_table$Description=str_replace(kegg_table$Description,' - Homo sapiens \\(human\\)','')
# write.table(kegg_table,'D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/kegg_table.csv',sep=',',row.names = FALSE)

kegg_table$Description
kegg_table <- kegg_table[kegg_table$Count > 2,]

kegg_table$metabs_mean <- apply(kegg_table, 1, path_avelog2FC)

# calculate fold enrichment
kegg_table$FoldEnrich <- apply(kegg_table, 1, 
                               function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                 as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                 as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
colnames(kegg_table)


kegg_table$geneID
kegg_table[1:5,]
-log10(kegg_table[['p.adjust']])
-log10(kegg_table[['pvalue']])
kegg_table=kegg_table[order(kegg_table$FoldEnrich,decreasing = FALSE),]
kegg_table=kegg_table[order(kegg_table$p.adjust,decreasing = TRUE),]
kegg_table
p2 = ggplot(kegg_table, aes(-log10(pvalue), Description)) +
  geom_point(aes(fill = FoldEnrich, size = -log10(p.adjust)), color = "black", shape = 21) +
  scale_y_discrete(limits = kegg_table[['Description']])+
  scale_size(range = c(1, 3.2), breaks = c(2,5,8)) +
# scale_color_gradient(low="blue",high = "red")+
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.5, 1.5, 2.5)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.7) +
  labs(x = "", y = "", title = "", 
       fill = "abs. Log2 (Fold change)", size = "-Log10 (P value)") +
  theme(axis.text.x = element_text(size = 8, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 8, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 8, face = "plain", colour = "black"), 
        legend.title = element_text(size = 8, face = "plain", colour = "black"))
p2


##### KEGG #####

df.name_map=read.csv("D:/repositories/liver-cancer/tasks/Tissue/data/20230830_delete_outliners/name_map_1.csv", header= TRUE, check.names=F)
df.path.metab <- read.csv('D:/repositories/liver-cancer/tasks/Tissue/data/kegg_hsa_pathway_map.csv', header = T)
# df.path.metab$cpd_id=toupper(df.path.metab$cpd_id)
df.path.metab[1:3,]
dim(df.path.metab)
dim(df.name_map)
typeof(df.name_map)
# unique(df.path.metab$cpd_id)
# unique(df.name_map$KEGG)
df.new=df.path.metab[df.path.metab$cpd_id %in% df.name_map$KEGG,]
colnames(df.new)
dim(df.new)
length(unique(df.new$cpd_id))
df.new$query=NA
dim(df.new)
# df.new$query[df.new$cpd_id %in% df.name_map$KEGG]=df.name_map$Query
df.new$query=df.name_map$Query[match(df.new$cpd_id, df.name_map$KEGG)]
df.new[1:5,]
df.new=df.new[!is.na(df.new$cpd_id),]
# df.new$analyst_name=df.path.metab[df.path.metab$cpd_id %in% df.name_map$KEGG,]
# df.name_map[df.name_map$KEGG %in% df.path.metab$cpd_id,]
# df.new
length(unique(df.new$cpd_id))
length(unique(df.name_map$KEGG))
# unmatch_metabolite_ids=(setdiff(unique(df.name_map$KEGG),intersect(unique(df.name_map$KEGG),unique(df.new$cpd_id))))
# print(df.name_map[df.name_map$KEGG %in% unmatch_metabolite_ids,]$Match)
# df.new
# colnames(df.path.metab)[colnames(df.path.metab) == "cpd_id"] ="KEGG"
length(unique(df.new$cpd_id))
# df.new$query <- unlist(lapply(df.new$cpd_id, function(x) {
#   df.name_map$Query[which(df.name_map$KEGG==x)]
# }))
# df.new$query

# read metab and select
LABEL_NUM=4
filepath.metab="D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv"
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
typeof(df.metab)
dim(df.metab)
# Temp
# df.metab.raw_data=df.metab[,unique(df.new$query)]
df.metab.raw_data=df.metab[,1:(ncol(df.metab)-LABEL_NUM)]
metab_num=dim(df.metab.raw_data)[2]
metab_num
# df.metab=cbind(df.metab.raw_data,df.metab[,(dim(df.metab)[2]-LABEL_NUM):dim(df.metab)[2]])
class_label="k_cluster_3"
colnames(df.metab)
# which(colnames(df.metab)==class_label)
df.metab[class_label]
df.metab=cbind(df.metab.raw_data,df.metab[class_label])
dim(df.metab)
type1=0
type2=1
type3=2
df.test=rbind(df.metab[which(df.metab$k_cluster_3 == type1),], 
                 df.metab[which(df.metab$k_cluster_3 == type2),],
                 df.metab[which(df.metab$k_cluster_3 == type3),])

# Build df.test of 3 Clusters
df.test.cluster01=rbind(df.metab[which(df.metab$k_cluster_3 == type1),], 
                 df.metab[which(df.metab$k_cluster_3 == type2),])
df.test.cluster02=rbind(df.metab[which(df.metab$k_cluster_3 == type1),], 
                 df.metab[which(df.metab$k_cluster_3 == type3),])
df.test.cluster12=rbind(df.metab[which(df.metab$k_cluster_3 == type2),], 
                 df.metab[which(df.metab$k_cluster_3 == type3),])
df.test

# # Temp calculate pvalue
# tmp_p_value<-function(df.test,metab_num,types=c()){
#     df.test.results=data.frame(Metabolites = (colnames(df.test)[1:metab_num]))
#     df.test.results$wilcox=apply(df.test[, 1:metab_num], 2,
#                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.test[['k_cluster_3']], data = df.test, exact = FALSE)[3]))
#     df.test.results$wilcox_BH=p.adjust(df.test.results$wilcox,method="BH")
#     df.test.results$FC <- apply(df.test[,1:metab_num], 2, 
#                                 function(x) 
#                                 mean(as.numeric(x[which(df.test$k_cluster_3 == types[1])]))/
#                                 mean(as.numeric(x[which(df.test$k_cluster_3 == types[2])])))
#     df.test.results$log2FC<- log2(df.test.results$FC)
#     return(df.test.results)
# }

# df.test.results.cluster01=tmp_p_value(df.test.cluster01,metab_num,c(0,1))
# df.test.results.cluster01[1:5,1:5]
# write.table(df.test.results.cluster01,'D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/full_cpd_01.csv',sep=',',row.names = FALSE)
types=c(1,2)
apply(df.test[,1:metab_num], 2, 
                                function(x) 
                                mean(as.numeric(x[which(df.test$k_cluster_3 == types[1])]))/
                                mean(as.numeric(x[which(df.test$k_cluster_3 == types[2])])))
df.test.results.cluster01=kegg_before_process(df.test.cluster01,metab_num,c(0,1))
df.test.results.cluster02=kegg_before_process(df.test.cluster02,metab_num,c(0,2))
df.test.results.cluster12=kegg_before_process(df.test.cluster12,metab_num,c(1,2))
df.test.results.cluster01
df.test.results.cluster02
df.test.results.cluster12

# KEGG
setwd("D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/")

# write.csv(df.test.results.cluster01,'test_results_cluster01.csv')
# write.csv(df.test.results.cluster02,'test_results_cluster02.csv')
# write.csv(df.test.results.cluster12,'test_results_cluster12.csv')
kegg_table.cluster01=kegg_after_process(df.test.results.cluster01)
kegg_table.cluster01$comporison='C_0/C_1'
kegg_table.cluster01
kegg_table.cluster02=kegg_after_process(df.test.results.cluster02)
kegg_table.cluster02$comporison='C_0/C_2'
kegg_table.cluster02
kegg_table.cluster12=kegg_after_process(df.test.results.cluster12,0.05,0,0)
kegg_table.cluster12$comporison='C_1/C_2'
kegg_table.cluster12
queries.cluster12=df.new[match(unlist(strsplit(kegg_table.cluster12[['geneID']],"[/]")),df.new$cpd_id),'query']
queries.cluster12
#### Test
pvalue_cutoff=5e-2
fc_up_cutoff=0
fc_down_cutoff=0
# df.test.results.cluster12$is_use=(df.test.results.cluster12$wilcox_BH<=pvalue_cutoff)&((df.test.results.cluster12$FC>=fc_up_cutoff)|(df.test.results.cluster12$FC<=fc_down_cutoff))
df.test.results.cluster12$is_use=(df.test.results.cluster12$wilcox_BH<=pvalue_cutoff)
used_cpds=df.test.results.cluster12$cpd_id[df.test.results.cluster12$is_use]
used_cpds
print(paste0(used_cpds,collapse = '\n'))
# write.csv(df.test.results.cluster12,'test_results_cluster12.csv')
ggplot(kegg_table.cluster12, aes(-log10(p.adjust), Description)) +
  geom_point(aes(fill = metabs_mean, size = -log10(p.adjust)), color = "black", shape = 21) +
#   scale_y_discrete(limits = kegg_all[['Description']])+
  scale_size(range = c(3, 12), breaks = c(3,5,7,8)) +
# scale_color_gradient(low="blue",high = "red")+
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.2,0.5, 1, 1.5)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.4) +
  labs(x = "", y = "", title = "", 
       fill = "Abs.Log2 (FC)", size = "-Log10 (P-value)") +
  theme(axis.text.x = element_text(size = 20, face = "plain", colour = "black",angle = 45,hjust = 1), 
        axis.text.y = element_text(size = 18, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))
####

# Concat 3 kegg_table and Draw
kegg_all <- rbind(kegg_table.cluster01, kegg_table.cluster02, kegg_table.cluster12)
# kegg_all <- rbind(kegg_table.cluster01, kegg_table.cluster02)
kegg_all
# write.table(kegg_all,'D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/kegg_table.csv',sep=',',row.names = FALSE)

# kegg_all=read.csv('D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/kegg_table.csv',check.names=FALSE)
ggplot(kegg_all, aes(comporison, Description)) +
  geom_point(aes(fill = metabs_mean, size = -log10(p.adjust)), color = "black", shape = 21) +
#   scale_y_discrete(limits = kegg_all[['Description']])+
  scale_size(range = c(3, 12), breaks = c(3,5,7,8)) +
# scale_color_gradient(low="blue",high = "red")+
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.5,0.6,0.7,0.8, 0.9, 1.5)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.4) +
  labs(x = "", y = "", title = "", 
       fill = "Abs.Log2 (FC)", size = "-Log10 (P-value)") +
  theme(axis.text.x = element_text(size = 20, face = "plain", colour = "black",angle = 45,hjust = 1), 
        axis.text.y = element_text(size = 18, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))
p2 

tmp_kegg_table=kegg_all
ggplot(tmp_kegg_table, aes(comporison, Description)) +
  geom_point(aes(fill = metabs_mean, size = -log10(p.adjust)), color = "black", shape = 21) +
  scale_size(range = c(3, 12), breaks = c(2,7,12)) +
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.5, 0.8, 1.1)) +
  #theme_dendro() + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  coord_fixed(ratio = 0.4) +
  labs(x = "", y = "", title = "", 
       fill = "abs. Log2 (Fold change)", size = "-Log10 (P value)") +
  theme(axis.text.x = element_text(size = 6, face = "plain", colour = "black"), 
        axis.text.y = element_text(size = 6, face = "plain", colour = "black"),
  ) +
  theme(legend.text = element_text(size = 6, face = "plain", colour = "black"), 
        legend.title = element_text(size = 6, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))

tmp_kegg_table$query=apply(tmp_kegg_table,1,get_cpd_names)
# write.table(tmp_kegg_table,'D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/kegg_table.csv',sep=',',row.names = FALSE)

##### HeatMap ######  

LABEL_NUM=4
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv'
# filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/data/metab_raw.csv'
# filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab_transcript.csv'
# filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab_codex.csv'
# filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab_type.csv'
# filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/data/SCLC_tissue/Lipid/combat_corrected.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
dim(df.metab)
df=df.metab
row.names(df)
print(rownames(df))
# df <- df[!(row.names(df) %in% c("8195","2236","210012-T")),]
# df <- df[df[['k_cluster_3']] %in% c(1,2),]
# df <- df[df[['batch']]==1,]
#df.metab['codex_try']
df.metab=df
colnames(df.metab)
rownames(df.metab)
(ncol(df.metab)-1)
class_label='k_cluster_3'
!is.na(df.metab[[class_label]])
df.metab[[class_label]]
df.metab=df.metab[!df.metab[[class_label]]=='',]
# class_label='batch'
df.metab[class_label]
metab_num=ncol(df.metab)-LABEL_NUM
pvalue_cutoff=1000000
use_row_ha=FALSE
col_split=TRUE
draw_single_heatmap(df.metab,metab_num,class_label,pvalue_cutoff,use_row_ha )

#### Test
LABEL_NUM=4
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv'
metab_num=ncol(df.metab)-LABEL_NUM
df.raw=df.metab[,1:metab_num]
df.raw[which(df.metab<=0, arr.ind = T)]=1
df.raw=log2(df.raw)
mat=data.matrix(df.raw)
rownames(mat)=rownames(df.metab)
mat=t(mat)
label_col='k_cluster_3'
col_ha=df.metab[,label_col]
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2.5,by=0.01))
ComplexHeatmap::pheatmap(mat,
                                #col = col_fun,
                                name = "Relative level",
                                column_split = col_ha,
                                # row_split = row_ha,
                                # top_annotation=top_ha,
                                cluster_row_slices=F,
                                fontsize_row = 7, 
                                cluster_cols = T,
                                row_title_rot = 0,
                                fontsize_col = 5,
                                treeheight_row = 20,
                                row_title_gp = gpar(fontsize = 7),
                                legend_breaks=seq(-4,4,2),
                                scale="row",
                                breaks = bk
)
p<-ComplexHeatmap::pheatmap(mat,
                                #col = col_fun,
                                name = "Relative level",
                                # column_split = col_ha,
                                # row_split = row_ha,
                                # top_annotation=top_ha,
                                cluster_row_slices=F,
                                fontsize_row = 7, 
                                cluster_cols = T,
                                row_title_rot = 0,
                                fontsize_col = 5,
                                treeheight_row = 20,
                                row_title_gp = gpar(fontsize = 7),
                                legend_breaks=seq(-4,4,2),
                                scale="row",
                                breaks = bk
)
p
### Test Log2
df.metab.raw=df.metab[,1:metab_num]
df.metab.log=log2(df.metab.raw)
idx <- which(df.metab.raw <= 0, arr.ind = TRUE)
print(df.metab.raw[idx])

###


##### CPD Change Line #####

# Prepare Kegg Data
filepath.kegg='D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/kegg_table.csv'
# kegg_table<-read.csv(filepath.kegg, header= TRUE, check.names=F,row.names=1)
kegg_all=read.csv('D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/kegg_table.csv',check.names=FALSE)
kegg_table=kegg_all
kegg_table
kegg_table=kegg_table[order(abs(kegg_table$metabs_mean),decreasing = TRUE),]
kegg_table=kegg_table[order(kegg_table$p.adjust,decreasing = FALSE),]
kegg_table[1:4,]
df.cpd_kegg=kegg_table[kegg_table['comporison']=='C_0/C_2',c('Description','query')]
df.cpd_kegg=df.cpd_kegg[!duplicated(df.cpd_kegg),]
nrow(df.cpd_kegg)
rownames(df.cpd_kegg)=c(1:nrow(df.cpd_kegg))
df.cpd_kegg
# Prepare Metab Data
LABEL_NUM=4
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
df.metab[1:5,metab_num:(metab_num+4)]
metab_num=dim(df.metab)[2]-LABEL_NUM
df.metab.raw_data=df.metab[,1:metab_num]+abs(min(df.metab[,1:metab_num]))+1
class_label="k_cluster_3"
df.metab=cbind(df.metab.raw_data,df.metab[class_label])

# Scale data to make it positive
# df.metab.raw_data[df.metab.raw_data==0]=10
df.metab.log=log2(df.metab.raw_data)
idx <- which(df.metab.raw_data <= 0, arr.ind = TRUE)
print(df.metab.raw_data[idx])

df.metab.log=cbind(df.metab.log,df.metab[class_label])
# Try to draw figures
df.use=df.metab.log
draw_num=nrow(df.cpd_kegg)
dirpath='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/full_figure/cpd/'
# color_palette=rainbow(19)
# barplot(1:19,col=color_palette)
color_palette <- brewer.pal(n = 9, name = "Set1") 
color_palette
barplot(1:length(color_palette),col=color_palette)

color_palette


## Check 
df.use[,'Oxaloacetate']
df.test.results.cluster02[df.test.results.cluster02['Metabolites']=='Oxaloacetate']
df.test.results.cluster01[df.test.results.cluster01['Metabolites']=='Oxaloacetate']
tapply(df.test.cluster02[['Oxaloacetate']], df.test.cluster02[['k_cluster_3']], mean)
tapply(df.use[['Oxaloacetate']], df.use[['k_cluster_3']], mean)
df.test.cluster02[,'Oxaloacetate']
##

tapply(df.use[['Xanthine']], df.use[['k_cluster_3']], mean)
tmp_sds=tapply(df.use[['Xanthine']], df.use[['k_cluster_3']], get_sem)
tmp_sds
df.use=df.metab.log
df.use=rbind(df.use[df.use['k_cluster_3']==0,],df.use[df.use['k_cluster_3']==2,])
# df.use
for (fig_num in c(1:draw_num)){
    pathway_name=df.cpd_kegg[fig_num,'Description']
    print(pathway_name)
    # gene_ids = unlist(strsplit(df.cpd_kegg[i,'query'],'/'))
    # match(gene_ids,df.new$cpd_id)
    # queries=df.new[match(gene_ids,df.new$cpd_id),'query']
    queries=unlist(strsplit(df.cpd_kegg[fig_num,'query'],'/'))
    print(queries)
    # df.selected=filter_by_pvalue(df.use,queries,class_label,0.05)
    # print(dim(df.selected))
    # df.use=df.selected
    # queries=colnames(df.selected)[1:length(colnames(df.selected))-1]
    df.draw=data.frame(matrix(nrow = 0,ncol=4))
    colnames(df.draw)=c('Metabolite','k_cluster_3','mean','sd')
    for (row_index in c(1:length(queries))){
        metabolite=queries[row_index]
        # print(metabolite)
        tmp_means=tapply(df.use[[metabolite]], df.use[[class_label]], mean)
        # tmp_sds=tapply(df.use[[metabolite]], df.use[[class_label]], sd)
        tmp_sems=tapply(df.use[[metabolite]], df.use[[class_label]], get_sem)

        for (i in c(1:2)){
            # print(i)
            tmp_name=names(tmp_means)[i]
            df.draw[nrow(df.draw)+1,]=c(metabolite,tmp_name,as.numeric(tmp_means[[tmp_name]]),as.numeric(tmp_sds[[tmp_name]]))
        }
    }
    df.draw[,2:4]=lapply(df.draw[,2:4],smartConvert)
    df.draw
    # color_palette <- brewesr.pal(n = 20, name = "Set2") 
    color_palette
    p=ggplot(data = df.draw) +
        geom_line(aes(x = k_cluster_3, y = mean, colour = Metabolite))+ 
        geom_errorbar(aes(x = k_cluster_3,y=mean,ymin = mean - sd, ymax = mean + sd,col=Metabolite), width = 0.1)+
        # scale_y_continuous(limits = c(5,45))+
        scale_x_continuous(breaks=c(0,1,2),labels=c('0','1','2'))+
        scale_color_manual(values = color_palette)+
        ggtitle(pathway_name)+
        theme(plot.title = element_text(hjust = 0.5,size=20), legend.position = "right")
    print(p)
    pathway_name=str_replace_all(pathway_name,'/','_')
    filepath.cpd_fig=paste(dirpath,fig_num,'_',pathway_name,'.pdf',sep='')
    filepath.cpd_fig=str_replace_all(filepath.cpd_fig,' ','_')
    # ggsave(filepath.cpd_fig,p)
}
p
typeof(p)
filepath.cpd_fig=paste(dirpath,pathway_name,'.pdf',sep='')
filepath.cpd_fig=str_replace_all(filepath.cpd_fig,' ','_')
    
filepath.cpd_fig

### old
gene_ids = strsplit(df.cpd_kegg[1,'geneID'],'/')[[1]]
match(gene_ids,df.new$cpd_id)
queries=df.new[match(gene_ids,df.new$cpd_id),'query']
queries
LABEL_NUM=3
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20230902_cluster_lifeline/csvs/metab.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)

metab_num=dim(df.metab)[2]-LABEL_NUM
df.metab.raw_data=df.metab[,1:metab_num]
class_label="k_cluster_3"
df.metab=cbind(df.metab.raw_data,df.metab[class_label])

type1=0
type2=1
type3=2
df.test=rbind(df.metab[which(df.metab$k_cluster_3 == type1),], 
                 df.metab[which(df.metab$k_cluster_3 == type2),],
                 df.metab[which(df.metab$k_cluster_3 == type3),])
# df.test[,-1]                
df.test.results=data.frame(Metabolites = (colnames(df.test)[1:metab_num]))            
# df.test.results
df.test.results$kruskal=apply(df.test[, 1:metab_num], 2, function(x) unlist(kruskal.test(as.numeric(x), df.test$k_cluster_3,exact=TRUE)[3]))
df.test.results$kruskal_BH=p.adjust(df.test.results$kruskal,method="BH")
df.test.results$Metabolites[df.test.results$kruskal_BH<=5e-2]
# soft zero process
df.metab.raw_data[df.metab.raw_data==0]<-10
df.metab.raw_data[[tmp_cpd_names[47]]]

df.metab.log=log2(df.metab.raw_data)
df.metab.log[['Aconitate']]
df.metab.log[,1:3]
df.metab.log=cbind(df.metab.log,df.metab$k_cluster_3)
colnames(df.metab.log)[ncol(df.metab.log)]='k_cluster_3'
tmp_cpd_names=df.test.results$Metabolites[df.test.results$kruskal_BH<5e-2]
tmp_cpd_names
tmp_df=df.metab.log[,c(tmp_cpd_names,'k_cluster_3')]
tmp_df
df.metab.log$GTP
colnames(df.metab.log)
# Generate Drawing Dataframe
df.draw=data.frame(matrix(nrow = 0,ncol=4))
colnames(df.draw)=c('Metabolite','k_cluster_3','mean','sd')
df.draw
df.use=df.metab.log
# df.use
for (row_index in c(1:length(tmp_cpd_names))){
    metabolite=tmp_cpd_names[row_index]
    print(metabolite)
    tmp_means=tapply(df.use[[metabolite]], df.use[[class_label]], mean)
    tmp_sds=tapply(df.use[[metabolite]], df.use[[class_label]], sd)

    for (i in c(1:3)){
        # print(i)
        tmp_name=names(tmp_means)[i]
        df.draw[nrow(df.draw)+1,]=c(metabolite,tmp_name,as.numeric(tmp_means[[tmp_name]]),as.numeric(tmp_sds[[tmp_name]]))
    }

}
df.draw
# is.na(df.draw)
which(is.na(df.draw))

tmp_cpd_names[48]
df.metab[[tmp_cpd_names[47]]]
df.metab.log[[tmp_cpd_names[47]]]
df.draw[,2:4]=lapply(df.draw[,2:4],smartConvert)
typeof(df.draw[[4]])
df.draw[141,4]
df.draw[[4]]
df.draw
is.numeric(as.numeric(df.draw$sd))
is.numeric(df.draw$sd)
colnames(df.draw)
max <- max(df.draw$mean)
min <- min(df.draw$mean)
range <- 3*(max - min)
min=min-range
max=max+range
ggplot(data = df.draw[1:27,]) +
  geom_line(aes(x = k_cluster_3, y = mean, colour = Metabolite))+ 
  geom_errorbar(aes(x = k_cluster_3,y=mean,ymin = mean - sd, ymax = mean + sd,col=Metabolite), width = 0.1)+
  scale_y_continuous(limits = c(10,40))+
  scale_x_continuous(breaks=c(0,1,2),labels=c('0','1','2'))


df.draw[1:10,]

##### CPD Change Box #####

# Prepare Kegg Data
filepath.kegg='D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/kegg_table.csv'
# kegg_table<-read.csv(filepath.kegg, header= TRUE, check.names=F,row.names=1)
kegg_all=read.csv('D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/kegg_table.csv',check.names=FALSE)
kegg_table=kegg_all
kegg_table
kegg_table=kegg_table[order(abs(kegg_table$metabs_mean),decreasing = TRUE),]
kegg_table=kegg_table[order(kegg_table$p.adjust,decreasing = FALSE),]
kegg_table[1:4,]
df.cpd_kegg=kegg_table[kegg_table['comporison']=='C_0/C_2',c('Description','query')]
df.cpd_kegg=df.cpd_kegg[!duplicated(df.cpd_kegg),]
nrow(df.cpd_kegg)
rownames(df.cpd_kegg)=c(1:nrow(df.cpd_kegg))
df.cpd_kegg
# Prepare Metab Data
LABEL_NUM=4
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/csvs/metab.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
df.metab[1:5,metab_num:(metab_num+4)]
metab_num=dim(df.metab)[2]-LABEL_NUM
df.metab.raw_data=df.metab[,1:metab_num]+abs(min(df.metab[,1:metab_num]))+1
class_label="k_cluster_3"
df.metab=cbind(df.metab.raw_data,df.metab[class_label])

# Scale data to make it positive
# df.metab.raw_data[df.metab.raw_data==0]=10
df.metab.log=log2(df.metab.raw_data)
idx <- which(df.metab.raw_data <= 0, arr.ind = TRUE)
print(df.metab.raw_data[idx])

df.metab.log=cbind(df.metab.log,df.metab[class_label])
# Try to draw figures
df.use=df.metab.log
draw_num=nrow(df.cpd_kegg)
dirpath='D:/repositories/liver-cancer/tasks/Tissue/results/20231103/full_figure/cpd/'
# color_palette=rainbow(19)
# barplot(1:19,col=color_palette)
color_palette <- brewer.pal(n = 9, name = "Set1") 
color_palette
barplot(1:length(color_palette),col=color_palette)

color_palette


tapply(df.use[['Xanthine']], df.use[['k_cluster_3']], mean)
tmp_sds=tapply(df.use[['Xanthine']], df.use[['k_cluster_3']], get_sem)
tmp_sds
df.use=df.metab.log
df.use=rbind(df.use[df.use['k_cluster_3']==0,],df.use[df.use['k_cluster_3']==2,])
# df.use
for (fig_num in c(1:draw_num)){
    pathway_name=df.cpd_kegg[fig_num,'Description']
    print(pathway_name)
    # gene_ids = unlist(strsplit(df.cpd_kegg[i,'query'],'/'))
    # match(gene_ids,df.new$cpd_id)
    # queries=df.new[match(gene_ids,df.new$cpd_id),'query']
    queries=unlist(strsplit(df.cpd_kegg[fig_num,'query'],'/'))
    print(queries)
    # df.selected=filter_by_pvalue(df.use,queries,class_label,0.05)
    # print(dim(df.selected))
    # df.use=df.selected
    # queries=colnames(df.selected)[1:length(colnames(df.selected))-1]
    df.draw=data.frame(matrix(nrow = 0,ncol=4))
    colnames(df.draw)=c('Metabolite','k_cluster_3','mean','sd')
    for (row_index in c(1:length(queries))){
        metabolite=queries[row_index]
        # print(metabolite)
        tmp_means=tapply(df.use[[metabolite]], df.use[[class_label]], mean)
        # tmp_sds=tapply(df.use[[metabolite]], df.use[[class_label]], sd)
        tmp_sems=tapply(df.use[[metabolite]], df.use[[class_label]], get_sem)

        for (i in c(1:2)){
            # print(i)
            tmp_name=names(tmp_means)[i]
            df.draw[nrow(df.draw)+1,]=c(metabolite,tmp_name,as.numeric(tmp_means[[tmp_name]]),as.numeric(tmp_sds[[tmp_name]]))
        }
    }
    df.draw[,2:4]=lapply(df.draw[,2:4],smartConvert)
    df.draw
    # color_palette <- brewesr.pal(n = 20, name = "Set2") 
    color_palette
    p=ggplot(data = df.draw,aes_string(x = 'k_cluster_3', y = metabolite)) +
        geom_boxplot(aes(x = k_cluster_3, y = mean, colour = Metabolite))+ 
        geom_jitter(shape=16,position = position_jitter(0.2))+
        # scale_y_continuous(limits = c(5,45))+
        # scale_x_continuous(breaks=c(0,1,2),labels=c('0','1','2'))+
        scale_color_manual(values = color_palette)+
        ggtitle(pathway_name)+
        theme(plot.title = element_text(hjust = 0.5,size=20), legend.position = "right")
    print(p)
    pathway_name=str_replace_all(pathway_name,'/','_')
    filepath.cpd_fig=paste(dirpath,fig_num,'_',pathway_name,'.pdf',sep='')
    filepath.cpd_fig=str_replace_all(filepath.cpd_fig,' ','_')
    # ggsave(filepath.cpd_fig,p)
}
p
typeof(p)
filepath.cpd_fig=paste(dirpath,pathway_name,'.pdf',sep='')
filepath.cpd_fig=str_replace_all(filepath.cpd_fig,' ','_')
    
filepath.cpd_fig

