library(stringr)
# library(eoffice)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
# library(rJava)
# library(xlsx)
library(RColorBrewer)
library(circlize)
library(umap)
library(Rtsne)
library(ggrepel)
library(ggplot2)
library(clusterProfiler)
library(KEGGREST)
# library(ConsensusClusterPlus)
library(viridis)
library(survminer)
library(survival)
library(mixOmics)
library(patchwork)

packageVersion("rlang")    #查看指定R包版本
R.version
getwd()
setwd("/home/suh/workstation/MetabSubtype/tasks/MultiOmics/notebooks/")
source('utils.R')
##### Metab #####
# filepath.rna_counts='D:/repositories/liver-cancer/tasks/Tissue/results/20231202/csvs/lipid.csv'
filepath.rna_counts='/home/suh/workstation/MetabSubtype/tasks/MultiOmics/data/rna/counts.csv'
filepath.sample='/home/suh/workstation/MetabSubtype/tasks/Subtype/results/20240406/lipid_cluster_7e-2.csv'
df.rna_counts<-read.csv(filepath.rna_counts, header= TRUE, check.names=F,row.names=1)
df.sample<-read.csv(filepath.sample, header= TRUE, check.names=F,row.names=1)
df.rna_counts[1:5,1:5]
dim(df.rna_counts)
length(unique((df.rna_counts[,1])))
unique(df.rna_counts[,1])
df.rna_counts[,1]
df.rna_counts[df.rna_counts=='']=NA
df.sample[df.sample=='']=NA
df.sample[df.sample=='Neg']=NA
# df.sample=df.sample[!is.na(df.sample['Sample Name']),c('Sample Name','主要分型','生存时间分组','MT2结构')]
dim(df.sample)
df.sample=df.sample[rownames(df.rna_counts),]
df.sample=na.omit(df.sample)
dim(df.sample)
##### Data Cleaning #####
smallestGroupSize=min(nrow(df.sample[df.sample['nmf_2_clusters']==1,]),nrow(df.sample[df.sample['nmf_2_clusters']==2,]))
keep_cols=colnames(df.raw_counts)[colSums(df.raw_counts>20)>=smallestGroupSize]
df.rna_counts=df.rna_counts[,keep_cols]
metab_num=ncol(df.rna_counts)
df.raw_counts=df.rna_counts[,keep_cols]
df.rna_counts.log=log2(df.raw_counts)
dim(df.rna_counts)

##### Pathway #####
# MT2 pathway 

dirpath="D:/repositories/liver-cancer-tasks/Tissue/results/20240202/heatmap/MT2_pathway/"
pathway_name="Citrate cycle (TCA cycle)"
pathway_name="Purine metabolism"
pathway_name="Pyrimidine metabolism"
df.new[df.new$cpd_id=='C00655',]
df.new[df.new$query=='inosine',]
for (pathway_name in pathway_names.mt2){
  # Get Pathway CPDs
  metabs.kegg.tca=df.new[df.new$pathway_name==pathway_name,'query']
  metabs.kegg.tca
  metabs.ms.tca=get_ms_metabs(metabs.kegg.tca)
  
  metabs.ms.tca
  metabs.ms.tca=metabs.ms.tca[metabs.ms.tca!='NA']
  mat=as.matrix(df.pathway_fc[,'log2fc'])
  mat=as.matrix(df.pathway_fc[metabs.ms.tca,'log2fc',drop=FALSE])
  mat=mat[order(mat),,drop=FALSE]
  bk=c(-1,0,1)
  col_fun<-colorRamp2(
    bk,
    c("#1D91C0", "white", "#E31A1C"))
  p=ComplexHeatmap::pheatmap(mat,
                           name = pathway_name,
                           cluster_cols = F,
                           cluster_rows = F,
                           cellheight = 25,
                           cellwidth = 40,
                           col = col_fun
                           # scale="col"
                           
  )
  # print(p)
  # Sys.sleep(3)
  pathway_name=str_replace_all(pathway_name,'/','_')
  pathway_name=str_replace_all(pathway_name,',','_')
  pathway_name=str_replace_all(pathway_name,'\\(','_')
  pathway_name=str_replace_all(pathway_name,'\\)','_')
  filepath=paste(dirpath,pathway_name,'.pdf',sep='')
  filepath=str_replace_all(filepath,'  ',' ')
  pdf(filepath)
  draw(p)
  dev.off()

}
pathway_names.mt2
pathway_names.too_many=c("Purine metabolism","Pyrimidine metabolism")
for (pathway_name in pathway_names.too_many){
  # Get Pathway CPDs
  metabs.kegg.tca=df.new[df.new$pathway_name==pathway_name,'query']
  metabs.kegg.tca
  metabs.ms.tca=get_ms_metabs(metabs.kegg.tca)
  
  metabs.ms.tca
  metabs.ms.tca=metabs.ms.tca[metabs.ms.tca!='NA']
  mat=as.matrix(df.pathway_fc[,'log2fc'])
  mat=as.matrix(df.pathway_fc[metabs.ms.tca,'log2fc',drop=FALSE])
  mat=mat[order(mat),,drop=FALSE]
  bk=c(-1,0,1)
  col_fun<-colorRamp2(
    bk,
    c("#1D91C0", "white", "#E31A1C"))
  p=ComplexHeatmap::pheatmap(mat,
                             name = pathway_name,
                             cluster_cols = F,
                             cluster_rows = F,
                             cellheight = 25,
                             cellwidth = 40,
                             col = col_fun
                             # scale="col"
                             
  )
  print(p)
  Sys.sleep(30)
  pathway_name=str_replace_all(pathway_name,'/','_')
  pathway_name=str_replace_all(pathway_name,',','_')
  pathway_name=str_replace_all(pathway_name,'\\(','_')
  pathway_name=str_replace_all(pathway_name,'\\)','_')
  filepath=paste(dirpath,pathway_name,'.pdf',sep='')
  filepath=str_replace_all(filepath,'  ',' ')
  # pdf(filepath)
  draw(p)
  dev.off()
  
}
mat=as.matrix(df.pathway_fc[,'fc'])
mat=as.matrix(df.pathway_fc[metabs.ms.tca,'fc'])
bk=c(0,1,3)
col_fun<-colorRamp2(
    bk,
    c("#1D91C0", "white", "#E31A1C"))
# col_fun <- c(rep('#05C0C5', 100), colorRampPalette(c("#05C0C5", "white"))(50)[1:50], 
    # colorRampPalette(c("white", "#FF6C6C"))(50)[4:50], rep('#FF6C6C', 100))
dim(mat)
mat
confused_metabs=c('Uridine','Adenine','GMP','IMP','dAMP')
mat[confused_metabs,]
dirpath()
ComplexHeatmap::pheatmap(mat,
                         name = pathway_name,
                         cluster_cols = F,
                         cluster_rows = F,
                        col = col_fun
                            # scale="col"
                        
)
df.confused=cbind(df.use[,confused_metabs],df.use[,class_label,drop=FALSE])
draw_single_heatmap(df.confused,length(confused_metabs),class_label,1,FALSE )

ComplexHeatmap::pheatmap(log2(df.raw_counts[,confused_metabs]),
                            # col = col_fun,
                            name = pathway_name,
                            top_annotation=HeatmapAnnotation(df=df.sample[,class_label],col = color_list,
                                                             simple_anno_size=unit(4, "mm"),
                                                             #border = T,
                                                             gap = unit(1, "points"),
                                                             show_annotation_name = TRUE),
                         
                            cluster_row_slices=F,
                            fontsize_row = 7, 
                            cluster_cols = T,
                            row_title_rot = 0,
                            fontsize_col = 5,
                            treeheight_row = 20,
                            treeheight_col = 20,
                            row_title_gp = gpar(fontsize = 7),
                            scale="col",
                            # breaks = bk
                        )

##### Mean Aggregate #####  
dirpath="D:/repositories/liver-cancer-tasks/Tissue/results/20240202/heatmap/"
# setwd(dirpath)
class_label='CODEX主要亚型'
df.use_sample=df.sample
df.use=merge(df.raw_counts,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
df.rna_counts.mean=df.use %>% group_by(!!sym(class_label)) %>% summarise_all(mean)
df.rna_counts.mean=as.data.frame(df.rna_counts.mean)
df.rna_counts.mean=df.rna_counts.mean[!is.na(df.rna_counts.mean[,class_label]),]
rownames(df.rna_counts.mean)=df.rna_counts.mean[[class_label]]
df.rna_counts.mean=cbind(df.rna_counts.mean[,2:(1+metab_num)],df.rna_counts.mean[,1,drop=FALSE])
dim(df.rna_counts.mean)
metab_num
pvalue_cutoff=1
use_row_ha=FALSE
col_split=TRUE
unique(df.new$pathway_name)
for (pathway_name in pathway_names){
  pathway_metabs=df.new[df.new$pathway_name==pathway_name,'query']
  pathway_metabs=get_ms_metabs(pathway_metabs)
  pathway_metabs=pathway_metabs[pathway_metabs!='NA']
  using_metabs=c(pathway_metabs,class_label)
  p=draw_single_heatmap(df.rna_counts.mean[,using_metabs],length(pathway_metabs),class_label,pvalue_cutoff,use_row_ha )
  print(p)
  pathway_name=str_replace_all(pathway_name,'/','_')
  pathway_name=str_replace_all(pathway_name,',','_')
  pathway_name=str_replace_all(pathway_name,'\\(','_')
  pathway_name=str_replace_all(pathway_name,'\\)','_')
  filepath=paste(dirpath,pathway_name,'.pdf',sep='')
  filepath=str_replace_all(filepath,' ','')
  pdf(filepath)
  draw(p)
  dev.off()
  # pdf(filepath)
  # break
  # Sys.sleep(10)
  
}

draw_single_heatmap(df.rna_counts.mean,metab_num,class_label,pvalue_cutoff,use_row_ha)
Sys.sleep(3)

##### Heatmap #####

# Relabel Codex

class_label='nmf_2_clusters'
df.use_sample=df.sample
df.use_sample=df.use_sample[!is.na(df.use_sample[class_label]),]
df.use=merge(df.raw_counts,df.use_sample[,class_label,drop=FALSE],by='row.names')
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

df.use=df.use[df.use[[class_label]]!='Neg',]
pvalue_cutoff=5e-2
use_row_ha=FALSE
col_split=TRUE
draw_single_heatmap(df.use,metab_num,class_label,pvalue_cutoff,use_row_ha )
table(df.use_sample[class_label])

dim(df.test.results[df.test.results$wilcox_BH<5e-2,])

#####  Volcano #####

# Relabel Survival
# class_label='MT2结构'
# get significant quantification
volcano.test.results=pvalue_by_wilcox(df.use,metab_num,class_label)
volcano.test.results
class_label
# type1=target_group
# type2='others'
type1='high'
type2='low'
volcano.test.results$FC=apply(df.use[,1:metab_num], 2, 
                          function(x) 
                            mean(as.numeric(x[which(df.use[class_label] == type1)]))/
                            mean(as.numeric(x[which(df.use[class_label] == type2)])))
volcano.test.results
head(volcano.test.results)
# Log2 (FC)
volcano.test.results$log2FC <- log2(volcano.test.results[,'FC'])

dim(volcano.test.results[volcano.test.results$pvalue<5e-2,])
dim(volcano.test.results[(volcano.test.results$pvalue<5e-2)&(abs(volcano.test.results$log2FC)>log2(1.4)),])
##标记上调下调
pvalue_cutoff=5e-2
fc_up_cutoff=1.4
fc_down_cutoff=1/fc_up_cutoff

draw_volcano(volcano.test.results,pvalue_cutoff,fc_up_cutoff)


##### PCA ######

df.use=merge(df.raw_counts,df.sample[,c('CODEX主要亚型','MT2结构','codex_new'),drop=FALSE],by='row.names')
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
group_col='CODEX主要亚型'
for (target_group in c('A','N','Y','P')){
  col_name=paste(target_group,'_others',sep='')
  col_name
  df.use[col_name]=apply(df.use,1,get_codex_label, group_col=group_col,target_group=target_group)
  table(df.use$tmp_codex)
  
  draw_pca(df.use, col_name, metab_num, "Metab", "")
  Sys.sleep(3)
  draw_plsda(df.use, col_name, metab_num, "Metab", "")
  Sys.sleep(3)
  
}
target_group='A'
col_name=paste(target_group,'_others',sep='')
df.use[col_name]=apply(df.use,1,get_codex_label, group_col=group_col,target_group=target_group)
draw_pca(df.use, col_name, metab_num, "Metab", "")
draw_plsda(df.use, col_name, metab_num, "Metab", "")

draw_pca(df.use, "codex_new", metab_num, "Metab", '')
draw_plsda(df.use, "codex_new", metab_num, "Metab", '')
draw_plsda(df.use, "CODEX主要亚型", metab_num, "Metab", '')


dirpath='D:/repositories/liver-cancer-tasks/Tissue/results/20240202/pca/'
filepath=paste(dirpath,'codex_pca','.pdf',sep = '')
draw_pca(df.use, "CODEX主要亚型", metab_num, "Metab", filepath)
filepath=paste(dirpath,'codex_plsda','.pdf',sep = '')
draw_plsda(df.use, "CODEX主要亚型", metab_num, "Metab", filepath)
filepath=paste(dirpath,'mt2_pca','.pdf',sep = '')
draw_pca(df.use, "MT2结构", metab_num, "Metab", filepath)
filepath=paste(dirpath,'mt2_plsda','.pdf',sep = '')
draw_plsda(df.use, "MT2结构", metab_num, "Metab", filepath)

table(df.sample['codex_new'])

# Check PCA DF
df.raw=df.use
using_num=metab_num
df.pure <- df.raw[, 1:using_num]
df.pure <- df.pure %>% replace(is.na(.), 0)
colnames(df.pure)[1:3]
df_pca <- prcomp(df.pure, scale. = TRUE)
df_pca$x
class(df_pca)
df.pca_cpd=merge(df_pca$x,df.sample[,c('CODEX主要亚型','MT2结构','codex_new'),drop=FALSE],by='row.names')
write.csv(df.pca_cpd,'D:/repositories/liver-cancer-tasks/Tissue/results/20240202/pca/metab_pca_compounds.csv')

# Check PLS-DA
target_col='CODEX主要亚型'
df.raw=df.use
using_num=metab_num
df.pure <- df.raw[, 1:using_num]
df.pure <- df.pure %>% replace(is.na(.), 0)
group <- df.raw[[target_col]]
# group<-df.target$stage
df_plsda <- plsda(df.pure, group, ncomp = 4)
df_plsda
df_plsda$prop_expl_var$X

df.plsda_cpd=merge(df_plsda$variates,df.sample[,c('CODEX主要亚型'),drop=FALSE],by='row.names')



##### KEGG #####

# Build Dataset
class_label="kmeans_2_clusters"
df.use_sample=df.sample
df.use_sample=df.use_sample[!is.na(df.use_sample[class_label]),]

df.use=merge(df.raw_counts,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

# df.test.results.mt2=kegg_before_process(df.use,metab_num,class_label,types=c('1','2'))

# Before KEGG Analysis
df.test=df.use
types=c('1','2')
df.test.results=data.frame(geneid = (colnames(df.test)[1:metab_num]))
df.test.results$wilcox=apply(df.test[, 1:metab_num], 2,
                             function(x) unlist(wilcox.test(as.numeric(x) ~ df.test[[class_label]], data = df.test, exact = FALSE)[3]))
df.test.results$wilcox_BH=p.adjust(df.test.results$wilcox,method="BH")
df.test.results$FC <- apply(df.test[,1:metab_num], 2, 
                            function(x) 
                              mean(as.numeric(x[which(df.test[class_label] == types[2])]))/
                              mean(as.numeric(x[which(df.test[class_label] == types[1])])))
df.test.results$log2FC<- log2(df.test.results$FC)
dim(df.test.results)
# Read Stat Results and Do KEGG
filepath.stat.counts='/home/suh/workstation/MetabSubtype/tasks/MultiOmics/results/20240507/counts_stat.csv'
df.test.results<-read.csv(filepath.stat.counts, header= TRUE, check.names=F,row.names=1)
colnames(df.test.results)[which(names(df.test.results) =="geneid")] <- "cpd_id"


pvalue_cutoff=5e-2
fc_up_cutoff=1.2
fc_down_cutoff=1.2

df.test.results$is_use=(df.test.results$wilcox<=pvalue_cutoff)&((df.test.results$FC>=fc_up_cutoff)|(df.test.results$FC<=fc_down_cutoff))
df.test.results.up=df.test.results[df.test.results$log2FC>0,]
df.test.results.down=df.test.results[df.test.results$log2FC<0,]

# Choose Up or Down
used_ens_ids=df.test.results.up$cpd_id[df.test.results.up$is_use]
used_ens_ids=df.test.results.down$cpd_id[df.test.results.down$is_use]
#

print(paste0('used metab number is ',length(used_ens_ids)))
columns(org.Hs.eg.db)

used_entrz_ids=mapIds(org.Hs.eg.db,keys=used_ens_ids,keytype="ENSEMBL",column="ENTREZID")
kegg_result=enrichKEGG(
  used_entrz_ids,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
kegg_result
kegg_table <- na.omit(as.data.frame(kegg_result))
kegg_table$Description
kegg_table$Description=str_replace(kegg_table$Description,' - Homo sapiens \\(human\\)','')
# write.table(kegg_table,'D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/kegg_table.csv',sep=',',row.names = FALSE)
kegg_table$Description

kegg_table <- kegg_table[kegg_table$Count > 2,]

kegg_table$metabs_mean <- apply(kegg_table, 1, path_avelog2FC,df.test.results=df.test.results)
# calculate fold enrichment
kegg_table[1:5,3:8]
dim(kegg_table)
kegg_table$FoldEnrich <- apply(kegg_table, 1, 
                               function(x) as.numeric(unlist(strsplit(x[['GeneRatio']], '[/]'))[1])/
                                 as.numeric(unlist(strsplit(x[['GeneRatio']], '[/]'))[2])*
                                 as.numeric(unlist(strsplit(x[['BgRatio']], '[/]'))[2])/
                                 as.numeric(unlist(strsplit(x[['BgRatio']], '[/]'))[1]))

# Generate Draw Table
kegg_table_draw=kegg_table
unwanted_pathways=c('Glycerolipid metabolism','Butanoate metabolism','Lipoic acid metabolism','Arginine biosynthesis','Nitrogen metabolism')
kegg_table_draw=kegg_table_draw[!kegg_table_draw$Description %in% unwanted_pathways,]
names(kegg_table_draw)[names(kegg_table_draw) == "p.adjust"] <- "FDR"

pathway_names=unique(kegg_table_draw$Description)

kegg_table_draw

# Figure using FDR as X
kegg_table_draw=kegg_table_draw[order(kegg_table_draw$pvalue,decreasing = TRUE),]
kegg_table_draw=kegg_table_draw[(nrow(kegg_table_draw)-15):nrow(kegg_table_draw),]
dirpath='~/workstation/MetabSubtype/tasks/MultiOmics/results/20240504_singleomic/'
# write.csv(kegg_table.mt2,paste0(dirpath,'counts_pathway_up.csv',collapse = '/'))

kegg_table_draw$FDR
kegg_table_draw$FoldEnrich
p=ggplot(kegg_table_draw, aes(-log10(FDR), Description)) +
  geom_point(aes(fill = -log10(FDR), size = FoldEnrich), color = "black", shape = 21) +
  scale_y_discrete(limits = kegg_table_draw[['Description']])+
  scale_size(range = c(3, 15), breaks = c(seq(1,3,0.2))) +
  # scale_color_gradient(low="blue",high = "red")+
  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(seq(1,3,0.5)))+
  #  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.2,0.5, 1, 1.5)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  # coord_fixed(ratio = 0.4) +
  labs(x = "-Log10(FDR)", y = "", title = "RNA Upregulate", 
       fill = "-Log10(FDR)", size = "Enrich Ratio") +
  theme(axis.text.x = element_text(size = 20, face = "plain", colour = "black",angle = 45,hjust = 1), 
        axis.text.y = element_text(size = 18, face = "plain", colour = "black")) +
  theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
        legend.title = element_text(size = 10, face = "plain", colour = "black"), 
        legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))

# pdf(paste(dirpath,'rna_counts_up.pdf',sep = ''),12,7)
print(p)
dev.off()




##### CPD Line #####
# Prepare Metab

# # Norm
# metab_mean=apply(df.rna_counts.log,1,mean)
# metab_mean
# length(metab_mean)
# metab_std=apply(df.rna_counts.log,1,sd)
# df.scaled.mannul=(df.rna_counts.log-metab_mean)/metab_std
# # Check
# apply(df.scaled.fun,1,mean)
# apply(df.scaled.mannul,1,mean)
# # 
# df.rna_counts.scaled=df.scaled.mannul

idx <- which(df.raw_counts <= 0, arr.ind = TRUE)
print(df.raw_counts[idx])
# Generate df.use
class_label='CODEX主要亚型'
df.use_sample=df.sample
df.use=merge(df.rna_counts,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
df.use=df.use[!is.na(df.use[class_label]),]
# Prepare KEGG
df.use
# desied_pathways.codex=c('Citrate cycle (TCA cycle)','Arginine biosynthesis','Alanine, aspartate and glutamate metabolism')
desied_pathways.codex=c('Pyruvate metabolism','Pantothenate and CoA biosynthesis','Glyoxylate and dicarboxylate metabolism')
df.cpd_kegg=kegg_table.codex[((kegg_table.codex[['pvalue']]<=5e-2)),]
df.cpd_kegg$Description
# df.cpd_kegg=df.cpd_kegg[df.cpd_kegg$Description %in% desied_pathways.codex,]
# df.cpd_kegg=df.cpd_kegg[df.cpd_kegg$pvalue<=5e-2,]
# df.cpd_kegg=df.cpd_kegg[df.cpd_kegg$pvalue<=0.1,]

# Directly Use KEGG Draw Table
df.cpd_kegg=kegg_table.mt2
query_cell.pentose=paste('6-phospho-D-gluconate', 'D-erythrose-4-phosphate', 'Glucose-6-phosphate', 'D-glyceraldehdye-3-phosphate',sep='/')
query_cell.purine=paste('GMP', 'IMP', 'dAMP', 'Adenine',sep='/')
df.cpd_kegg=combine_and_deduplicate_queries(df.cpd_kegg )
df.cpd_kegg[df.cpd_kegg$Description=='Purine metabolism','query']=query_cell.purine
df.cpd_kegg[df.cpd_kegg$Description=='Pentose phosphate pathway','query']=query_cell.pentose

df.cpd_kegg[,c('Description','query')]
draw_num=nrow(df.cpd_kegg)
pathway_name='Citrate cycle (TCA cycle)'
df.cpd_kegg[['Description']]
which(df.cpd_kegg[['Description']]==pathway_name)
draw_kegg_cpd<-function(df.use,df.cpd_kegg,pathway_name){
  color_palette <- brewer.pal(n = 9, name = "Set1")
  print(pathway_name)
  queries=unlist(strsplit(df.cpd_kegg[which(df.cpd_kegg[['Description']]==pathway_name),'query'],'/'))
  queries=unname(unlist(sapply(queries,function(search_metabolite){(colnames(df.raw_counts)[grep(search_metabolite,colnames(df.raw_counts),ignore.case = TRUE)][1])})))
  print(queries)
  # Try to norm in loop
  
  df.rna_counts.scaled=df.rna_counts.log
  df.rna_counts.scaled=df.rna_counts.scaled[,queries]
  metab_mean=apply(df.rna_counts.log,1,mean)
  metab_mean
  length(metab_mean)
  metab_std=apply(df.rna_counts.log,1,sd)
  df.rna_counts.scaled=(df.rna_counts.scaled-metab_mean)/metab_std
  
  class_label='CODEX主要亚型'
  df.use_sample=df.sample
  df.use=merge(df.rna_counts.scaled,df.use_sample[,class_label,drop=FALSE],by="row.names")
  row.names(df.use)=df.use[[1]]
  df.use=df.use[,2:ncol(df.use)]
  df.use=df.use[!is.na(df.use[class_label]),]
  #
  df.draw=data.frame(matrix(nrow = 0,ncol=4))
  colnames(df.draw)=c('Metabolite','k_cluster_3','mean','sem')
  for (row_index in c(1:length(queries))){
    metabolite=queries[row_index]
    # print(metabolite)
    tmp_means=tapply(df.use[[metabolite]], df.use[[class_label]], mean)
    # tmp_sds=tapply(df.use[[metabolite]], df.use[[class_label]], sd)
    tmp_sems=tapply(df.use[[metabolite]], df.use[[class_label]], get_sem)
    
    for (i in c(1:4)){
      # print(i)
      tmp_name=names(tmp_means)[i]
      df.draw[nrow(df.draw)+1,]=c(metabolite,tmp_name,as.numeric(tmp_means[[tmp_name]]),as.numeric(tmp_sems[[tmp_name]]))
    }
  }
  df.draw[,2:4]=lapply(df.draw[,2:4],smartConvert)
  df.draw
  # color_palette <- brewesr.pal(n = 20, name = "Set2") 
  color_palette
  blank_margin=0.01
  ylim_down=min(df.draw$mean)-blank_margin
  ylim_up=max(df.draw$mean)+blank_margin
  p=ggplot(data=df.draw,aes(x = k_cluster_3, y = mean, colour = Metabolite,group=Metabolite),
           position=position_dodge(0.05))+
    geom_point(aes(fill=Metabolite),size=4,pch=16,position=position_dodge(0.05))+
    geom_line(size=1.1,position=position_dodge(0.05))+
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem,colour = Metabolite,group=Metabolite), width=.2,size=1,
                  position=position_dodge(0.05))+
    # scale_x_discrete(labels=c('high','low'),breaks=c('high','low'))+
    ylim(ylim_down,ylim_up)+
    theme_bw() +
    labs(x='',y='Normalized')+
    ggtitle(pathway_name)+
    theme(plot.title = element_text(hjust = 0.5,size=20), legend.position = "right")+
    ##固定长宽比例不随着改变而改变
    theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1, panel.grid = element_blank())+
    ##固定字体大小  
    theme(axis.title.x = element_text(size = 17, hjust = 0.5, face = "plain"),
          axis.title.y = element_text(size = 17, face = "plain"), 
          axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 15, face = "plain", colour = "black")) 
    ##调整legend字体  
    theme(legend.text = element_text(size = 15), 
          legend.title = element_text(size = 15), 
          legend.spacing.y = unit(0.2, 'cm'), 
          legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm')) 
  # theme(panel.border = element_rect(size=2))
  print(p)
}
draw_kegg_cpd(df.use,df.cpd_kegg,'Purine metabolism')
# df.cpd_kegg=kegg_table.mt2
df.cpd_kegg
for (fig_num in c(1:draw_num)){
  pathway_name=df.cpd_kegg[fig_num,'Description']
  draw_kegg_cpd(df.use,df.cpd_kegg ,pathway_name )
  Sys.sleep(7)
  
}
df.cpd_kegg
dirpath='D:/repositories/liver-cancer-tasks/Tissue/results/20240202/line/mt2/'
color_palette <- brewer.pal(n = 9, name = "Set1")
class_label='MT2结构'
df.cpd_kegg

df.test.results.mt2=kegg_before_process(df.use,ncol(df.use)-1,class_label,types=c('high','low'))
fig_num=nrow(df.cpd_kegg)
ps=list()
for (fig_num in c(1:nrow(df.cpd_kegg))){
    pathway_name=df.cpd_kegg[['Description']][fig_num]
    print(pathway_name)
    queries=unlist(strsplit(df.cpd_kegg[which(df.cpd_kegg[['Description']]==pathway_name),'query'],'/'))
    queries=unname(unlist(sapply(queries,function(search_metabolite){(colnames(df.raw_counts)[grep(search_metabolite,colnames(df.raw_counts),ignore.case = TRUE)][1])})))
    print(queries)
    # Try to norm in loop
    
    df.rna_counts.scaled=df.rna_counts.log
    df.rna_counts.scaled=df.rna_counts.scaled[,queries,drop=FALSE]
    metab_mean=apply(df.rna_counts.log,1,mean)
    metab_mean
    length(metab_mean)
    metab_std=apply(df.rna_counts.log,1,sd)
    df.rna_counts.scaled=(df.rna_counts.scaled-metab_mean)/metab_std
    
    df.use_sample=df.sample
    df.use=merge(df.rna_counts.scaled,df.use_sample[,class_label,drop=FALSE],by="row.names")
    row.names(df.use)=df.use[[1]]
    df.use=df.use[,2:ncol(df.use)]
    df.use=df.use[!is.na(df.use[class_label]),]
    #
    df.draw=data.frame(matrix(nrow = 0,ncol=4))
    colnames(df.draw)=c('Metabolite','k_cluster_3','mean','sem')
    for (row_index in c(1:length(queries))){
      metabolite=queries[row_index]
      # print(metabolite)
      tmp_means=tapply(df.use[[metabolite]], df.use[[class_label]], mean)
      # tmp_sds=tapply(df.use[[metabolite]], df.use[[class_label]], sd)
      tmp_sems=tapply(df.use[[metabolite]], df.use[[class_label]], get_sem)
      
      for (i in c(1:2)){
        # print(i)
        tmp_name=names(tmp_means)[i]
        df.draw[nrow(df.draw)+1,]=c(metabolite,tmp_name,as.numeric(tmp_means[[tmp_name]]),as.numeric(tmp_sems[[tmp_name]]))
      }
    }
    df.draw[,2:4]=lapply(df.draw[,2:4],smartConvert)
    df.draw
    # color_palette <- brewesr.pal(n = 20, name = "Set2") 
    color_palette
    blank_margin=0.02
    ylim_down=min(df.draw$mean)-max(df.draw$sem)-blank_margin
    ylim_up=max(df.draw$mean)+max(df.draw$sem)+blank_margin
    p=ggplot(data=df.draw,aes(x = k_cluster_3, y = mean, colour = Metabolite,group=Metabolite),
             position=position_dodge(0.05))+
      geom_point(aes(fill=Metabolite),size=3,pch=16,position=position_dodge(0.05))+
      geom_line(size=1.1,position=position_dodge(0.05))+
      geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem,colour = Metabolite,group=Metabolite), width=.2,size=1.1,
                    position=position_dodge(0.05))+
      # scale_x_discrete(labels=c('high','low'),breaks=c('high','low'))+
      ylim(ylim_down,ylim_up)+
      theme_bw() +
      labs(x='',y='Normalized')+
      ggtitle(str_wrap(pathway_name,width = 30))+
      theme(plot.title = element_text(hjust = 0.5,size=17), legend.position = "right")+
      ##固定长宽比例不随着改变而改变
      theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"), aspect.ratio = 1, panel.grid = element_blank())+
      ##固定字体大小  
      theme(axis.title.x = element_text(size = 16, hjust = 0.5, face = "plain"),
            axis.title.y = element_text(size = 16, face = "plain"), 
            axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
            axis.text.x = element_text(size = 15, face = "plain", colour = "black")) +
      ##调整legend字体  
      theme(legend.text = element_text(size = 15), 
            legend.title = element_text(size = 15), 
            legend.spacing.y = unit(0.2, 'cm'), 
            legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm')) 
      # theme(panel.border = element_rect(size=2))
    ps[[fig_num]]=p
    # print(p)
    # Sys.sleep(10)
    pathway_name=str_replace_all(pathway_name,'/','_')
    pathway_name=str_replace_all(pathway_name,',','_')
    filepath.cpd_fig=paste(dirpath,pathway_name,'.pdf',sep='')
    filepath.cpd_fig=str_replace_all(filepath.cpd_fig,' ','_')
    '
    pdf(filepath.cpd_fig)
    print(p)
    dev.off()
    '
}
res<-wrap_plots(ps,nrow=4,guides='keep')
pdf(paste0(dirpath,'/all_in_one.pdf'),25,25)
print(res)
dev.off()
res
##### Specific CPD Line #####

class_label='MT2结构'
df.use_sample=df.sample
df.use_sample=df.use_sample[df.use_sample[['CODEX主要亚型']] %in% c('A'),]
df.use=merge(df.rna_counts.scaled,df.use_sample[,class_label,drop=FALSE],by='row.names')
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

pathway_names=c('Pentose Phosphate Pathway','Purine and Pyrimidine')
pathway_1_cpds=c('6-phospho-D-gluconate', 'D-erythrose-4-phosphate', 'Glucose-6-phosphate', 'D-glyceraldehdye-3-phosphate')
pathway_2_cpds=c('Uridine', 'Adenine', 'GMP', 'IMP', 'Cytidine', 'dAMP')
pathway_cpds=list(pathway_1_cpds,pathway_2_cpds)

## Test
df.rna_counts.scaled=df.rna_counts.log
df.rna_counts.scaled=df.rna_counts.scaled[,queries]
df.use_sample=df.sample
df.use=merge(df.rna_counts.scaled,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
df.use=df.use[!is.na(df.use[class_label]),]
df.use[order(df.use$MT2结构),]
for(metabolite in queries){
  print(tapply(df.use[[metabolite]], df.use[[class_label]], mean))
}
##

# For Loop 
pathway_index=1
queries=pathway_cpds[[pathway_index]]
queries
# Scale
df.rna_counts.scaled=df.rna_counts.log
# df.rna_counts.scaled=df.rna_counts.scaled[,queries]
metab_mean=apply(df.rna_counts.scaled,1,mean)
metab_mean
length(metab_mean)
metab_std=apply(df.rna_counts.scaled,1,sd)
df.rna_counts.scaled=(df.rna_counts.scaled-metab_mean)/metab_std
# Generate DF.use
df.use_sample=df.sample
df.use=merge(df.rna_counts.scaled,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
df.use=df.use[!is.na(df.use[class_label]),]
df.use[order(df.use$MT2结构),]
# generate df.draw
df.draw=data.frame(matrix(nrow = 0,ncol=4))
colnames(df.draw)=c('Metabolite','k_cluster_3','mean','sem')
for (row_index in c(1:length(queries))){
  metabolite=queries[row_index]
  # print(metabolite)
  tmp_means=tapply(df.use[[metabolite]], df.use[[class_label]], mean)
  # tmp_sds=tapply(df.use[[metabolite]], df.use[[class_label]], sd)
  tmp_sems=tapply(df.use[[metabolite]], df.use[[class_label]], get_sem)
  
  for (i in c(1:2)){
    # print(i)
    tmp_name=names(tmp_means)[i]
    df.draw[nrow(df.draw)+1,]=c(metabolite,tmp_name,as.numeric(tmp_means[[tmp_name]]),as.numeric(tmp_sems[[tmp_name]]))
  }
}
df.draw[,2:4]=lapply(df.draw[,2:4],smartConvert)
df.draw
# color_palette <- brewesr.pal(n = 20, name = "Set2") 
color_palette
blank_margin=0.01
ylim_down=min(df.draw$mean)-max(df.draw$sem)-blank_margin
ylim_up=max(df.draw$mean)+max(df.draw$sem)+blank_margin
p=ggplot(data=df.draw,aes(x = k_cluster_3, y = mean, colour = Metabolite,group=Metabolite),
         position=position_dodge(0.05))+
  geom_point(aes(fill=Metabolite),size=4,pch=16,position=position_dodge(0.05))+
  geom_line(size=1.1,position=position_dodge(0.05))+
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem,colour = Metabolite,group=Metabolite), width=.2,size=1,
                position=position_dodge(0.05))+
  # scale_x_discrete(labels=c('high','low'),breaks=c('high','low'))+
  ylim(ylim_down,ylim_up)+
  theme_bw() +
  labs(x='',y='Normalized',title = pathway_names[pathway_index])+
  theme(plot.title = element_text(hjust = 0.5,size=20), legend.position = "right")+
  ##固定长宽比例不随着改变而改变
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1, panel.grid = element_blank())+
  ##固定字体大小  
  theme(axis.title.x = element_text(size = 17, hjust = 0.5, face = "plain"),
        axis.title.y = element_text(size = 17, face = "plain"), 
        axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
        axis.text.x = element_text(size = 15, face = "plain", colour = "black")) +
  ##调整legend字体  
  theme(legend.text = element_text(size = 15), 
        legend.title = element_text(size = 15), 
        legend.spacing.y = unit(0.2, 'cm'), 
        legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm')) 
# theme(panel.border = element_rect(size=2))
pdf('D:/repositories/liver-cancer-tasks/Tissue/results/20240202/line/mt2/specific_pentose.pdf',9,9)
print(p)
dev.off()
##### Specific CPD BOX #####
pathway_names=c('Pentose Pathway','Purine and Pyrimidine')
pathway_1_cpds=c('6-phospho-D-gluconate', 'D-erythrose-4-phosphate', 'Glucose-6-phosphate', 'D-glyceraldehdye-3-phosphate')
pathway_2_cpds=c('Uridine', 'Adenine', 'GMP', 'IMP', 'Cytidine', 'dAMP')


codex_comporisons= list( c("A", "N"), c("A", "P"), c("A", "Y"), c("N", "P"), c("N", "Y"), c("Y", "P") )

for(feature in significant_metabolites){
  p=ggplot(df.use, aes(x = !!sym(class_label), y = !!sym(feature))  )+
    geom_boxplot(width = 0.4,aes_string(color=class_label)) +  # Remove outlier points
    geom_jitter(width = 0.1,aes_string(color=class_label))+
    # stat_compare_means()+
    stat_compare_means(comparisons = codex_comporisons,label = "p.signif",method='wilcox.test')+
    # theme_bw() +
    labs(x = "", y = "Normalized", title = feature,color='') +
    theme(axis.title.x = element_text(size = 17, hjust = 0.5, face = "plain"),
          axis.title.y = element_text(size = 17, face = "plain"), 
          axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 15, face = "plain", colour = "black")) +
    theme(plot.title = element_text(size=20,hjust = 0.5))+
    ##调整legend字体  
    theme(legend.text = element_text(size = 15), 
          legend.title = element_text(size = 15), 
          legend.spacing.y = unit(0.2, 'cm'), 
          legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'))+
    theme(legend.position="none") # remove legend
  print(p)
  Sys.sleep(10)
  
}


##### CPD Box #####

# Df
class_label='CODEX主要亚型'
df.use_sample=df.sample
df.use=merge(df.rna_counts.scaled,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
df.use=df.use[!is.na(df.use[class_label]),]

# Selection
df.cpd_kegg
metabolites=c()
for (queries in unlist(df.cpd_kegg$query)){
  queries=unlist(strsplit(queries,'/'))
  queries=unname(unlist(sapply(queries,function(search_metabolite){(colnames(df.use)[grep(search_metabolite,colnames(df.use),ignore.case = TRUE)][1])})))
  metabolites=append(metabolites,unlist(strsplit(queries,'/')))
}
metabolites
metabolite_counts=table(metabolites)
significant_metabolites=names(metabolite_counts[metabolite_counts>1])
significant_metabolites

# Draw
codex_comporisons= list( c("A", "N"), c("A", "P"), c("A", "Y"), c("N", "P"), c("N", "Y"), c("Y", "P") )
# write.csv(df.use,'D:/repositories/liver-cancer-tasks/Tissue/results/20231202/csvs/tmp_use.csv')
for(feature in significant_metabolites){
  p=ggplot(df.use, aes(x = !!sym(class_label), y = !!sym(feature))  )+
    geom_boxplot(width = 0.4,aes_string(color=class_label)) +  # Remove outlier points
    geom_jitter(width = 0.1,size=1.5,aes_string(color=class_label))+
    # stat_compare_means()+
    stat_compare_means(aes(alpha=p.value),comparisons = codex_comporisons,label = "p.signif",method='wilcox.test',hide.ns = TRUE)+
    theme_bw() +
    labs(x = "", y = "Normalized", title = feature,color='') +
    theme(axis.title.x = element_text(size = 17, hjust = 0.5, face = "plain"),
          axis.title.y = element_text(size = 17, face = "plain"), 
          axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 15, face = "plain", colour = "black")) +
    theme(plot.title = element_text(size=20,hjust = 0.5))+
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), aspect.ratio = 1, panel.grid = element_blank())+
    ##调整legend字体  
    theme(legend.text = element_text(size = 15), 
          legend.title = element_text(size = 15), 
          legend.spacing.y = unit(0.2, 'cm'), 
          legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'))+
    theme(legend.position="none") # remove legend
  print(p)
  Sys.sleep(10)
  
}


##### Grep Analysis #####
dim(df.raw_counts)
metab_num
df.grep_lipid=get_grep_df(df.raw_counts,metab_num)
lipid_head_num=ncol(df.grep_lipid)-(ncol(df.rna_counts)-metab_num)
headgroups=colnames(df.grep_lipid)
headgroups
"
# Select Batch
df.use_sample=df.sample
df.use_sample=df.use_sample[df.use_sample[['CODEX主要亚型']] %in% c('A'),]
# df.use_sample=df.use_sample[df.use_sample[['TMN']] %in% c('II'),]

low_cutoff=24
high_cutoff=60
df.use_sample$tmp_survival=apply(df.use_sample[,c('os','oss')],1,get_survival_label, low_cutoff=low_cutoff,high_cutoff=high_cutoff)
df.use_sample[(df.use_sample['oss']==1)&(df.use_sample['os']<=12),]
"

# Relabel Codex

df.use_sample=df.sample
df.use_sample=df.use_sample[!is.na(df.use_sample['CODEX主要亚型']),]
df.use_sample['CODEX主要亚型']
group_col='CODEX主要亚型'

target_group='N'
df.use_sample$tmp_codex=apply(df.use_sample,1,get_codex_label, group_col=group_col,target_group=target_group)
table(df.use_sample$tmp_codex)

# Relabel Survival
class_label='tmp_codex'


df.use=merge(df.grep_lipid,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
colnames(df.use)

df.use=df.use[!is.na(df.use[[class_label]]),]
df.use=cbind(scale(df.use[,1:use_num]),df.use[,(use_num+1):ncol(df.use),drop=FALSE])
use_num=ncol(df.use)-1

# Calculate p-value
pvalues=sapply(X=(colnames(df.use)[1:use_num]),FUN=function(feature){
    wilcox.test(df.use[[feature]] ~ df.use[[class_label]],exact=FALSE)$p.value
})
pvalues
significant_headgroups=names(pvalues[pvalues<=5e-2])
headgroups=names(pvalues[pvalues<=1e-1])
significant_headgroups
headgroups
draw_boxplot(df.use,significant_headgroups[1],class_label)
draw_boxplot(df.use,significant_headgroups[2],class_label)
draw_boxplot(df.use,headgroups[1],class_label)

for (headgroup in headgroups){
  draw_boxplot(df.use,headgroup,class_label)
  Sys.sleep(3)
  
}

# CODEX with Signif
class_label='MT2结构'
class_label='CODEX主要亚型'

df.use_sample=df.sample
df.use=merge(df.grep_lipid,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
colnames(df.use)
df.use=df.use[!is.na(df.use[class_label]),]
df.use=df.use[df.use[[class_label]]!='0',]
df.use=df.use[df.use[[class_label]]!='Neg',]
use_num=ncol(df.use)-1

# Scale Data
df.use=cbind(scale(df.use[,1:use_num]),df.use[,(use_num+1):ncol(df.use),drop=FALSE])


df.use
codex_comporisons= list( c("A", "N"), c("A", "P"), c("A", "Y"), c("N", "P"), c("N", "Y"), c("Y", "P") )
 
feature=colnames(df.use)[1]
colnames(df.use)
feature
for(feature in colnames(df.use)[1:use_num]){
  p=ggplot(df.use, aes_string(x = class_label, y = feature) )+
    geom_boxplot(width = 0.4,aes_string(color=class_label)) +  # Remove outlier points
    geom_jitter(width = 0.1,aes_string(color=class_label))+
    # stat_compare_means()+
    stat_compare_means(comparisons = codex_comporisons,label = "p.signif",method='wilcox.test')+
    # theme_bw() +
    labs(x = "", y = "Normalized", title = feature,color='') +
    theme(axis.title.x = element_text(size = 17, hjust = 0.5, face = "plain"),
            axis.title.y = element_text(size = 17, face = "plain"), 
            axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
            axis.text.x = element_text(size = 15, face = "plain", colour = "black")) +
    theme(plot.title = element_text(size=20,hjust = 0.5))+
    ##调整legend字体  
    theme(legend.text = element_text(size = 15), 
            legend.title = element_text(size = 15), 
            legend.spacing.y = unit(0.2, 'cm'), 
            legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'))+
    theme(legend.position="none") # remove legend
  print(p)
  Sys.sleep(3)

}
