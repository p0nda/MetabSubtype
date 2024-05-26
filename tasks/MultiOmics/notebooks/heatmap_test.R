
##### new #####
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
source('~/workstation/MetabSubtype/tasks/MultiOmics/notebooks/utils.R')

##### RNA 100  #####

LABEL_NUM=3
path.metab='D:/repositories/liver-cancer/tasks/Tissue/results/RNA/20231109/batch_1.csv'
path.I_pvalue='D:/repositories/liver-cancer/tasks/Tissue/results/RNA/20231109/I_pvalue.csv'
path.Cluster2_pvalue='D:/repositories/liver-cancer/tasks/Tissue/results/RNA/20231109/Cluster2_pvalue.csv'

df.metab<-read.csv(path.metab, header= TRUE, check.names=F,row.names=1)
df.I_pvalue<-read.csv(path.I_pvalue, header= TRUE, check.names=F,row.names=1)
df.Cluster2_pvalue<-read.csv(path.Cluster2_pvalue, header= TRUE, check.names=F,row.names=1)

I_significant_cols=rownames(df.I_pvalue)[which(df.I_pvalue$pvalue<=5e-2)]
Cluster2_significant_cols=rownames(df.Cluster2_pvalue)[which(df.Cluster2_pvalue$pvalue<=5e-2)]
(ncol(df.metab)-1)
colnames(df.metab)
class_label='nmf_cluster_4'
df.metab[class_label]
# Select cols
dim(df.metab)
df.metab=df.metab[,c(intersect(I_significant_cols,Cluster2_significant_cols),class_label)]
# df.metab=df.metab[,c(I_significant_cols,class_label)]
# df.metab=df.metab[,c(Cluster2_significant_cols,class_label)]
dim(df.metab)
#
metab_num=ncol(df.metab)-LABEL_NUM
pvalue_cutoff=1
use_row_ha=FALSE
draw_single_heatmap(df.metab,metab_n`um,class_label,pvalue_cutoff,use_row_ha )

######## draw single heatmap #########
df.lipid<-read.csv("D:/repositories/liver-cancer/tasks/eno_gra/data/batch_2/tissue/lipid.csv", header= TRUE, check.names=F,row.names=1)
class_label='molecular'
df.lipid[class_label]
lipid_num=ncol(df.lipid)-LABEL_NUM
pvalue_cutoff=0.05
use_row_ha=TRUE
draw_single_heatmap(df.lipid,lipid_num,class_label,pvalue_cutoff,use_row_ha )



######## Batch2 Data #########
source_data_path="D:/repositories/liver-cancer/tasks/Tissue/data/20230830_delete_outliners/metab.csv"
df.metab<-t(read.csv(source_data_path, header= TRUE, check.names=F,row.names=1))
# df.metab=t(df.metab)
df.sample<-read.csv("D:/repositories/liver-cancer/tasks/Tissue/data/SCLC_tissue/第二批/sample_info.csv", header= TRUE, check.names=F,row.names=2)
df.sample=df.sample['TMN']
rownames(df.metab)
rownames(df.sample)
colnames(df.sample)
dim(df.metab)
df.metab=merge(df.metab,df.sample,"row.names")
dim(df.metab)
rownames(df.metab)=df.metab[[1]]
df.metab=df.metab[,-1]
dim(df.metab)

print(colnames(df.metab))
LABEL_NUM=1
class_label='TMN'
df.metab[class_label]
metab_num=ncol(df.metab)-LABEL_NUM
pvalue_cutoff=0.05
use_row_ha=FALSE
draw_single_heatmap(df.metab,metab_num,class_label,pvalue_cutoff,use_row_ha )

##### 20231209 #####

LABEL_NUM=3
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231202/csvs/metab.csv'
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
class_label='MT2结构'
!is.na(df.metab[[class_label]])
df.metab[[class_label]]
df.metab=df.metab[!df.metab[[class_label]]=='',]
# class_label='batch'
df.metab[class_label]
metab_num=ncol(df.metab)-LABEL_NUM
pvalue_cutoff=5e-2
use_row_ha=FALSE
col_split=TRUE
draw_single_heatmap(df.metab,metab_num,class_label,pvalue_cutoff,use_row_ha )

##### 20231217 Mutiple Label #####

LABEL_NUM=3
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231202/csvs/metab.csv'
filepath.sample='D:/repositories/liver-cancer/tasks/Tissue/data/20231202_batch_sum/new_sample.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
df.sample<-read.csv(filepath.sample, header= TRUE, check.names=F)
df.metab[df.metab=='']=NA
df.sample[df.sample=='']=NA
df.sample=df.sample[!is.na(df.sample['Sample Name']),c('Sample Name','主要分型','生存时间分组','MT2结构')]
df.sample['Sample Name']
row.names(df.sample)=df.sample[['Sample Name']]
df.sample=df.sample[rownames(df.metab),]
dim(df.metab)
metab_num=ncol(df.metab)-LABEL_NUM
df.raw_metab=df.metab[,1:metab_num]

ha_cols=c("主要分型","MT2结构")
df.ha=df.sample[,ha_cols]
df.ha=df.ha[rowSums(is.na(df.ha[,ha_cols]))==0,,drop=FALSE]
df.ha=df.ha[df.ha[['主要分型']] %in% c('A'),]
df.ha=df.ha[,c('MT2结构'),drop=FALSE]
df.ha
df.use=merge(df.raw_metab,df.ha,by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
dim(df.use)
dim(df.metab)
colnames(df.use)

(ncol(df.use)-1)
class_label='主要分型'
!is.na(df.use[[class_label]])
df.use[[class_label]]
# class_label='batch'
pvalue_cutoff=5
use_row_ha=FALSE
col_split=TRUE
draw_single_heatmap(df.use,metab_num,class_label,pvalue_cutoff,use_row_ha )

##### 20231222 分型TMN #####
# df.sample$survival=apply(df.sample[,c('os','oss')],1,get_survival_label, low_cutoff=12,high_cutoff=48)
##

LABEL_NUM=3
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231202/csvs/lipid.csv'
filepath.sample='D:/repositories/liver-cancer/tasks/Tissue/data/20231202_batch_sum/new_sample.csv'
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
dim(df.metab)
metab_num=ncol(df.metab)-LABEL_NUM
df.raw_metab=df.metab[,1:metab_num]


# Generate using DF

class_label='MT2结构'

# Select Batch
df.use_sample=df.sample
# df.use_sample=df.use_sample[df.use_sample[['CODEX主要亚型']] %in% c('P'),]
# df.use_sample=df.use_sample[df.use_sample[['TMN']] %in% c('II'),]
df.use_sample

# Relabel Codex
"
df.use_sample=df.sample
df.use_sample=df.use_sample[!is.na(df.use_sample['CODEX主要亚型']),]
df.use_sample['CODEX主要亚型']
group_col='CODEX主要亚型'

target_group='P'
df.use_sample$tmp_codex=apply(df.use_sample,1,get_codex_label, group_col=group_col,target_group=target_group)
table(df.use_sample$tmp_codex)
"

# Relabel Survival
"
low_cutoff=24
high_cutoff=60
df.use_sample$tmp_survival=apply(df.use_sample[,c('os','oss')],1,get_survival_label, low_cutoff=low_cutoff,high_cutoff=high_cutoff)
df.use_sample[(df.use_sample['oss']==1)&(df.use_sample['os']<=12),]
table(df.use_sample$tmp_survival)
"
# merge
df.use=merge(df.raw_metab,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
!is.na(df.use[[class_label]])
df.use[[class_label]]
# class_label='batch'
df.use=df.use[df.use[[class_label]]!='0',]
df.use=df.use[df.use[[class_label]]!='Neg',]
pvalue_cutoff=5e-2
use_row_ha=FALSE
col_split=TRUE
draw_single_heatmap(df.use,metab_num,class_label,pvalue_cutoff,use_row_ha )
table(df.use_sample$MT2结构)

##### 20240105 Lipid #####

LABEL_NUM=3
filepath.metab='D:/repositories/liver-cancer/tasks/Tissue/results/20231202/csvs/lipid.csv'
filepath.sample='D:/repositories/liver-cancer/tasks/Tissue/data/20231202_batch_sum/new_sample.csv'
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
dim(df.metab)
metab_num=ncol(df.metab)-LABEL_NUM
df.raw_metab=df.metab[,1:metab_num]


# Generate using DF



# Relabel Codex
"
df.use_sample=df.sample
df.use_sample=df.use_sample[!is.na(df.use_sample['CODEX主要亚型']),]
df.use_sample['CODEX主要亚型']
group_col='CODEX主要亚型'

target_group='P'
df.use_sample$tmp_codex=apply(df.use_sample,1,get_codex_label, group_col=group_col,target_group=target_group)
table(df.use_sample$tmp_codex)
"

# Relabel Survival
class_label='tmp_survival'

# Select Batch
df.use_sample=df.sample
# df.use_sample=df.use_sample[df.use_sample[['CODEX主要亚型']] %in% c('A'),]
df.use_sample=df.use_sample[df.use_sample[['TMN']] %in% c('II'),]
df.use_sample
low_cutoff=36
high_cutoff=60
df.use_sample$tmp_survival=apply(df.use_sample[,c('os','oss')],1,get_survival_label, low_cutoff=low_cutoff,high_cutoff=high_cutoff)
df.use_sample[(df.use_sample['oss']==1)&(df.use_sample['os']<=12),]
table(df.use_sample$tmp_survival)

# merge
df.use=merge(df.raw_metab,df.use_sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
!is.na(df.use[[class_label]])
df.use[[class_label]]
# class_label='batch'
df.use=df.use[df.use[[class_label]]!='0',]
df.use=df.use[df.use[[class_label]]!='Neg',]
pvalue_cutoff=5e-2
use_row_ha=TRUE
col_split=TRUE
draw_single_heatmap(df.use,metab_num,class_label,pvalue_cutoff,use_row_ha )
table(df.use_sample$tmp_survival)


##### Modify Test ######

LABEL_NUM=0
filepath.metab='~/workstation/MetabSubtype/tasks/MultiOmics/data/Using/metab.csv'
filepath.sample='~/workstation/MetabSubtype/tasks/MultiOmics/data/Using/sample.csv'
df.metab<-read.csv(filepath.metab, header= TRUE, check.names=F,row.names=1)
df.sample<-read.csv(filepath.sample, header= TRUE, check.names=F)
df.metab[df.metab=='']=NA
df.sample[df.sample=='']=NA
df.sample=df.sample[!is.na(df.sample['Sample Name']),c('Sample Name','主要分型','生存时间分组','MT2结构')]
df.sample['Sample Name']
row.names(df.sample)=df.sample[['Sample Name']]
df.sample=df.sample[rownames(df.metab),]
df.sample
rownames(df.metab)
rownames(df.sample)
# # Can Delete 
# df.metab=df.metab[rowSums(is.na(df.metab[,ha_cols]))==0,,drop=FALSE]
ha_cols=c("主要分型","生存时间分组")
df.use=df.metab
df.ha=df.sample[,ha_cols]
df.ha=df.ha[rowSums(is.na(df.ha[,ha_cols]))==0,,drop=FALSE]
df.ha=df.ha[df.ha[['主要分型']] %in% c('A','N'),]
df.ha
candidate_colors <- c("#8ea3c2",  "#8ccdbf", "#edb17f",  "#f1df82","#FEA3A2","#8E8BFE")
type_1_colors=candidate_colors[1:length(unique(df.ha[,1]))]
type_2_colors=candidate_colors[1:length(unique(df.ha[,2]))]
unique(df.ha[[ha_cols[2]]])
type_2_colors
type_1_colors=setNames(type_1_colors,unique(df.ha[[ha_cols[1]]]))
type_2_colors=setNames(type_2_colors,unique(df.ha[[ha_cols[2]]]))

color_list=list(type_1_colors,type_2_colors)
names(color_list)=ha_cols
color_list

# color_list[2]
top_ha_1 = HeatmapAnnotation(df=df.ha[1],col = color_list[1],
                             simple_anno_size=unit(4, "mm"),
                             #border = T,
                             gap = unit(1, "points"),
                             show_annotation_name = TRUE)
top_ha_1
top_ha_2 = HeatmapAnnotation(df=df.ha[2],col = color_list[2],
                             simple_anno_size=unit(4, "mm"),
                             #border = T,
                             gap = unit(1, "points"),
                             show_annotation_name = TRUE)
top_ha_2


class_label ='主要分型'
col_ha=df.ha[,class_label]
col_ha
length(col_ha)
metab_num=ncol(df.metab)-LABEL_NUM
metab_num
loaddata=df.use[rownames(df.ha),1:metab_num]
mat=data.matrix(loaddata)
rownames(mat)=rownames(loaddata)
mat=t(mat)

# bk <- c(seq(0,10000,by=0.01))
bk=c(-5,0,1.5)
col_fun<-colorRamp2(
  bk,
  c("#1D91C0", "white", "#E31A1C"))
ComplexHeatmap::pheatmap(mat,
                         col = col_fun,
                         name = "Relative level",
                         column_split = df.ha,
                         #row_split = row_ha,
                         top_annotation=c(top_ha_1,top_ha_2),
                         cluster_row_slices=T,
                         fontsize_row = 7, 
                         cluster_cols = T,
                         row_title_rot = 0,
                         fontsize_col = 5,
                         treeheight_row = 20,
                         row_title_gp = gpar(fontsize = 7),
                         legend_breaks=bk,
                         scale="row",
                         breaks = bk
)

class(c(top_ha_1,top_ha_2))
# ComplexHeatmap::pheatmap(mat,
#                                 #col = col_fun,
#                                 name = "Relative level",
#                                 column_split = list('MT2',
#                                 #row_split = row_ha,
#                                 top_annotation=c(top_ha_1,top_ha_2),
#                                 cluster_row_slices=F,
#                                 fontsize_row = 7, 
#                                 cluster_cols = T,
#                                 row_title_rot = 0,
#                                 fontsize_col = 5,
#                                 treeheight_row = 20,
#                                 row_title_gp = gpar(fontsize = 7),
#                                 legend_breaks=seq(-4,4,2),
#                                 scale="row",
#                                 breaks = bk
# )
# I mean I want to split columns by "MT2" and "Survival", but if I changed code to `column_split = c('MT2','Survival')`
# There is error "rror in as.data.frame.default(x[[i]], optional = TRUE) : 
#   cannot coerce class 'structure("HeatmapAnnotation", package = "ComplexHeatmap")' to a data.frame"

# # Headannotation Color Mapping
# ###
# ha_col='克隆'
# df.ha=df.metab[,ha_col,drop=FALSE]
# df.ha=df.ha[!is.na(df.ha),ha_col,drop=FALSE]
# df.ha
# ha_color_mapping=c()
# ha_color_mapping=candidate_colors[1:length(unique(df.ha[[ha_col]]))]
# ha_color_mapping= setNames(ha_color_mapping,unique(df.ha[[ha_col]]))
# color_list=NULL
# color_list[[ha_col]]=ha_color_mapping
# color_list
# top_ha=NULL
# df.ha
# top_ha = HeatmapAnnotation(df=df.ha,col = color_list,
#     simple_anno_size=unit(4, "mm"),
#     #border = T,
#     gap = unit(1, "points"),
#     show_annotation_name = TRUE)
# color_list
# top_ha
###

