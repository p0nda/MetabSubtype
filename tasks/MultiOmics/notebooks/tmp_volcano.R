


##### Prepare Data ######
##### Metab #####
LABEL_NUM=0
filepath.metab='~/workstation/MetabSubtype/tasks/Subtype/data/Using/metab.csv'
filepath.sample='~/workstation/MetabSubtype/tasks/Subtype/data/Using/sample.csv'

filepath.cluster_result='/home/suh/workstation/MetabSubtype/tasks/MultiOmics/data/metab/cluster_result.csv'
df.cluster_result<-read.csv(filepath.cluster_result, header= TRUE, check.names=F,row.names=1)
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
dim(df.sample)
dim(merge(df.metab,df.sample,by.x="row.names",by.y="Sample Name"))


# Normalize Data
df.raw_metab.log=log2(df.raw_metab)
df.raw_metab.scaled=df.raw_metab.log
raw_metab_mean=apply(df.raw_metab.log,1,mean)
raw_metab_std=apply(df.raw_metab.log,1,sd)
df.raw_metab.scaled=(df.raw_metab.scaled-raw_metab_mean)/raw_metab_std
df.raw_metab.scaled=as.data.frame(nneg(as.matrix(df.raw_metab.scaled),method='min'))
colnames(df.raw_metab.scaled)
dim(df.raw_metab.scaled)

##### Draw ######
class_label='kmeans_2_clusters'
df.use_sample=df.sample
df.use=merge(df.raw_metab,df.cluster_result[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]

volcano.test.results=pvalue_by_wilcox(df.use,metab_num,class_label)

type1='1'
type2='2'
volcano.test.results$FC=apply(df.use[,1:metab_num], 2, 
                              function(x) 
                                mean(as.numeric(x[which(df.use[class_label] == type2)]))/
                                mean(as.numeric(x[which(df.use[class_label] == type1)])))
head(volcano.test.results)
# Log2 (FC)
volcano.test.results$log2FC <- log2(volcano.test.results[,'FC'])

dim(volcano.test.results[volcano.test.results$pvalue<5e-2,])
dim(volcano.test.results[(volcano.test.results$pvalue<5e-2)&(abs(volcano.test.results$log2FC)>log2(1.4)),])
##标记上调下调
pvalue_cutoff=5e-2
fc_up_cutoff=1.2
fc_down_cutoff=1/fc_up_cutoff

draw_volcano(volcano.test.results,pvalue_cutoff,fc_up_cutoff)