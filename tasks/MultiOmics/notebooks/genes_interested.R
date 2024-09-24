##### Prepare Data


##### Metab #####
# filepath.rna_counts='D:/repositories/liver-cancer/tasks/Tissue/results/20231202/csvs/lipid.csv'
filepath.raw_counts='/home/suh/workstation/MetabSubtype/tasks/MultiOmics/data/rna/counts.csv'
filepath.sample='/home/suh/workstation/MetabSubtype/tasks/MultiOmics/data/Using/sample.csv'
filepath.cluster_result='/home/suh/workstation/MetabSubtype/tasks/MultiOmics/data/metab/cluster_result.csv'
df.raw_counts<-read.csv(filepath.raw_counts, header= TRUE, check.names=F,row.names=1)
df.sample<-read.csv(filepath.sample, header= TRUE, check.names=F)
df.cluster_result<-read.csv(filepath.cluster_result, header= TRUE, check.names=F,row.names=1)
df.raw_counts[1:5,1:5]
dim(df.raw_counts)
length(unique((df.raw_counts[,1])))
unique(df.raw_counts[,1])
df.raw_counts[,1]
df.raw_counts[df.raw_counts=='']=NA

# Sample DF Process
df.sample[df.sample=='']=NA
df.sample[df.sample=='Neg']=NA
df.sample=df.sample[!is.na(df.sample['Sample Name']),]
df.sample['Sample Name']
row.names(df.sample)=df.sample[['Sample Name']]
df.sample[match('1520',rownames(df.sample)),'oss']=0
# df.sample=df.sample[!is.na(df.sample['Sample Name']),c('Sample Name','主要分型','生存时间分组','MT2结构')]
dim(df.sample)
# df.sample=df.sample[rownames(df.raw_counts),]
# df.sample=na.omit(df.sample,'all')
dim(df.cluster_result)
df.raw_counts[,ncol(df.raw_counts),drop=FALSE]
df.sample_cluster=merge(df.sample,df.cluster_result[,4,drop=FALSE],by.x="Sample Name",by.y="row.names",all=TRUE)

##### DGE Process
df.using=merge(df.raw_counts,df.cluster_result,by="row.names")
rownames(df.using)=df.using[,1]
df.using=df.using[,2:ncol(df.using)]
x.raw_counts=df.using[,1:ncol(df.raw_counts)]
x.cluster_result=df.using[,(ncol(df.raw_counts)+1):ncol(df.using)]
for (col in colnames(x.cluster_result)){
  x.cluster_result[col]=as.factor(x.cluster_result[[col]])
}
str(x.cluster_result)

class_label='kmeans_2_clusters'
group=factor(x.cluster_result[[class_label]])
batch=factor(x.cluster_result[['batch']])
length(group)