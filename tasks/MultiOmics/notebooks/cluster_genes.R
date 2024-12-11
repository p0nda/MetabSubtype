# Load necessary library
library(stats)
library(ggplot2)
library(Rtsne)
library(Hmisc)

# Function to perform hierarchical clustering on RNA expression data
perform_hierarchy_clustering <- function(rna_df, num_clusters = 3, method = "average") {
  
  rna_df = t(rna_df)
  # Calculate the distance matrix between genes based on expression values
  dist_matrix <- dist(rna_df, method = "euclidean")
  
  # Perform hierarchical clustering using the specified method
  hc <- hclust(dist_matrix, method = method)
  
  # Cut the dendrogram to create the specified number of clusters
  gene_clusters <- cutree(hc, k = num_clusters)
  
  # Create a dataframe to store the gene IDs and their cluster assignments
  cluster_df <- data.frame(
    GeneID = rownames(rna_df),
    Cluster = paste0("cluster", gene_clusters),
    stringsAsFactors = FALSE
  )
  
  # Return the clustering result as a dataframe
  return(cluster_df)
}

##### Prepare Data #####
filepath.edger='/home/suh/workstation/MetabSubtype/tasks/MultiOmics/data/rna/edger.csv'
df.edger<-read.csv(filepath.edger, header= TRUE, check.names=F,row.names=1)

num_clusters = 10

df.using=merge(df.raw_counts,df.sample[,c("batch",'metab_cluster','lipid_cluster')],by="row.names")
rownames(df.using)=df.using[,1]
df.using=df.using[,2:ncol(df.using)]
x.raw_counts=df.using[,1:ncol(df.raw_counts)]
x.cluster_result=df.using[,(ncol(df.raw_counts)+1):ncol(df.using)]
for (col in colnames(x.cluster_result)){
  x.cluster_result[col]=as.factor(x.cluster_result[[col]])
}

class_label='lipid_cluster'
group=factor(x.cluster_result[[class_label]])
batch=factor(x.cluster_result[['batch']])

# DGE Building
dge <- DGEList(counts=t(df.using[,1:(ncol(df.using)-3)]), group=group)

dge <- calcNormFactors(dge)

design <- model.matrix(~group)

# Build voom
v <- voom(dge, design)
dim(design)
dim(dge$counts)
# 使用 removeBatchEffect 去除批次效应
v$E <- removeBatchEffect(v$E, batch=batch, design=design)
dge <- estimateDisp(dge, design)

v_corrected.edger=as.data.frame(t(v$E))

##### Best Cluster #####
dim(v$E)
rcc = ConsensusClusterPlus(as.matrix(v$E[df.test_results[df.test_results$CODING,'ENSEMBL'],]),
                           maxK = 15, 
                           reps = 100, 
                           pItem = 0.8, 
                           pFeature = 1,
                           title = 'RNA',
                           clusterAlg = "km",
                           plot="pdf")
##### Method 1: Heatmap #####
mat = df.edger
p = Heatmap(as.matrix(t(mat)) , width = unit(10, "cm"),
            km = num_clusters, 
            cluster_columns = F, 
            show_row_names = F, 
            row_title_rot = 0, 
            row_gap = unit(3, "mm"), 
            name = "Log2FC", 
            column_title = "Test Clusters",
            column_title_gp = gpar(fontfamily = "sans", fontsize = 28), 
            column_names_gp = gpar(fontfamily = "sans", fontsize = 20),
            # col = mat[,class_label]
          )
p = draw(p)
rcl.list <- row_order(p)
table(rcl.list)
clu_df.heatmap <- do.call(rbind, lapply(1:length(rcl.list), function(i){
  out <- data.frame(GeneID = colnames(mat[rcl.list[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}))
dim(clu_df.heatmap)
table(clu_df.heatmap$Cluster)

# write.csv(clu_df.heatmap,'~/workstation/MetabSubtype/tasks/MultiOmics/results/20241114_cluster_genes/heatmap_10.csv')

##### Method 2: Kmeans Clustering #####

clu_df.km <- kmeans(t(mat), centers = num_clusters)$cluster
clu_df.km <- data.frame(
  GeneID = names(clu_df.km),
  Cluster = paste0("cluster", clu_df.km),
  stringsAsFactors = FALSE
)
table(clu_df.km$Cluster)

##### Method 3: Hierarchy #####
clu_df.hierachy = perform_hierarchy_clustering(mat, num_clusters)

table(clu_df.hierachy$Cluster)

##### Visualisation ######
mat = mat %>% distinct()
dim(mat)
set.seed(42)
tsne_mat = t(mat)
tsne_mat = as.data.frame(tsne_mat) %>% distinct()

tsne_input <- as.matrix(tsne_mat)
tsne_result <- Rtsne(tsne_input, dims = 2, perplexity = 30, verbose = TRUE)
tsne_data <- as.data.frame(tsne_result$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")
tsne_data$GeneID <- rownames(tsne_mat)

tsne_data$KM = clu_df.km[match(tsne_data$GeneID,clu_df.km$GeneID),'Cluster']
tsne_data$Heatmap = clu_df.heatmap[match(tsne_data$GeneID,clu_df.heatmap$GeneID),'Cluster']

p1 = ggplot(tsne_data, aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = KM), size = 0.1) +
  labs(
    title = "t-SNE of KM",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme_minimal() +
  scale_color_manual(values = rainbow(length(unique(tsne_data$KM))))


p2 = ggplot(tsne_data, aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = Heatmap), size = 0.1) +
  labs(
    title = "t-SNE of Heatmap",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2"
  ) +
  theme_minimal() +
  scale_color_manual(values = rainbow(length(unique(tsne_data$Heatmap))))

combined_plot = p1 + p2
print(combined_plot)


##### Correlation Clustering #####
# Step 1: Build the correlation matrix
build_correlation_matrix <- function(lipid_data, rna_data) {
  # Ensure the samples match between the two datasets
  common_samples <- intersect(rownames(lipid_data), rownames(rna_data))
  
  # Subset the data based on common samples
  lipid_data <- lipid_data[common_samples, , drop = FALSE]
  rna_data <- rna_data[common_samples, , drop = FALSE]
  
  # Compute correlation matrix (lipids as rows, RNAs as columns)
  cor_matrix <- cor(t(lipid_data), t(rna_data), use = "pairwise.complete.obs", method = "pearson")
  return(cor_matrix)
}

dim(df.edger)
dim(df.lipid)

filepath.lipid='~/workstation/MetabSubtype/tasks/Subtype/data/Using/lipid.csv'
df.lipid<-read.csv(filepath.lipid, header= TRUE, check.names=F,row.names=1)
dim(df.lipid)
length(unique((df.lipid[,1])))
df.lipid[df.lipid=='']=NA
dim(df.lipid)

dim(df.lipid)
lipid_num=ncol(df.lipid)
df.raw_lipid=df.lipid[,1:lipid_num]

df.raw_lipid=drop_odd_chain_cols(df.raw_lipid)
dim(df.raw_lipid)
lipid_num=dim(df.raw_lipid)[2]
# Normalize Data
df.raw_lipid.log=log2(df.raw_lipid)
df.raw_lipid.scaled=df.raw_lipid.log
raw_lipid_mean=apply(df.raw_lipid.log,1,mean)
raw_lipid_std=apply(df.raw_lipid.log,1,sd)
df.raw_lipid.scaled=(df.raw_lipid.scaled-raw_lipid_mean)/raw_lipid_std
# df.raw_lipid.scaled=scale(df.raw_lipid.log, center = TRUE, scale = TRUE)

# df.raw_lipid.scaled=as.data.frame(nneg(as.matrix(df.raw_lipid.scaled),method='min'))

df.cor = build_correlation_matrix(df.raw_lipid.scaled, df.edger)


lipid_data = df.raw_lipid.scaled
rna_data = df.edger
common_samples <- intersect(rownames(lipid_data), rownames(rna_data))

# Subset the data based on common samples
lipid_data <- lipid_data[common_samples, , drop = FALSE]
rna_data <- rna_data[common_samples, , drop = FALSE]

# Compute correlation matrix (lipids as rows, RNAs as columns)
dim(lipid_data)
dim(rna_data)

# Compute correlations lipid by lipid

n_cores <- detectCores()

cor_matrix <- mclapply(seq_len(ncol(lipid_data)), function(i) {
  apply(rna_data, 2, function(rna) cor(lipid_data[, i], rna, use = "pairwise.complete.obs"))
}, mc.cores = n_cores)

length(cor_matrix)
cor_matrix <- do.call(rbind, cor_matrix)
cor_matrix = t(cor_matrix)
dim(cor_matrix)

cluster_and_visualize(cor_matrix, 5)

k = 5
# Perform k-means clustering on the correlation matrix
cor_matrix_unique <- cor_matrix[!duplicated(cor_matrix), ]
kmeans_result <- kmeans(cor_matrix_unique, centers = k, nstart = 25)

cor_df = as.data.frame(cor_matrix_unique)
dim(cor_df)
colnames(cor_matrix_unique)
# Perform t-SNE dimensionality reduction
tsne_result <- Rtsne(cor_matrix_unique, dims = 2, perplexity = 30, verbose = FALSE)

# Prepare a data frame for visualization
tsne_df <- data.frame(
  X = tsne_result$Y[, 1],
  Y = tsne_result$Y[, 2],
  Cluster = factor(kmeans_result$cluster),
  RNA = rownames(cor_matrix_unique)
)

# Draw the t-SNE figure using ggplot2
tsne_plot <- ggplot(tsne_df, aes(x = X, y = Y, color = Cluster, label = RNA)) +
  geom_point(size = 1, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE of RNA Clusters", x = "t-SNE 1", y = "t-SNE 2") +
  theme(legend.title = element_blank())
tsne_plot


# Heatmap
candidate_colors <- c(
  "#8ea3c2",  # Original
  "#8ccdbf",  # Original
  "#edb17f",  # Original
  "#f1df82",  # Original
  "#FEA3A2",  # Original
  "#8E8BFE",  # Original
  "#97b7d2",  # New: lighter shade of blue
  "#7eb8af",  # New: variation of teal
  "#e7a169",  # New: darker orange
  "#f6e094",  # New: lighter yellow
  "#f88f8e",  # New: softer red
  "#7a76e6",  # New: deeper purple
  "#adc0d8"   # New: very light blue-gray
)
# candidate_colors=c("#FEA3A2","#8E8BFE","#")
mat = cor_df

col_annotation <- data.frame(Cluster = factor(kmeans_result$cluster))
rownames(col_annotation) <- colnames(cor_matrix_unique)
# Map cluster colors
col_color_mapping <- candidate_colors[1:length(unique(col_annotation$Cluster))]
col_color_mapping <- setNames(col_color_mapping, unique(col_annotation$Cluster))
col_color_list <- list(Cluster = col_color_mapping)

# Create HeatmapAnnotation for columns
top_col_ha <- HeatmapAnnotation(df = col_annotation,
                                col = col_color_list,
                                simple_anno_size = unit(4, "mm"),
                                gap = unit(1, "points"),
                                show_annotation_name = TRUE)



lipid_classes <- str_extract(colnames(df.raw_lipid.scaled), '^[A-Za-z]+')

# Annotate lipid classes with "others" for non-target groups
target_headgroups <- c('Cer', 'AcCa', 'LPC', 'LPE', 'SM', 'PC', 'PE', 'PG', 'PI', 'PS', 'TAG', 'DAG', 'TG')
annotated_classes <- ifelse(lipid_classes %in% target_headgroups, lipid_classes, "others")

# Create the row annotation
row_annotation <- data.frame(Headgroup = annotated_classes)

# Map headgroup colors
row_color_mapping <- candidate_colors[1:length(unique(row_annotation$Headgroup))]
row_color_mapping <- setNames(row_color_mapping, unique(row_annotation$Headgroup))
row_color_list <- list(Headgroup = row_color_mapping)

# Create HeatmapAnnotation for rows
row_ha <- rowAnnotation(df = row_annotation,
                        col = row_color_list,
                        simple_anno_size = unit(4, "mm"),
                        gap = unit(1, "points"),
                        show_annotation_name = TRUE)




break_low_boundary=-2
break_high_boundary=3
break_step_length=0.01
# bk <- c(seq(break_low_boundary,0,by=0.01),seq(0,break_high_boundary,by=0.01))
bk <- c(seq(break_low_boundary,break_high_boundary,by=break_step_length))
bk=c(-1.5,0,1.5)
col_fun<-colorRamp2(
  bk,
  c("#1D91C0", "white", "#E31A1C"))
if(col_split==F){
  col_ha=NULL
}
  
dim(mat)
row_ha
top_col_ha

col_split = as.vector(col_annotation$Cluster)
row_split = as.vector(row_annotation$Headgroup)

p<-ComplexHeatmap::pheatmap(t(mat),
                              col = col_fun,
                              name = "Relative level",
                              column_split = col_split,
                              row_split = row_split,
                              top_annotation=top_col_ha,
                              left_annotation = row_ha,
                              cluster_row_slices=F,
                              fontsize_row = 7, 
                              cluster_cols = T,
                              row_title_rot = 0,
                              fontsize_col = 5,
                              treeheight_row = 20,
                              treeheight_col = 20,
                              row_title_gp = gpar(fontsize = 7),
                              legend_breaks=seq(break_low_boundary,break_high_boundary,(break_high_boundary-break_low_boundary)/4),
                              scale="row",
                              breaks = bk,
                            # show_row_names = FALSE,  # Hide row names
                            # show_column_names = FALSE,  # Hide column names
  )

all_na_columns <- names(t(mat))[colSums(is.na(t(mat))) == nrow(t(mat))]

table(col_split)
dim(mat)
10377+9727+10158+16630+10540
pdf('~/workstation/MetabSubtype/tasks/MultiOmics/results/20240930_essay_figures/lipid/cor_heatmap.pdf')
draw(
  p,  
  heatmap_legend_side = "left", 
  annotation_legend_side = "left"
)
dev.off()

selected_rows <- t(mat)[rowSums(abs(t(mat)) > 0.5) == ncol(t(mat)),]
length(selected_rows)
