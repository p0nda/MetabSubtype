# Load necessary library
library(stats)
library(ggplot2)
library(Rtsne)
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
