##################################################################
################ Integrated Subtype Analysis Figures #############
##################################################################

###### Setup & Data Preparation ######
### Load Libraries ###
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Rtsne)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(viridis)
library(colorRamp2)

source("/home/suh/workstation/MetabSubtype/tasks/Figure/app/utils.R")

preprocess_data <- function(data_matrix) {
    # Wash
    data_matrix=data_matrix[,colnames(data_matrix)!="s-adenosyl-l-homocysteine_pos"]
    data_matrix[data_matrix=='']=NA

    # log
    data_matrix = log2(data_matrix)
    # Scaling
    raw_metab_mean=apply(data_matrix,1,mean)
    raw_metab_std=apply(data_matrix,1,sd)
    data_matrix=(data_matrix-raw_metab_mean)/raw_metab_std
    
    return(data_matrix)
}
preprocess_sample <- function(data_matrix){
    data_matrix[data_matrix=='']=NA
    data_matrix[data_matrix=='Neg']=NA
    data_matrix=data_matrix[!is.na(data_matrix['Sample Name']),]
    row.names(data_matrix)=data_matrix[['Sample Name']]
    return(data_matrix)
}
### Load Data ###
# Load your data files
filepath.metab <- "/home/suh/workstation/MetabSubtype/tasks/Figure/data/Using/metab.csv"
filepath.lipid <- "/home/suh/workstation/MetabSubtype/tasks/Figure/data/Using/lipid.csv"
filepath.sample <- "/home/suh/workstation/MetabSubtype/tasks/Figure/data/Using/whole_sample_info.csv"

df.metab <- read.csv(filepath.metab, row.names = 1, check.names = F)
df.lipid <- read.csv(filepath.lipid, row.names = 1, check.names = F)
df.sample <- read.csv(filepath.sample, row.names = 1, check.names = F)

df.metab = preprocess_data(df.metab)
df.lipid = preprocess_data(df.lipid)
df.sample = preprocess_sample(df.sample)

metab_num = ncol(df.metab)

class_label = "metab_cluster"
###### Figure 1: Metabolite-Based Subtypes ######
##### Figure 1a: Consensus Matrix #####
metab_data <- data.matrix(df.metab)
metab_cc <- ConsensusClusterPlus(metab_data,
                                 maxK = 5, 
                                 reps = 100,
                                 title = "Metab_Consensus",
                                 plot = "pdf")

##### Figure 1b: Survival Analysis #####
surv_data <- Surv(df.sample$os, df.sample$oss)
metab_clusters <- df.sample[,'metab_cluster']
surv_fit <- survfit(surv_data ~ metab_clusters)
surv_plot <- ggsurvplot(surv_fit, data = df.sample,
                        pval = TRUE, risk.table = TRUE,
                        palette = c("#E69F00", "#56B4E9"))

##### Figure 1c: Heatmap #####


df.use=merge(df.metab,df.sample[,class_label,drop=FALSE],by="row.names")
row.names(df.use)=df.use[[1]]
df.use=df.use[,2:ncol(df.use)]
df.use = df.use[!is.na(df.use[,class_label]),]
selected_metabs <- filter_by_wilcox(df.use, colnames(df.metab), class_label , 0.05)
metab_heatmap <- draw_heatmap(df.use[,c(selected_metabs,class_label)], 
                             feature_cols = selected_metabs,class_label = class_label,
                             ha_col = class_label, use_row_ha = FALSE, col_split = TRUE)

##### Figure 1d: Volcano Plot #####
de_results <- pvalue_by_wilcox(df.use, ncol(df.metab), class_label)
volcano_plot <- draw_volcano(de_results, 0.05, 1.5)

##### Figure 1e: KEGG Pathways #####
kegg_results <- kegg_after_process(
    kegg_before_process(df.metab, ncol(df.metab), "cluster", c("1","2")),
    0.05, 1.5
)
kegg_bubble <- dotplot(kegg_results, showCategory=15) + 
               theme(axis.text.y = element_text(size=8))

###### Figure 2: Lipid-Based Subtypes ######  
##### Figure 2a: Lipid Consensus Matrix #####
lipid_data <- t(scale(log2(df.lipid + 1)))
lipid_cc <- ConsensusClusterPlus(lipid_data,
                                 maxK = 5,
                                 title = "Lipid_Consensus",
                                 plot = "pdf")

##### Figure 2b: Lipid Class Distribution #####
lipid_grep <- get_grep_df(df.lipid, ncol(df.lipid))
lipid_class_plot <- lipid_grep %>%
    gather(Class, Value, -cluster) %>%
    ggplot(aes(x=Class, y=Value, fill=cluster)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle=45, hjust=1))

##### Figure 2c: Lipid Heatmap #####
lipid_sel <- filter_by_wilcox(df.lipid, colnames(df.lipid), "cluster", 0.05)
lipid_heatmap <- draw_heatmap(log2(df.lipid[lipid_sel]), 
                             lipid_sel, "cluster", use_row_ha=TRUE)

##### Figure 2d: Lipid Volcano #####
lipid_de <- pvalue_by_wilcox(df.lipid, ncol(df.lipid), "cluster")
lipid_volcano <- draw_volcano(lipid_de, 0.05, 1.5)

##### Figure 2e: Boxplots #####
top_lipids <- head(arrange(lipid_de, pvalue)$Metabolites, 6)
lipid_boxplots <- wrap_plots(lapply(top_lipids, function(x) 
    draw_boxplot(df.lipid, x, "cluster")))

###### Integrate All Figures ######
figure1 <- (surv_plot$plot | (volcano_plot / kegg_bubble)) /
           (wrap_elements(grid.grabExpr(draw(metab_heatmap))))

figure2 <- (lipid_class_plot | lipid_volcano) /
           (wrap_elements(grid.grabExpr(draw(lipid_heatmap))) | lipid_boxplots)

ggsave("Figure1.pdf", figure1, width=16, height=12)
ggsave("Figure2.pdf", figure2, width=16, height=12)
