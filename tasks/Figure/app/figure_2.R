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

### Load Data ###
# Load your data files
filepath.metab <- "path/to/metab.csv"
filepath.lipid <- "path/to/lipid.csv"
filepath.sample <- "path/to/sample.csv"

df.metab <- read.csv(filepath.metab, row.names = 1, check.names = F)
df.lipid <- read.csv(filepath.lipid, row.names = 1, check.names = F)
df.sample <- read.csv(filepath.sample, row.names = 1)

###### Figure 1: Metabolite-Based Subtypes ######
##### Figure 1a: Consensus Matrix #####
metab_data <- t(scale(log2(df.metab + 1)))
metab_cc <- ConsensusClusterPlus(metab_data,
                                 maxK = 5, 
                                 reps = 100,
                                 title = "Metab_Consensus",
                                 plot = "pdf")

##### Figure 1b: Survival Analysis #####
surv_data <- Surv(df.sample$OS, df.sample$Status)
metab_clusters <- metab_cc[[2]]$consensusClass
surv_fit <- survfit(surv_data ~ metab_clusters)
surv_plot <- ggsurvplot(surv_fit, data = df.sample,
                        pval = TRUE, risk.table = TRUE,
                        palette = c("#E69F00", "#56B4E9"))

##### Figure 1c: Heatmap #####
selected_metabs <- filter_by_wilcox(df.metab, colnames(df.metab), "cluster", 0.05)
metab_heatmap <- draw_heatmap(log2(df.metab[selected_metabs]), 
                             feature_cols = selected_metabs,
                             class_label = "cluster")

##### Figure 1d: Volcano Plot #####
de_results <- pvalue_by_wilcox(df.metab, ncol(df.metab), "cluster")
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
