# Install packages from CRAN

# stringr
if (!require(stringr)) {
  install.packages("stringr")
}
install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
install.packages("stringr", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

library(stringr)

# dplyr
if (!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# ComplexHeatmap
if (!require(ComplexHeatmap)) {
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)

# ggplot2
if (!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)


install.packages("survminer")
install.packages("survival")
install.packages("umap")
install.packages("tsne")


if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/survminer", build_vignettes = FALSE)

# rJava
if (!require(rJava)) {
  install.packages("rJava")
}
library(rJava)

# xlsx
if (!require(xlsx)) {
  install.packages("xlsx",dependencies = TRUE)
}
library(xlsx)

# RColorBrewer
if (!require(RColorBrewer)) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)

# circlize
if (!require(circlize)) {
  install.packages("circlize")
}
library(circlize)

# umap
if (!require(umap)) {
  install.packages("umap")
}
library(umap)

# Rtsne
if (!require(Rtsne)) {
  install.packages("Rtsne")
}
library(Rtsne)

# ggrepel
if (!require(ggrepel)) {
  install.packages("ggrepel")
}
library(ggrepel)

# clusterProfiler
BiocManager::install("clusterProfiler")
BiocManager::install("KEGGREST")
BiocManager::install("survminer",dependencies = TRUE, INSTALL_opts = '--no-lock')
library(KEGGREST)
BiocManager::install("clusterProfiler", dependencies = TRUE, INSTALL_opts = '--no-lock')
library(clusterProfiler)
BiocManager::install("pathview")
BiocManager::install("DESeq2")
BiocManager::install("KEGGREST")
BiocManager::install("ComBat")
BiocManager::install("sva")


BiocManager::install("org.Hs.eg.db")
library('org.Hs.eg.db')
## Test 
library('org.Hs.eg.db')
columns(org.Hs.eg.db)
ids=colnames(df.metab)[1:10000]
ids
org.Hs.eg.db$ENSEMBL
gene=mapIds(org.Hs.eg.db,keys=ids,keytype="ENSEMBL",column="ENTREZID")
kegg_result=enrichKEGG(
  gene,
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
##
# KEGGREST
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# ConsensusClusterPlus
if (!require(ConsensusClusterPlus)) {
  install.packages("ConsensusClusterPlus")
}
library(ConsensusClusterPlus)

# viridis
if (!require(viridis)) {
  install.packages("viridis")
}
library(viridis)

# survminer
if (!require(survminer)) {
  install.packages("survminer")
}
library(survminer)

# survival
if (!require(survival)) {
  install.packages("survival")
}
library(survival)

# mixOmics
BiocManager::install("mixOmics")


# patchwork
if (!require(patchwork)) {
  install.packages("patchwork")
}
library(patchwork)


##### Notebooks #####
install.packages("IRkernel")
IRkernel::installsepc()

install.packages("rmarkdown") # 写rmarkdown必需的包
install.packages("knitr") # 导出文件必需的包
install.packages("tinytex") # TeX的轻量级发行版，用于PDF文件的导出
install.packages("rticles") # 配合中文导出PDF，有很多不同的文档模板可供使用

# Bioconductor packages (replace with the specific ones you need)
# Example:
# if (!require(org.Hs.eg.db)) {
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("org.Hs.eg.db")
# }
# library(org.Hs.eg.db)

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
library(mixOmics)
library(patchwork)

