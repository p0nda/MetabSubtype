# Install packages from CRAN

# stringr
if (!require(stringr)) {
  install.packages("stringr")
}
library(stringr)

# dplyr
if (!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)

# ComplexHeatmap
if (!require(ComplexHeatmap)) {
  install.packages("ComplexHeatmap")
}
library(ComplexHeatmap)

# ggplot2
if (!require(ggplot2)) {
  install.packages("ggplot2")
}
library(ggplot2)

# rJava
if (!require(rJava)) {
  install.packages("rJava")
}
library(rJava)

# xlsx
if (!require(xlsx)) {
  install.packages("xlsx")
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

library(clusterProfiler)

# KEGGREST
install.packages("BiocManager") 
BiocManager::install("KEGGREST")

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

