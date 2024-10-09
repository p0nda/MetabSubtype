
#__________________________
# Heatmap
#__________________________

##### Cluster #####

perform_clustering <- function(data,metab_num, nmf_k, kmeans_k) {
  df.feature <- data[,1:metab_num]  # Assuming the first column is not a feature
  dim(df.feature)
  
  for (k in c(3,2)){
    # set.seed(10403)
    
    # Perform NMF clustering
    nmf_object <- nmf(t(df.feature), rank = k, method = "brunet")
    nmf_result <- predict(nmf_object)
    nmf_clusters=as.numeric(nmf_result)
    names(nmf_clusters)=names(nmf_result)
    col_name=paste("nmf",k,"clusters",sep='_')
    data[col_name]=nmf_clusters
    # Perform k-means clustering
    # set.seed(10403)
    col_name=paste("kmeans",k,"clusters",sep='_')
    # set.seed(66)
    kmeans_clusters <- kmeans(df.feature, centers = k)$cluster
    # Add cluster assignments as new columns to the dataframe
    data[col_name]=kmeans_clusters
  }
  
  df.cluster_result=data[,(ncol(data)-3):ncol(data)]
  rownames(df.cluster_result)=rownames(data)
  return(df.cluster_result)
}

# cluster_label_correction(df.cluster_result)

cluster_label_correction<-function(df,reference_col,correction_cols){
  reference_labels=df[,reference_col]
  for (col_name in correction_cols) {
    
    # Create a table to store sample counts for each reference label
    label_counts <- table(reference_labels, df[, col_name])
    
    # Find the most similar reference cluster for each cluster in the current column
    corrected_labels <- apply(label_counts, 1, function(x) names(which.max(x)))
    
    # Replace original labels with corrected labels in the current column
    df[, col_name] <- corrected_labels[df[, col_name]]
  }
  return(df)
}

##### Lipid #####

get_grep_df <- function(lipid_df, lipid_num) {
    # Select relevant columns
    colnames <- colnames(lipid_df)[1:lipid_num]
    pure_df <- lipid_df[, colnames]
    # Group columns by lipid class
    lipid_dict <- list()
    for (col in colnames) {
        lipid_class <- str_extract(col,'([A-Z]*[a-z]*)*')  # Extract lipid class using regex
        lipid_dict[[lipid_class]] <- c(lipid_dict[[lipid_class]], col)
    }

    # Create new DataFrame with averages for each lipid class
    grep_df <- data.frame(row.names = rownames(pure_df))
    # pure_df <- as.data.frame(lapply(pure_df, as.numeric),check.names=FALSE)  # Ensure numeric columns
    for (lipid_class in names(lipid_dict)) {
        grep_df[[lipid_class]] <- rowSums(pure_df[, lipid_dict[[lipid_class]], drop = FALSE], na.rm = TRUE)
    }
    if(lipid_num!=ncol(lipid_df)){
        grep_df <- cbind(grep_df, lipid_df[, ((lipid_num+1):ncol(lipid_df))])
    }
    return(grep_df)
}

drop_odd_chain_cols <- function(df) {
    odd_chainlength_cols <- character(0)
    df=df.raw_metab
    for(col in colnames(df)){
        chain_lengths <- as.integer(regmatches(col, gregexpr("\\d+", col))[[1]])
        chain_length_index=1
        while(chain_length_index<=length(chain_lengths)){
            if(chain_length_index%%2==0){
                chain_length_index=chain_length_index+1
                next
            }
            ele=chain_lengths[chain_length_index]
            if ((ele) %% 2 != 0) {   
                odd_chainlength_cols <- c(odd_chainlength_cols, col)
            }
            chain_length_index=chain_length_index+1

        }
    }
    df <- df[, setdiff(colnames(df), odd_chainlength_cols)]

    return(df)
}
##### Metab Name Process #####

get_kegg_metabs=function(metab_list){
  search_result=c()
  target_db=df.new[['query']]
  for (search_metabolite in metab_list){
    if (length(grep(paste('\\b',search_metabolite,'\\b',sep=''),target_db,ignore.case = TRUE))==0){
      search_result=append(search_result,'NA')
    }
    else{
      search_result=append(search_result,target_db[grep(paste('\\b',search_metabolite,'\\b',sep=''),target_db,ignore.case = TRUE)[1]])
    }
  }
  return(search_result)
}
get_ms_metabs=function(metab_list){
  search_result=c()
  target_db=colnames(df.raw_metab)
  search_metabolite='inosine'
  for (search_metabolite in metab_list){
    if (length(grep(paste('\\b',search_metabolite,'\\b',sep=''),target_db,ignore.case = TRUE))==0){
      search_result=append(search_result,'NA')
    }
    else{
      search_result=append(search_result,target_db[grep(paste('\\b',search_metabolite,'\\b',sep=''),target_db,ignore.case = TRUE)[1]])
    }
    
  }
  return(search_result)
}

##### PCA #####

draw_pca <- function(df.raw, target_col, using_num, figure_title, figure_save_path) {
    df.raw <- df.raw[!is.na(df.raw[target_col]), ]
    #candidate_colors <- c("#00BFC4", "#F8766D", "#39d68f")
    candidate_colors <- c("#8ea3c2",  "#FEA3A2","#f1df82", "#8ccdbf","#39d68f","#edb17f",  "#8E8BFE")
    
    df.pure <- df.raw[, 1:using_num]
    using_colors <- candidate_colors[1:(length(unique(df.raw[, target_col])))]
    using_colors
    df.pure <- df.pure %>% replace(is.na(.), 0)
    colnames(df.pure)[1:3]
    df_pca <- prcomp(df.pure, scale. = TRUE)
    df_pcs <- data.frame(df_pca$x, target_col = df.raw[target_col],check.names = FALSE)
    dim(df_pcs)
    colnames(df_pcs)
    percentage <- round((df_pca$sdev)^2 / sum((df_pca$sdev)^2) * 100, 2)
    percentage <- paste(colnames(df_pcs), "(", paste(as.character(percentage), "%", ")", sep = ""))
    df_pcs[target_col]
    PCA_FONTSIZE <- 18
    p <- ggplot(df_pcs, aes_string(x = "PC1", y = "PC2", fill = c(target_col))) +
        ggtitle(figure_title) +
        scale_fill_manual(values = using_colors) +
        geom_point(size = 3, alpha = 1, shape = 21, color = "black") +
        xlab(percentage[1]) +
        ylab(percentage[2]) +
        stat_ellipse(level = 0.93, show.legend = F, linetype = 2) +
        theme_bw() +
        theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1) +
        theme(
            panel.border = element_rect(colour = "black", size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
        ) +
        theme(
            plot.title = element_text(size = 20, vjust = 6, hjust = 0.5, face = "plain"),
            axis.title.x = element_text(size = PCA_FONTSIZE, vjust = -2, hjust = 0.55, face = "plain"),
            axis.title.y = element_text(size = PCA_FONTSIZE, vjust = 2, face = "plain"),
            axis.text.x = element_text(size = PCA_FONTSIZE, face = "plain", colour = "black"),
            axis.text.y = element_text(size = PCA_FONTSIZE, face = "plain", colour = "black"),
            legend.title = element_text(size = PCA_FONTSIZE, face = "plain", colour = "black"),
            legend.text = element_text(size = PCA_FONTSIZE),
            legend.key.size = unit(0.4, "cm"),
            legend.spacing = unit(0.4, "cm")
        )
    if (figure_save_path != "") {
        ggsave(figure_save_path)
    } else {
        print(p)
    }
}

draw_plsda <- function(df.raw, target_col, using_num, figure_title, figure_save_path) {
    # candidate_colors <- c("#00BFC4", "#F8766D", "#39d68f")
    candidate_colors <- c("#8ea3c2",  "#FEA3A2","#f1df82", "#8ccdbf","#39d68f","#edb17f",  "#8E8BFE")
    
    df.raw <- df.raw[!is.na(df.raw[target_col]), ]
    using_colors <- candidate_colors[1:(length(unique(df.raw[, target_col])))]
    df.pure <- df.raw[, 1:using_num]
    df.pure <- df.pure %>% replace(is.na(.), 0)

    df.target <- df.pure
    df.target[target_col] <- df.raw[target_col]

    group <- df.target[[target_col]]
    # group<-df.target$stage
    df_plsda <- plsda(df.pure, group, ncomp = 4)
    df_plsda
    plsda_result_eig <- {
        df_plsda$prop_expl_var$X
    }[1:2]
    plsda_result_eig
    sample_site <- data.frame(df_plsda$variates)[1:2]
    sample_site
    sample_site$names <- rownames(sample_site)
    names(sample_site)[1:2] <- c("plsda1", "plsda2")
    sample_site$group <- group

    GROUP_FONTSIZE <- 15
    p <- ggplot(sample_site, aes(x = plsda1, y = plsda2, fill = group)) +
        ggtitle(figure_title) +
        scale_fill_manual(values = using_colors) +
        geom_point(size = 3, alpha = 1, shape = 21, color = "#030504") +
        stat_ellipse(level = 0.97, show.legend = F, linetype = 2) +
        labs(
            x = paste("PLS-DA axis1 ( ", round(100 * plsda_result_eig[1], 2), "% )", sep = ""),
            y = paste("PLS-DA axis2 ( ", round(100 * plsda_result_eig[2], 2), "% )", sep = "")
        ) +
        theme_bw() +
        theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1) +
        theme(
            panel.border = element_rect(colour = "black", size = 0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")
        ) +
        theme(
            plot.title = element_text(size = 20, vjust = 6, hjust = 0.5, face = "plain"),
            axis.title.x = element_text(size = GROUP_FONTSIZE, vjust = -2, hjust = 0.55, face = "plain"),
            axis.title.y = element_text(size = GROUP_FONTSIZE, vjust = 2, face = "plain"),
            axis.text.x = element_text(size = GROUP_FONTSIZE, face = "plain", colour = "black"),
            axis.text.y = element_text(size = GROUP_FONTSIZE, face = "plain", colour = "black"),
            legend.title = element_text(size = GROUP_FONTSIZE, face = "plain", colour = "black"),
            legend.text = element_text(size = GROUP_FONTSIZE),
            legend.key.size = unit(0.4, "cm"),
            legend.spacing = unit(0.4, "cm")
        )
    if (figure_save_path != "") {
        ggsave(figure_save_path)
    } else {
        print(p)
    }
}
##### BoxPlot #####
draw_boxplot<-function(df.use,feature,class_label){
    p=ggplot(df.use, aes_string(x = class_label, y = feature) )+
        geom_boxplot(width = 0.4,aes_string(color=class_label)) +  # Remove outlier points
        geom_jitter(width = 0.1,aes_string(color=class_label))+
        # stat_compare_means()+
        stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = 0.5,
                        label.y = max(df.use[[feature]]) + 1,size=8) +
        theme_bw() +
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
}

##### Label Process #####

get_survival_label<-function(sample_info,low_cutoff,high_cutoff){
  if (sample_info[["os"]] <= low_cutoff) {
    if (sample_info[["oss"]] == 1) {
      return("STS")
    } else {
      return(NA)  # Use NA for missing values in R
    }
  } else if (sample_info[["os"]] >= high_cutoff) {
    return("LTS")
  } else {
    return(NA)
  }
}
get_codex_label<-function(sample_info,group_col,target_group){
    if (sample_info[[group_col]] == target_group) {
        return(target_group)
    }
    else {
        return('others')
  }
}

##### Heatmap #####

filter_by_wilcox<-function(df.test,feature_cols,class_label,threshold=100){
    df.test=df.test[!is.na(df.test[class_label]),]
    df.test.results=data.frame(Metabolites = feature_cols)                        
    df.test.results$wilcox <- apply(df.test[, feature_cols], 2,
                            function(x) unlist(wilcox.test(as.numeric(x) ~ df.test[,class_label], data = df.test, exact = FALSE)[3]))
    # min(df.test.results[2])
    sig_df<-which(df.test.results$wilcox<threshold)
    a<-df.test.results[sig_df,]
    sig_list<-a$Metabolites    
    df.selected<-df.test[,c(sig_list,class_label),drop=FALSE]
    print(dim(df.selected))
    return(df.selected)      
}
pvalue_by_wilcox<-function(df.test,metab_num,class_label){
    df.test=df.test[!is.na(df.test[class_label]),]
    feature_cols=colnames(df.test)[1:metab_num]
    df.test.results=data.frame(Metabolites = feature_cols)                        
    df.test.results$pvalue <- apply(df.test[, feature_cols], 2,
                            function(x) unlist(wilcox.test(as.numeric(x) ~ df.test[,class_label], data = df.test, exact = FALSE)[3]))

    return(df.test.results)      
}
pvalue_by_kruskal<-function(df.test,metab_num,class_label){
  df.test=df.test[!is.na(df.test[class_label]),]
  feature_cols=colnames(df.test)[1:metab_num]
  df.test.results=data.frame(Metabolites = feature_cols)                        
  df.test.results$kruskal <- apply(df.test[, feature_cols], 2,
                                   function(x) unlist(kruskal.test(as.numeric(x) ~ df.test[,class_label], data = df.test)[3]))

  return(df.test.results)      
}

filter_by_ttest<-function(df.test,feature_cols,class_label,threshold=100){
    df.test=df.test[!is.na(df.test[class_label]),]
    df.test=na.omit(df.test[,c(feature_cols,class_label)])
    df.test.results=data.frame(Metabolites = feature_cols)    
    print(dim(df.test[, feature_cols]))            
    print(head(df.test[, 1]))
        
    df.test.results$ttest <- apply(df.test[, feature_cols], 2,
                            function(x) {
                                if(sum(!is.na(x)) <= 1 || sd(x, na.rm = TRUE) == 0) {
                                    return(NA)
                                } else {
                                    return(unlist(t.test(as.numeric(x) ~ df.test[,class_label], data = df.test, exact = FALSE)[3]))
                                }
                            })
    sig_df<-which(df.test.results$ttest<threshold)
    a<-df.test.results[sig_df,]
    sig_list<-a$Metabolites    
    df.selected<-df.test[,c(sig_list,class_label),drop=FALSE]
    print(dim(df.selected))
    return(df.selected)      
}

filter_by_kruskal<-function(df.test,feature_cols,class_label,threshold=100){
    df.test.results=data.frame(Metabolites = feature_cols)                        
    df.test.results$kruskal <- apply(df.test[, feature_cols], 2,
                            function(x) unlist(kruskal.test(as.numeric(x) ~ df.test[,class_label], data = df.test)[3]))
    # min(df.test.results[2])
    sig_df<-which(df.test.results$kruskal<=threshold)
    a<-df.test.results[sig_df,]
    sig_list<-a$Metabolites    
    df.selected<-df.test[,c(sig_list,class_label),drop=FALSE]
    print(dim(df.selected))
    return(df.selected)      
}

filter_by_pvalue<-function(df.test,feature_cols,class_label,threshold=100){
    label_type_num=length(unique(df.test[[class_label]]))
    
    if(label_type_num>2){
        df.selected=filter_by_kruskal(df.test,feature_cols,class_label,threshold)
        print(paste(class_label,'anova'))
    }else{
        df.selected=filter_by_wilcox(df.test,feature_cols,class_label,threshold)
        print(paste(class_label,'wilcox'))
    }
    print(df.selected[,1,drop=FALSE])
    return(df.selected)

}
draw_heatmap<-function(loaddata,feature_cols,class_label,ha_col,use_row_ha=FALSE,col_split=TRUE){
    candidate_colors <- c("#8ea3c2",  "#8ccdbf", "#edb17f",  "#f1df82","#FEA3A2","#8E8BFE")
    # candidate_colors=c("#FEA3A2","#8E8BFE","#")
    mat=data.matrix(loaddata[,feature_cols])
    rownames(mat)=rownames(loaddata)
    mat=t(mat)
    df.label_col=loaddata[,class_label]
    df.ha=loaddata[,ha_col,drop=FALSE]
    i=1
    ha_color_mapping=c()
    ha_color_mapping=candidate_colors[1:length(unique(df.ha[[ha_col]]))]
    ha_color_mapping= setNames(ha_color_mapping,unique(df.ha[[ha_col]]))
    color_list=NULL
    color_list[[ha_col]]=ha_color_mapping
    top_ha = HeatmapAnnotation(df=df.ha,col = color_list,
        simple_anno_size=unit(4, "mm"),
        #border = T,
        gap = unit(1, "points"),
        show_annotation_name = TRUE)
    col_ha=df.label_col
    break_low_boundary=-2
    break_high_boundary=3
    break_step_length=0.01
    # bk <- c(seq(break_low_boundary,0,by=0.01),seq(0,break_high_boundary,by=0.01))
    bk <- c(seq(break_low_boundary,break_high_boundary,by=break_step_length))
    bk=c(-1.5,0,1.5)
    col_fun<-colorRamp2(
        bk,
        c("#1D91C0", "white", "#E31A1C"))
    if(use_row_ha){
        row_ha=data.frame(lipid=feature_cols)
        row_ha$headgroup=str_extract(feature_cols,'([A-Z]*[a-z]*)*')
        row_ha=as.vector(row_ha$headgroup)
        p<-ComplexHeatmap::pheatmap(mat,
                                    col = col_fun,
                                    name = "Relative level",
                                    column_split = col_ha,
                                    row_split = row_ha,
                                    top_annotation=top_ha,
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
                                    breaks = bk
                                )
    }else{
        p<-ComplexHeatmap::pheatmap(mat,
                                col = col_fun,
                                name = "Relative level",
                                column_split = col_ha,
                                #row_split = row_ha,
                                top_annotation=top_ha,
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
                                breaks = bk#,
                                # column_split = col_ha
                                )
    }
    
    q=draw(
        p,  
        heatmap_legend_side = "left", 
        annotation_legend_side = "left"
    )
    return(p)
}

draw_single_heatmap<-function(df.raw,feature_num,class_label,pvalue_cutoff=0.05,use_row_ha=FALSE){
    print(class_label)
    df.raw.log=log2(df.raw[,1:feature_num])
    print(sum(is.na(df.raw.log)))

    # df.raw.log[class_label]=df.raw[,class_label]
    # df.raw.log=df.raw.log[!is.na(df.raw.log[class_label]),]
    # df.test=df.raw.log

    # df.raw[class_label]=df.raw[,class_label]
    df.raw=df.raw[!is.na(df.raw[class_label]),]
    df.test=df.raw
    print(sum(is.na(df.test)))

    feature_cols=colnames(df.test)[1:feature_num]
    # if pvalue_cutoff>=1{

    # }else {
    #     df.raw.selected=filter_by_pvalue(df.test,feature_cols ,class_label ,pvalue_cutoff)   

    # }
    df.raw.selected=filter_by_pvalue(df.test,feature_cols ,class_label ,pvalue_cutoff) 
    metab_num.selected=ncol(df.raw.selected)-1
    print(sum(is.na(df.test)))
    loaddata=df.raw.selected
    print(sum(is.na(loaddata)))

    feature_cols=colnames(df.raw.selected)[1:metab_num.selected]
    p=draw_heatmap(loaddata ,feature_cols ,class_label,class_label,use_row_ha  )
    draw(
        p,  
        heatmap_legend_side = "left", 
        annotation_legend_side = "left"
    )
    # print(q)
    print('finished')
}

##### KEGG #####

get_symbols<-function(kegg_table_row){
  gene_ids=unlist(strsplit(kegg_table_row[['geneID']],'[/]'))
  gene_ids
  keytypes(org.Hs.eg.db)
  symbols=bitr(gene_ids, #数据集
               fromType="ENTREZID", #输入为SYMBOL格式
               toType="SYMBOL",  # 转为ENTERZID格式
               OrgDb="org.Hs.eg.db") #人类 数据库
  symbols=symbols[['SYMBOL']]
  return(paste(symbols,collapse = '/'))
  
}


kegg_before_process<-function(df.test,metab_num,class_label,types=c()){
    df.test.results=data.frame(Metabolites = (colnames(df.test)[1:metab_num]))
    # df.test.results$cpd_id=df.new$cpd_id[match(df.test.results$Metabolites,df.new$query)]
    for(metabolite in df.test.results$Metabolites){
        search_metabolite=str_replace_all(metabolite,'_pos','')
        search_metabolite=str_replace_all(search_metabolite,'_neg','')
        # search_metabolite=str_replace_all(search_metabolite,'+','')
        print(grep(metabolite,df.new[['query']],ignore.case = TRUE))
        if (length(grep(search_metabolite,df.new[['query']],ignore.case = TRUE))==0){
                df.test.results[df.test.results$Metabolites==metabolite,'cpd_id']='NA'
        }
        else{
                print(grep(search_metabolite,df.new[['query']],ignore.case = TRUE))
                df.test.results[df.test.results$Metabolites==metabolite,'cpd_id']=df.new$cpd_id[grep(search_metabolite,df.new[['query']],ignore.case = TRUE)][1]
        }
    }
    df.test.results$wilcox=apply(df.test[, 1:metab_num], 2,
                            function(x) unlist(wilcox.test(as.numeric(x) ~ df.test[[class_label]], data = df.test, exact = FALSE)[3]))
    df.test.results$wilcox_BH=p.adjust(df.test.results$wilcox,method="BH")
    df.test.results$FC <- apply(df.test[,1:metab_num], 2, 
                                function(x) 
                                mean(as.numeric(x[which(df.test[class_label] == types[1])]))/
                                mean(as.numeric(x[which(df.test[class_label] == types[2])])))
    df.test.results$log2FC<- log2(df.test.results$FC)
    return(df.test.results)
}
path_avelog2FC <- function(x,df.test.results){
    metabs <- unlist(strsplit(x[[8]], "[/]"))
    mean(abs(df.test.results[match(metabs, df.test.results$cpd_id),'FC']))
}
kegg_after_process<-function(df.test.results,pvalue_cutoff=5e-2,fc_up_cutoff=3/2,fc_down_cutoff=2/3){
    path1 <- df.path.metab %>% dplyr::select(pathway_id, cpd_id) 
    path2 <- df.path.metab %>% dplyr::select(pathway_id, pathway_name)
    # 3 pairs of kegg
    # fc_up_cutoff=3/2
    # fc_down_cutoff=2/3
    df.test.results$is_use=(df.test.results$wilcox<=pvalue_cutoff)&((df.test.results$FC>=fc_up_cutoff)|(df.test.results$FC<=fc_down_cutoff))
    # df.test.results$is_use=(df.test.results$wilcox_BH<=pvalue_cutoff)&(abs(df.test.results$log2FC)>=log2fc_cutoff)
    used_cpds=df.test.results$cpd_id[df.test.results$is_use]
    print(paste0('used metab number is ',length(used_cpds)))
    print(used_cpds)
    x <- enricher(used_cpds,
                TERM2GENE = path1, 
                TERM2NAME = path2,
                minGSSize = 1, pvalueCutoff = 1, pAdjustMethod = "BH")
    x
    kegg_table <- na.omit(as.data.frame(x))
    kegg_table$Description=str_replace(kegg_table$Description,' - Homo sapiens \\(human\\)','')
    # write.table(kegg_table,'D:/repositories/liver-cancer/tasks/Tissue/results/tmp_R/kegg_table.csv',sep=',',row.names = FALSE)
    kegg_table$Description
    # kegg_table <- kegg_table[kegg_table$Count > 2,]
    if (nrow(kegg_table)==0){
        print('Nothing is enriched.')
        return(NULL)
    }
    kegg_table$metabs_mean <- apply(kegg_table, 1, path_avelog2FC,df.test.results=df.test.results)
    # calculate fold enrichment
    kegg_table$FoldEnrich <- apply(kegg_table, 1, 
                                function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                                    as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                                    as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                                    as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
    kegg_table$query=apply(kegg_table,1,get_cpd_names)
    return(kegg_table)

}
get_cpd_names <- function(x,df.test.results){
    metabs <- unlist(strsplit(x[[8]], "[/]"))
    queries=df.new[match(metabs,df.new$cpd_id),'query']
    # print(queries)
    # print(paste(queries,sep = '/',collapse = '/'))
    return(paste(queries,collapse = '/'))
}
draw_comporison_kegg<-function(kegg_table_draw,pvalue_form){
  size_range=c(6,15)
  size_breaks=c(seq(5,25,5))
  p=ggplot(kegg_table_draw, aes(comporison, Description)) +
    geom_point(aes(fill = -log10(!!sym(pvalue_form)), size = FoldEnrich), shape = 21,stroke=1,alpha=0.9) +
    # scale_y_discrete(limits = kegg_table_draw[['Description']])+
    scale_size(range = size_range, breaks = size_breaks) +
    # scale_color_gradient(low="blue",high = "red")+
    scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4 ,breaks = -log10(c(0.05,0.01,0.005,0.001)),labels=c(0.05,0.01,0.005,0.001) )+
    #  scale_fill_viridis_c(option = "A", direction = -1, begin = 0.4, breaks = c(0.2,0.5, 1, 1.5)) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    # coord_fixed(ratio = 0.4) +
    labs(x = expression("-Log"[10]*"(P-value)"), y = "", title = "", 
         fill = expression("P-value"), size = "Enrich Ratio") +
    theme(axis.text.x = element_text(size = 20, face = "plain", colour = "black",angle = 45,hjust = 1), 
          axis.text.y = element_text(size = 18, face = "plain", colour = "black")) +
    theme(legend.text = element_text(size = 10, face = "plain", colour = "black"), 
          legend.title = element_text(size = 10, face = "plain", colour = "black"),
          legend.key.height = unit(0.3, "cm"), legend.key.width = unit(0.3, "cm"))
  print(p)

}
##### CPD Change Line #####

"filter_by_wilcox<-function(df.test,feature_cols,class_label,threshold=100){
    df.test=df.test[!is.na(df.test[class_label]),]
    df.test.results=data.frame(Metabolites = feature_cols)                        
    df.test.results$wilcox <- apply(df.test[, feature_cols], 2,
                            function(x) unlist(wilcox.test(as.numeric(x) ~ df.test[,class_label], data = df.test, exact = FALSE)[3]))
    # min(df.test.results[2])
    sig_df<-which(df.test.results$wilcox<threshold)
    a<-df.test.results[sig_df,]
    sig_list<-a$Metabolites    
    df.selected<-df.test[,c(sig_list,class_label)]
    print(dim(df.selected))
    return(df.selected)      
}

filter_by_ttest<-function(df.test,feature_cols,class_label,threshold=100){
    df.test=df.test[!is.na(df.test[class_label]),]
    df.test=na.omit(df.test[,c(feature_cols,class_label)])
    df.test.results=data.frame(Metabolites = feature_cols)    
    print(dim(df.test[, feature_cols]))            
    print(head(df.test[, 1]))
        
    df.test.results$ttest <- apply(df.test[, feature_cols], 2,
                            function(x) {
                                if(sum(!is.na(x)) <= 1 || sd(x, na.rm = TRUE) == 0) {
                                    return(NA)
                                } else {
                                    return(unlist(t.test(as.numeric(x) ~ df.test[,class_label], data = df.test, exact = FALSE)[3]))
                                }
                            })
    sig_df<-which(df.test.results$ttest<threshold)
    a<-df.test.results[sig_df,]
    sig_list<-a$Metabolites    
    df.selected<-df.test[,c(sig_list,class_label)]
    print(dim(df.selected))
    return(df.selected)      
}
filter_by_kruskal<-function(df.test,feature_cols,class_label,threshold=100){
    df.test.results=data.frame(Metabolites = feature_cols)                        
    df.test.results$kruskal <- apply(df.test[, feature_cols], 2,
                            function(x) unlist(kruskal.test(as.numeric(x) ~ df.test[,class_label], data = df.test)[3]))
    # min(df.test.results[2])
    sig_df<-which(df.test.results$kruskal<threshold)
    a<-df.test.results[sig_df,]
    sig_list<-a$Metabolites    
    df.selected<-df.test[,c(sig_list,class_label)]
    print(dim(df.selected))
    return(df.selected)      
}
filter_by_pvalue<-function(df.test,feature_cols,class_label,threshold=100){
    label_type_num=length(unique(df.test[[class_label]]))
    
    if(label_type_num>2){
        df.selected=filter_by_kruskal(df.test,feature_cols,class_label,threshold)
        print(paste(class_label,'anova'))
    }else{
        df.selected=filter_by_wilcox(df.test,feature_cols,class_label,threshold)
        print(paste(class_label,'wilcox'))
    }
    return(df.selected)

}"


combine_and_deduplicate_queries <- function(df.cpd_kegg) {
  
  df.combined_queries <- df.cpd_kegg %>%
    group_by(Description) %>%
    summarise(
      query = paste(unique(unlist(strsplit(query,'/'))), collapse = "/")
    )
  df.combined_queries=as.data.frame(df.combined_queries)
  return(df.combined_queries)
}

smartConvert <-  function(x) {
  if(any(is.na(as.numeric(as.character(x))))) x else as.numeric(x)
    # return(as.numeric(x))
}

get_sem=function(x){
    sd(x)/sqrt(length(x))
}

# MetaboAnalyst Data Process
meta_result_process<-function(raw_result_df){
  processed_df=raw_result_df
  processed_df$FoldEnrich=processed_df$hits/processed_df$expected
  processed_df$Description=str_replace(row.names(processed_df),'metabolism','')
  processed_df['p.adjust']=processed_df['Holm p']
  processed_df['p_value']=processed_df['Raw p']
  return(processed_df)
}


##### Volcano #####

draw_volcano<-function(volcano.test.results,pvalue_cutoff,fc_up_cutoff,use_std=FALSE){
    fc_down_cutoff=1/fc_up_cutoff
        
    volcano.test.results <- volcano.test.results %>% 
    mutate(Expression = case_when(FC >= fc_up_cutoff & pvalue <= pvalue_cutoff ~ "Up",
                                    FC <= fc_down_cutoff & pvalue <= pvalue_cutoff ~ "Down",
                                    TRUE ~ "Not significant"))

    ##加上标记
    if(use_std){
      volcano.test.results $Metabolites=volcano.test.results$std_Metabolites
    }
    volcano.test.results$relabel <- NA
    volcano.test.results$relabel[volcano.test.results$Expression == "Up"] <- volcano.test.results $Metabolites[volcano.test.results$Expression == "Up"]
    volcano.test.results$relabel[volcano.test.results$Expression == "Down"] <- volcano.test.results $Metabolites[volcano.test.results$Expression == "Down"]
    head(volcano.test.results)
    # write.table(volcano.test.results, file="lipid_tissue_M.csv", sep=',', col.names=T)
    
    ##画图
    xlim_cutoff=max(abs(log2(min(volcano.test.results$FC))),abs(log2(max(volcano.test.results$FC))))
    p1=ggplot(volcano.test.results, aes(log2FC, -log(pvalue,10), fill = Expression, label=relabel)) +  
      # ggplot(volcano.test.results, aes(logFC, -log(pvalue,10), fill = Expression)) +  
      geom_point(colour = "gray" ,alpha=0.8, size = 2, shape = 21)   +
      xlim((-xlim_cutoff-0.5),(xlim_cutoff+0.5))+
      geom_text_repel (size = 2.5)   +
      scale_fill_manual(values = c("Up" = alpha("#F8766D", 1),  # Full opacity for "Up"
                                   "Down" = alpha("#009000", 1),  # Full opacity for "Down"
                                   "Not significant" = alpha("grey", 0.4))) +  # Semi-transparent for "Not significant"  xlab(expression("log"[2]*"FC")) + 
      xlab(expression("log"[2]*"FC")) +
      ylab(expression("-log"[10]*"P value")) +
      theme_bw() +
      ##固定长宽比例不随着改变而改变
      theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1, panel.grid = element_blank())      +
      ##加垂直线
      geom_vline(xintercept=c(log2(fc_down_cutoff), log2(fc_up_cutoff)), col="grey", linetype="dashed") +
      geom_hline(yintercept=-log10(pvalue_cutoff), col="grey", linetype="dashed")   +
      ##固定字体大小  
      theme(axis.title.x = element_text(size = 17, hjust = 0.5, face = "plain"),
            axis.title.y = element_text(size = 17, face = "plain"), 
            axis.text.y = element_text(size = 15, face = "plain", colour = "black"),
            axis.text.x = element_text(size = 15, face = "plain", colour = "black")) +
      ##调整legend字体  
      theme(legend.text = element_text(size = 15), 
            legend.title = element_text(size = 15), 
            legend.spacing.y = unit(0.2, 'cm'), 
            legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm')) +
      theme(panel.border = element_rect(size=2))
    print(p1)
}
