
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
    if(col_split==F){
      col_ha=NULL
    }
    print(col_ha)
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
    return(p)
}
