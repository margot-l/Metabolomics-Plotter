# Margot Lautens
# Date created: 20170926
# Last edited: 20170928
# Code for plotting heat maps in R
# Input is in the format of a tidy (longform) dataset of metabolite and sample
# Adapted from code by Soumaya Zlitni and Julia Hanchard

require("dplyr")
require("tidyr")
require("pheatmap")
require("grid")
require("purrr")
require("rlang")
require("cowplot")
require("ggplotify")
require("svglite")
require("dendsort")
require("ggplot2")

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

#Calculate Fold Change

makeHeatmapTable <- function(met_data,
                             met_table=c(),
                             value_y="foldChange",
                             value_x="featureName",
                             value_z="normMethod",
                             model,
                             norm_list=c(),
                             met_list=c(),
                             pathway_name=c(),
                             data_date=format(Sys.Date(),"%Y%m%d"),
                             con_names=c(),
                             data_sort=TRUE,
                             data_mean=FALSE,
                             sym_breaks=TRUE,
                             filetype=".png",
                             same_legend=F){
  
  if (is.null(met_list)&is.null(pathway_name)){
    met_list=met_data%>%
      distinct(featureName,KEGG_ID)
  }else if (!(is.null(pathway_name))){
    met_list<- met_data%>%
      inner_join(model%>%
                   filter(Pathway==pathway_name),by="KEGG_ID")%>%
      distinct(featureName)
  }else if(!(is.data.frame(met_list))){
    met_list=as.data.frame(met_list)%>%
      rename(featureName="met_list")
  }else{
    met_list=met_list%>%
      filter(KEGG_ID!="Missing")
  }
  
  if(!(is.null(norm_list))){
    norm_list=as.data.frame(norm_list)%>%
      rename(normMethod="norm_list")%>%
      rename()
  }else{
    norm_list=met_data%>%
      distinct(normMethod)
  }
  
  filteredPlots <- met_data%>%
    semi_join(met_list%>%
                drop_na())%>%
    semi_join(norm_list%>%
                drop_na())

  
  if(data_mean){
    data_type <- "_mean"
  }else{
    data_type <- NA
  }
  if (data_sort){
    sort_type <- "_sorted"
  }else{
    sort_type <- "_unsorted"
  }
  
  max_val <- 0
  min_val <- 0
  
  if((same_legend==T)&(sym_breaks==T)){
    max_val <- max(abs(filteredPlots%>%
                         pull(.data[[value_y]])))
    min_val <- -max(abs(filteredPlots%>%
                          pull(.data[[value_y]])))
  }else if((same_legend==T)&(sym_breaks==F)){
    max_val <- max(filteredPlots%>%
                     pull(.data[[value_y]]))
    min_val <- min(filteredPlots%>%
                     pull(.data[[value_y]]))
  }
  
  heatmap_data <<- filteredPlots%>%
    mutate(saveName1=paste0(value_z, sort_type, "_", "heatmap"))%>%
    group_by({{value_z}},saveName1)%>%
    group_nest()%>%
    mutate(heatMap=map2(data,.data[[value_z]],possibly(~makeHeatmap(.x,
                                                               factors1=.x$factors,
                                                               value_y=y,
                                                               value_x=x,
                                                               value_z=.y,
                                                               data_sort = data_sort,
                                                               data_mean = data_mean,
                                                               con_names = con_names,
                                                               sym_breaks = sym_breaks,
                                                               max_val=max_val,
                                                               min_val=min_val),
                                                  NA,quiet=FALSE)
    )
    )
  
  saveHeatmap(heatmap_table=heatmap_data,
              pathway_name=pathway_name,
              filetype=filetype,
              data_date=data_date,
              data_type=data_type,
              data_sort=sort_type)
  
  return(heatmap_data)
}

makeHeatmap <- function(met_data,
                        value_y,
                        value_x,
                        value_z,
                        factors1=c(),
                        con_names=c(),
                        data_sort,data_sort_y=data_sort,data_sort_x=data_sort,
                        data_mean=FALSE,
                        k=NA,
                        sym_breaks=TRUE,
                        max_val,min_val,...){
  

  Calc_pro_levels <<- met_data%>% 
    ungroup()%>%
    mutate(factors2 = factors1%>%
             map(~as.character(interaction(as.list(.), sep='_', lex.order = TRUE)))%>%
             unlist(),
           factors3 = factors1%>% 
             imap(~interaction(met_data[.y, match(.x, names(met_data))], sep="_", lex.order = TRUE))%>%
             unlist())
  
  Summ_pro_levels <- Calc_pro_levels %>%
    group_by(factors3,unlist(factors1),{{value_x}})%>%
    summarise({{value_y}}:=mean({{value_y}},na.rm = T))%>%
    ungroup()
  
  Summ_pro_levels <- Summ_pro_levels%>%
    rename(sample_group = factors3)
  
  
  Calc_pro_levels <- Calc_pro_levels %>% 
    mutate(sample_group = interaction(factors3,
                                      Replicate,sep="_",lex.order=TRUE))
  
  if (data_mean){
    analysis_data <- Summ_pro_levels
    data_type <- "_mean"
  } else {
    analysis_data <-  Calc_pro_levels
    data_type <- NA
  }
  
  
  
  # Hierarchical Clustering by Compound
  mat_pro_levels <- analysis_data%>%     
    select(sample_group,.data[[value_x]],.data[[value_y]])%>% 
    pivot_wider(sample_group,.data[[value_y]])%>%
    remove_missing()
  
  norm.data <<- mat_pro_levels
  #norm.data <- norm.data[rowSums(is.na(norm.data))==0, ]
  
  if(dim(norm.data)[1]>15000){
    k=15000
  }else if(dim(norm.data)[2]>250){
    k=250
  }else{
    k=NA
  }
  
  # save the Metabolite data as a vector
  met <- norm.data%>%
    pull(quo_name(value_x))
  
  # convert the data to a matrix for the heatmap.2() function. 
  #Have to do that only for the numeric data so the Metabolite column is saved as a vector
  # set the Metabolite names as row names of the matrix
  norm.data.matrix <- data.matrix(norm.data[, 2:ncol(norm.data)])
  rownames(norm.data.matrix) <- met
  
  if (data_sort){
    col_dendo <- "correlation"
    row_dendo <- "correlation"
    sort_type <- "_sorted"
    clust_type <- "single"
  } else if (data_sort_x){
    col_dendo <- "correlation"
    row_dendo <- NA
    sort_type <- "_sortedx"
    clust_type <- "single"
  } else if(data_sort_y) {
    col_dendo <- NA
    row_dendo <- "correlation"
    sort_type <- "_sortedy"
    clust_type <- "single"
  } else{
    col_dendo <- NA
    row_dendo <- NA
    sort_type <- "_unsorted"
    clust_type <- NA
  }
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  cluster_cols <- sort_hclust(hclust(dist(t(norm.data.matrix))))
  cluster_rows <- sort_hclust(hclust(dist(norm.data.matrix)))
  
  # set the color palette. n= needs to be at least 1 less than the total of breaks (if the breaks argument is used in the heatmap.2 function)
  heatmap.palette <- colorRampPalette(c("blue", "white", "orange"))(n = 300)
  
  if((max_val==0)&(min_val==0)){
    
    if (sym_breaks){
      myBreaks <- c(seq(-(max(abs(pull(analysis_data,.data[[value_y]])))), 0,
                        length.out=ceiling(300/2) + 1), 
                    seq(max(abs(pull(analysis_data,.data[[value_y]])))/300,
                        max(abs(pull(analysis_data,.data[[value_y]]))),
                        length.out=floor(300/2)))
    }else{
      myBreaks <- c(seq(min(pull(analysis_data,.data[[value_y]])), 0,
                        length.out=ceiling(300/2) + 1), 
                    seq(max(pull(analysis_data,.data[[value_y]]))/300,
                        max(pull(analysis_data,.data[[value_y]])),
                        length.out=floor(300/2)))
    }
  }else{
    myBreaks <- c(seq(min_val,0,length.out = ceiling(300/2) + 1),
                  seq(max_val/300,max_val,length.out = floor(300/2)))
  }
  
  
  
  fontsize=100/((nrow(norm.data.matrix)^0.15))
  
  annotation_col <<- as.data.frame(analysis_data%>%
                                     group_by_at(.vars=c("sample_group",unlist(factors1)))%>%
                                     summarise())
  
  rownames(annotation_col) <- annotation_col[,1]
  annotation_col[,1] <- NULL
  
  gaps_col <- c(ifelse(data_sort,50,0))
  gaps_row <- c(ifelse(data_sort,50,0))
  
  norm.data.heatmap <<- pheatmap(norm.data.matrix,
                                 #main=value_z,
                                 color=heatmap.palette,
                                 #cluster_cols      = cluster_cols,
                                 #cluster_rows      = cluster_rows,
                                 cluster_cols = data_sort_x,
                                 cluster_row = data_sort_y,
                                 clustering_distance_rows = row_dendo,
                                 clustering_distance_cols = col_dendo,
                                 clustering_method = clust_type,
                                 show_rownames = FALSE,
                                 show_colnames = FALSE,
                                 display_numbers = FALSE,
                                 kmeans_k = k,
                                 fontsize=12,
                                 #fontsize_row = fontsize,
                                 annotation_col = annotation_col,
                                 gaps_col = NULL,
                                 gaps_row = NULL,
                                 cutree_rows = 4,
                                 cutree_cols = 4,
                                 #cellheight = 20,
                                 #cellwidth=40,
                                 breaks = myBreaks,
                                 drop_levels = TRUE,
                                 border_color = FALSE
                                 
                                 
  )
  
  return(as.ggplot(norm.data.heatmap))
}

saveHeatmap <- function(heatmap_table,
                        pathway_name=c(),
                        filetype,
                        data_date,
                        data_type,
                        data_sort){
  p <- plot_grid(plotlist = heatmap_table$heatMap,ncol=2,nrow=ceiling((dim(heatmap_table)[1])/2),axis="t")
  
  if((!(is.null(pathway_name)))){
    title <- ggdraw()+
      draw_label(pathway_name,fontface = "bold")
    
    p2 <- plot_grid(title,p,ncol=1,rel_heights=c(0.1,1))
    
    save_plot(paste0(data_date,"_",
                     pathway_name,"_",data_type, 
                     data_sort, "_", "heatmap",filetype),
              p2,
              ncol = 2,
              nrow=ceiling((dim(heatmap_table)[1])/2),
              base_aspect_ratio = 1.3,
              limitsize = FALSE)
  }else{
    save_plot(paste0(data_date,
                     data_type, 
                     data_sort, "_", "heatmap",filetype),
              p,limitsize = FALSE,
              base_height=7*1,base_width = 10*1)
  }
}




