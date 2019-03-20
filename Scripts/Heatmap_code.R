# Margot Lautens
# Date created: 20170926
# Last edited: 20170928
# Code for plotting heat maps in R
# Input is in the format of a tidy (longform) dataset of metabolite and sample
# Adapted from code by Soumaya Zlitni and Julia Hanchard

require("dplyr")
require("tidyr")
require("gplots")
require("pheatmap")
require("grid")

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
                         y="foldChange",
                         norm_list=c(),
                         met_list=c(),
                         data_date=format(Sys.Date(),"%Y%m%d"),
                         con_names=c(),
                         data_sort=TRUE,
                         data_mean=FALSE){
  
  if (!(is.null(met_table))){
    met_data <- met_data%>%
      semi_join(met_table)
  }
  
  if(!(is.null(norm_list))){
    norm_list <- as.data.frame(norm_list)%>%
      rename(normMethod="norm_list")
    
    met_data <- met_data%>%
      semi_join(norm_list)
  }
  
  if(!(is.null(met_list))){
    met_list <- as.data.frame(met_list)%>%
      rename(featureName="met_list")
    
    met_data <- met_data%>%
      semi_join(met_list)
  }
  
  met_data%>%
    group_by(normMethod)%>%
    nest()%>%
    mutate(heatMap=map(data,possibly(~makeHeatmap(.,factors=unlist(.$factors),
                                         value=y,
                                         data_name = .$saveName,
                                         data_sort = data_sort,
                                         data_mean = data_mean,
                                         con_names = con_names),
                                     NA)
                       )
           )
}

makeHeatmap <- function(met_data,
                        value,
                        factors=c(),
                        data_name,
                        data_date=format(Sys.Date(),"%Y%m%d"),
                        con_names=c(),
                        data_sort,data_sort_y=data_sort,data_sort_x=data_sort,
                        data_mean=FALSE,
                        k=NA,...){

  value <- as.symbol(value)
  value <- enquo(value)
  
  Calc_pro_levels <- met_data%>% 
    ungroup()
  
  Summ_pro_levels <- Calc_pro_levels %>%
    group_by_at(.vars=(c("featureName",factors)))%>%
    summarise(!!(quo_name(value)):=mean((!!value),na.rm = T))%>%
    ungroup()
  
  Summ_pro_levels <- Summ_pro_levels %>% 
    unite_("sample_group",factors,remove=FALSE)
  
  if(!data_sort_x){
    Summ_pro_levels$sample_group <- factor(Summ_pro_levels$sample_group,
                                           levels=con_names)
  }

  Calc_pro_levels <- Calc_pro_levels %>% 
    unite_("sample_group",c(factors,"Replicate"),remove=FALSE)
  
  if (data_mean){
    analysis_data <- Summ_pro_levels
    data_type <- "_mean"
  } else {
    analysis_data <-  Calc_pro_levels
    data_type <- NA
  }
  
  if(dim(analysis_data)[1]>15000){
    k=15000
  }else{
    k=NA
  }
  
  # Hierarchical Clustering by Compound
  mat_pro_levels <- analysis_data%>%     
    select(sample_group,featureName,!!quo(quo_name(value)))%>% 
    spread_("sample_group",quo_name(value))
  
  norm.data <- mat_pro_levels
  norm.data <- norm.data[rowSums(is.na(norm.data))==0, ]
  
  # save the Metabolite data as a vector
  met <- norm.data$featureName
  
  # convert the data to a matrix for the heatmap.2() function. 
  #Have to do that only for the numeric data so the Metabolite column is saved as a vector
  # set the Metabolite names as row names of the matrix
  norm.data.matrix <- data.matrix(norm.data[, 2:ncol(norm.data)])
  rownames(norm.data.matrix) <- met
  
  if (data_sort){
    col_dendo <- "correlation"
    row_dendo <- "correlation"
    sort_type <- "_sorted"
  } else if (data_sort_x){
    col_dendo <- "correlation"
    row_dendo <- NA
    sort_type <- "_sortedx"
  } else if(data_sort_y) {
    col_dendo <- NA
    row_dendo <- "correlation"
    sort_type <- "_sortedy"
  } else{
    col_dendo <- NA
    row_dendo <- NA
    sort_type <- "_unsorted"
  }
  
  # set the color palette. n= needs to be at least 1 less than the total of breaks (if the breaks argument is used in the heatmap.2 function)
  heatmap.palette <- colorRampPalette(c("blue", "white", "orange"))(n = 300)
  
  myBreaks <- c(seq(min(pull(analysis_data,!!quo(quo_name(value)))), 0, length.out=ceiling(300/2) + 1), 
                seq(max(pull(analysis_data,!!quo(quo_name(value))))/300, max(pull(analysis_data,!!quo(quo_name(value)))), length.out=floor(300/2)))
  
  
  fontsize=100/((nrow(norm.data.matrix)^0.7))
  
  annotation_col <- as.data.frame(analysis_data%>%
    group_by_at(.vars=c("sample_group",factors))%>%
    summarise())
  
  rownames(annotation_col) <- annotation_col[,1]
  annotation_col[,1] <- NULL
  
  gaps_col <- c(ifelse(data_sort,50,0))
  gaps_row <- c(ifelse(data_sort,50,0))
  
  # create a png for the heat map
  png(paste0(data_name, data_type, sort_type, "_", data_date, "_", "heatmap.png"),       
      width = 8,        # 8 x 300 pixels
      height = 15,
      res = 500,            # 300 pixels per inch
      bg = "transparent",
      units = "in",
      pointsize = 8)        # smaller font size
  
  #par(mar = c(3,3,3,3))
  norm.data.heatmap <- pheatmap(norm.data.matrix,
                                color=heatmap.palette,
                                cluster_cols = data_sort_x,
                                cluster_row = data_sort_y,
                                clustering_distance_rows = row_dendo,
                                clustering_distance_cols = col_dendo,
                                clustering_method = "ward.D2",
                                show_rownames = TRUE,
                                show_colnames = TRUE,
                                display_numbers = FALSE,
                                kmeans_k = k,
                                fontsize_row = fontsize,
                                annotation_col = annotation_col,
                                gaps_col = NULL,
                                gaps_row = NULL,
                                cutree_rows = 2,
                                cutree_cols = 2,
                                breaks = myBreaks
    
  )
    
  dev.off()
  #dev.off()
  return(norm.data.heatmap)
}

