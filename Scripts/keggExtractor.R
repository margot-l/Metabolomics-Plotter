# Title: KEGG Extractor
# Author: Margot Lautens
# Data Created: March 14, 2019
# Date Modified: 
# Summary: This function will download and save a data table of all pathways for a KEGG organism 
#          of interest. 
# Adapted from "20170821_KEGG_ID_for_metabolites" by Soumaya Zlitni.

# Required libraries
require("KEGGREST")     # client interface to the KEGGREST server
require("magrittr")     # for syntax used in writing the function
require("tidyverse")    # for syntax used in writing the function

keggExtractor <- function(organism,
                          save=FALSE){
  
  #
  base.list <- keggList("pathway",organism)
  
  pathwayList <- as_tibble(base.list,rownames="Pathway.ID")%>%
    mutate(Pathway.ID=gsub("path:","",Pathway.ID))%>%
    pull(Pathway.ID)
  
  
  # Make a list of dataframes of compounds, one for each module provided in myPathways
  cmpd.list <- lapply(pathwayList, pathway.cmpds)
  
  # Make the full data frame.
  final.cmpd.df <- as.tibble(do.call(rbind, cmpd.list))%>%
    separate(Pathway,into=c("Pathway","Organism"),sep=" - ",extra="merge")
  
  # Save the final dataframe
  if (save){
    write.csv(final.cmpd.df, file = paste0("final.cmpd.df", ".csv", sep = ""), row.names = FALSE)
  }
  
  return(final.cmpd.df)
}

pathway.cmpds <- function(x) {
  if (!is.null(keggGet(x)[[1]]$COMPOUND)){
    df <- sapply(keggGet(x), "[[", "COMPOUND") %>%
      as.data.frame(.) %>%
      rownames_to_column(.) %>%
      cbind(., sapply(keggGet(x), "[[", "NAME"))
    colnames(df) <- c("KEGG_ID", "Compound", "Pathway")
    return(df)
  }
  
  
}