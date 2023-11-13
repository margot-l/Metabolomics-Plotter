require("ggthemes")
require("rlang")

cutoffFun <- function(dataFrame,quo_peak,sClass="blank",saveName,percentile){
  cutoff <- quantile(dataFrame%>%
                       filter(sampleClass==sClass&featureClass!="ref")%>%
                       pull(!!quo_peak),probs=c(percentile),na.rm=T)
  ggplot(dataFrame%>%
           filter(sampleClass==sClass&featureClass!="ref"),
         aes_(quo_peak))+
    geom_histogram(bins=100)+
    geom_vline(xintercept = cutoff)+
    theme_few()+
    scale_y_log10()+
    ylab("log10(count)")
  ggsave(paste0(format(Sys.Date(),"%Y%m%d"),"_",toString(quo_peak),saveName,".png"))
  
  #Check the percentile cutoff in the histogram. When you are happy with it, run the following.
  shortenedData <- dataFrame%>%
    select(-intensity)%>%
    group_by(featureName)%>%
    filter(sum((!!quo_peak)[sampleClass==sClass]>cutoff)==0)
 
  return(shortenedData)
}

newcutoffFun <- function(dataFrame,quo_peak,
                         sClass="blank",saveName,
                         percentile,
                         q1=1,q2=1){
  shortenedData <- dataFrame%>%
    group_by(featureName)%>%
    filter((mean((!!quo_peak)[featureClass!="ref"&sampleClass==sClass&
                                (!!quo_peak)<=quantile((!!quo_peak)[sampleClass==sClass&featureClass!="ref"],q1)])/
              mean((!!quo_peak)[featureClass!="ref"&sampleClass=="sample"&
                                  (!!quo_peak)>=quantile((!!quo_peak)[sampleClass=="sample"&featureClass!="ref"],q2)]))<percentile)%>%
    filter((!!quo_peak)!=0)
  
  ggplot(dataFrame%>%
           filter(featureClass!="ref")%>%
           group_by(featureName)%>%
           mutate(ratio=(mean((!!quo_peak)[sampleClass==sClass&
                                             (!!quo_peak)<=quantile((!!quo_peak)[sampleClass==sClass],q1)])/
                           mean((!!quo_peak)[sampleClass=="sample"&
                                               (!!quo_peak)>=quantile((!!quo_peak)[sampleClass=="sample"],q2)])))%>%
           filter(sampleClass=="sample"),
         aes(x=ratio))+
    geom_histogram(bins=100)+
    theme_few()+
    geom_vline(xintercept = percentile)
  ggsave(paste0(format(Sys.Date(),"%Y%m%d"),"_",toString(quo_peak),"_",saveName,".png"),
         width=4,height=4)
  
  return(shortenedData)
}

horizon <- function(dataFrame,cut_height,group_v){
  horizonSet <- dataFrame%>%
    filter(sampleClass=="sample"&featureClass=="met")%>%
    group_by_at(group_v)%>%
    summarise(max_i2=max(maxIntensity))%>%
    filter(max_i2<cut_height)
  
  cutSet <- dataFrame%>%
    anti_join(horizonSet,by=group_v)
  
  return(cutSet)
    
}

dataCleanup <- function(oDF,sampleTable,compounds,
                        peak=maxIntensity,
                        percentile=0.9,percentile2=percentile,
                        cut_height=NULL,
                        method="new",refNorm=FALSE,
                        factors=c(),
                        group_v=c('featureName'),
                        setexp=paste0(factors,collapse="_"),
                        KEGG=TRUE){
  quo_peak <- enquo(peak)
  
  if(KEGG==TRUE){
    oDF <- oDF%>%
      left_join(compounds%>%
                  select(featureName,featureClass,KEGG_ID)%>%
                  mutate(featureClass=tolower(featureClass),
                    KEGG_ID=forcats::fct_explicit_na(KEGG_ID)),
                by=c("featureName","featureClass"))
  }else{
    oDF <- oDF%>%
      left_join(compounds%>%
                  select(featureName,featureClass),
                by=c("featureName","featureClass"))
  }
  
  
  if(!(is.null(cut_height))){
    cutData <- horizon(oDF,cut_height,group_v)
  }
  else{
    cutData <- oDF
  }
  
  if (method=="old"){
    #Remove features which fall outside of a set distribution of the peaks in the blank.
    shortData <- cutoffFun(cutData,quo_peak,sClass="blank",
                           saveName="quant_runblanks",percentile=percentile)
    
    shorterData <- cutoffFun(shortData,quo_peak,sClass="filter",
                             saveName="quant_repblanks",percentile=percentile2)
    
  }else if (method=="new"){
    shortData <<- newcutoffFun(cutData,quo_peak,sClass="blank",
                              saveName="ratio_runblanks",percentile=percentile)

    shorterData <- newcutoffFun(shortData,quo_peak,sClass="filter",
                                saveName="ratio_repblanks",percentile=percentile2)
  }else{
    print ("No method of that name exists; try 'old' or 'new'.")
  }
  
  #Remove features which don't have features in all samples.
  maxSamples <- shorterData%>%
    ungroup()%>%
    filter(sampleClass=="sample")%>%
    group_by(featureName)%>%
    mutate(n=n())%>%
    ungroup()%>%
    summarise(max(n))%>%
    pull()

  sampleData <- shorterData%>%
    ungroup()%>%
    filter(sampleClass=="sample")%>%
    group_by(featureName)%>%
    mutate(n=n())%>%
    filter(sum(n<maxSamples)==0)%>%
    select(-n)%>%
    inner_join(sampleTable,by=c("sampleName","sampleClass"))%>%
    ungroup()%>%
    mutate_at(factors,funs('contrasts<-'(as.factor(.),,contr.treatment)))

  
  if (refNorm){
    sampleData <- sampleData%>%
      group_by(featureName,sampleName)%>%
      mutate(norm=((!!quo_peak)[featureClass!="ref"])/((!!quo_peak)[featureClass=="ref"]))%>%
      filter(featureClass!="ref")%>%
      ungroup()%>%
      mutate(norm=ifelse(norm<=0,1,norm))
    saveName <<- paste0(format(Sys.Date(),"%Y%m%d"),"_",setexp,"_",as_name(quo_peak),"_",method,"_norm_")
  }else{
    sampleData <- sampleData%>%
      filter(featureClass!="ref")%>%
      group_by(featureName,sampleName,featureClass)%>%
      mutate(norm=(!!quo_peak))%>%
      ungroup()%>%
      mutate(norm=ifelse(norm<0,1,norm))
    saveName <<- paste0(format(Sys.Date(),"%Y%m%d"),"_",setexp,"_",toString(quo_peak),"_",method,"_nonorm_")
  }
  
  filteredMets <<- anti_join(oDF,sampleData,by="featureName")%>%
    distinct(featureName)%>%
    pull(featureName)
  
  return(sampleData)
}