require("gdtools")
require("dplyr")
require("tidyr")
require("broom")
require("cowplot")
require("ggsignif")
require("rlang")
require("purrr")
require("gsubfn")
require("ggthemes")

signANOVA <- function(normData,sampleData,final.cmpd.df,
                      pVal=0.05,q=c(pVal,(1-pVal)),
                      comp=character(),
                      saveName){
  mets <<- normData%>%
    group_by(normMethod,featureName)%>%
    group_nest()%>%
    mutate(mets=map(data,~tidy(TukeyHSD(aov(as.formula(.$form),.),conf.level=(1-pVal)))))
  
  data1 <- mets%>%
    unnest(mets,.drop=FALSE)%>%
    group_by(term)%>%
    mutate(estimate_sign=ifelse(xor(estimate<quantile(estimate,q[1]),
                                    estimate>quantile(estimate,q[2]))
                                &adj.p.value<pVal,1,0))
  data1.cor <- data1%>%
    group_by(term,normMethod)%>%
    tally(estimate_sign==1)
  
  volcanoPlot <- ggplot(data1,
                        aes(x=sign(estimate)*sqrt(abs(estimate)),y=-log10(adj.p.value),color=estimate_sign))+
    geom_point()+
    facet_grid(normMethod~term)+
    theme_few()+
    geom_hline(yintercept=-log10(pVal),linetype="dashed")+
    geom_label(data=data1.cor,aes(x=Inf,y=Inf,label=n,hjust=1.25,vjust=1.25),
               inherit.aes=FALSE,parse=FALSE,label.size=0.15)+
    #geom_vline(xintercept = c(quantile(estimate,q[1]),quantile(estimate,q[2])),linetype="dashed")+
    theme(legend.position = "none",aspect.ratio=1)
  
  save_plot(filename=paste0(saveName,"volcano_01.png"),plot=volcanoPlot,
            ncol=dim(data1.cor%>%ungroup()%>%distinct(term))[1],
            nrow=dim(data1.cor%>%ungroup()%>%distinct(normMethod))[1],
            base_height = 2,base_aspect_ratio=1)
  
  sign_mets <- mets%>%
    unnest(mets,.drop=FALSE)%>%
    #filter(!(comparison%in% comp))%>%
    left_join(sampleData%>%
                select(featureName,rtCenter,mzCenter,KEGG_ID)%>%
                distinct(),by="featureName")%>%
    select(normMethod,featureName,term,KEGG_ID,mzCenter,rtCenter,comparison,estimate,adj.p.value)%>%
    ungroup()%>%
    group_by(term)%>%
    filter(estimate<quantile(estimate,q[1])|
             estimate>quantile(estimate,q[2]))%>%
    filter(adj.p.value<=pVal)%>%
    left_join(final.cmpd.df,by="KEGG_ID")
  
  return(sign_mets)
}

writeTable <- function(writeData,met_list,
                       ys=c("integratedIntensity","norm","foldChange"),
                       spreadFactor="Drug",factors=c("Replicate"),
                       saveName){
  for (i in ys){
    qi=enquo(i)
    qFactor=enquo(spreadFactor)
    col_vector=c("featureName",spreadFactor,factors,i)
    new_norm <- writeData%>%
      filter(featureName %in% met_list)%>%
      select_(.dots=col_vector)%>%
      spread(!!qFactor,!!qi)
    write.csv(new_norm,file=paste0(saveName,i,".csv"))
  }
}

plotPeaks <- function(plotData,integratedset,final.cmpd.df,
                      met_list=c(),pathway_list=c(),
                      yPlot="foldChange", fillPlot="Drug", fillVector=c("DMSO","12.5uM"),
                      sampleData,cBlacklist=c(),
                      plotP=FALSE,plotPeak=TRUE,plotPlot=TRUE,
                      ...){
  if (!(is.null(pathway_list))){
    met_list=plotData%>%
      inner_join(final.cmpd.df%>%
                   filter(Pathway%in%pathway_list),by="KEGG_ID")%>%
      pull(featureName)
  }
  if (!(is.null(met_list))){
    plotData <- plotData%>%
      filter(featureName%in%met_list)%>%
      ungroup()
  }
  
  plotData <- plotData%>%
    mutate(factors1 = factors%>%
             map(~as.character(interaction(as.list(.), sep=':', lex.order = TRUE)))%>%
             unlist(),
           factors2 = factors%>% 
             imap(~interaction(plotData[.y, match(.x, names(plotData))], sep=":", lex.order = TRUE))%>%
             unlist(),
           factors3 = factors%>% 
             imap(~interaction(plotData[.y, match(.x[.x!=fillPlot], names(plotData))], sep="\n", lex.order = TRUE))%>%
             unlist(),
           factors4 = gsub("(?<=\\+)[a-zA-Z]*","",
                           gsubfn("[1-9][a-zA-Z0-9]*", ~ paste(rep("+", substr(n,1,1)), collapse = ""),
                                  gsub("0[a-zA-Z0-9]*","-",factors3)),perl=TRUE),
           factors5 = factors%>% 
             map(~as.character(interaction(as.list(.[.!=fillPlot]), sep='\n', lex.order = TRUE)))%>%
             unlist()
           )
  if (plotPlot){
    plots <- plotData%>%
      ungroup()%>%
      group_by(normMethod,featureName,form,factors1,saveName)%>%
      group_nest()%>%
      mutate(savedPlot=
               pmap(list(x=data,
                         y=featureName,z=normMethod,
                         a=factors1,b=form),
                    function(x,y,z,a,b)
                      ggplot(data=x,
                             aes_string(y=yPlot,x="factors2"))+
                      geom_boxplot(aes_string(fill=fillPlot,colour=fillPlot),
                                   outlier.shape = NA,alpha=0.5)+
                      geom_jitter(aes_string(shape=fillPlot),width=0.2)+
                      theme_few(base_size=17)+
                      theme(plot.title = element_text(size = 30,
                                                      hjust=0.5),
                            legend.title = element_text(size=20),
                            legend.text = element_text(size = 20),
                            aspect.ratio=1.0,
                            axis.text =  element_text(size=25),
                            axis.title=element_text(size=27))+
                      geom_text(aes(label = x$factors5, x = -Inf, y = -Inf),size=8,vjust=1.25,hjust=1.1)+
                      coord_cartesian(clip="off")+
                      labs(title=gsub("_","-",gsub("(_?[DL]+_)","",gsub("(?<=[0-9])_(?=[0-9])",",",gsub("(ic.?(a|A)cid)", "ate",y),perl=TRUE))),
                           y=paste0("log2(",yPlot,")"),
                           x="",
                           caption=z
                      )
               )
      )
    if (plotPlot&plotP){
      anovaMets <- plotData%>%
        ungroup()%>%
        group_by(normMethod,featureName,factors1,saveName)%>%
        group_nest()%>%
        mutate(mets=map(data,~tidy(TukeyHSD(aov(as.formula(.$form),.)))))%>%
        unnest(mets,.drop=F)
      
      sumMets <- anovaMets%>%
        filter(adj.p.value<=0.1&term==factors1)%>%
        filter(comparisonFinder(comparison)==TRUE)%>%
        group_by(normMethod,featureName)%>%
        summarise(n())
      
      plots1 <- suppressWarnings(plots%>%
                                   semi_join(sumMets,by = c("featureName", "normMethod"))%>%
                                   ungroup()%>%
                                   mutate(savedPlot=pmap(list(x=data,
                                                              y=featureName,z=normMethod,
                                                              a=factors1,b=form,
                                                              p=savedPlot),
                                                         function(x,y,z,a,b,p)
                                                           (p+
                                                              geom_signif(data=tidy(TukeyHSD(aov(as.formula(b),data=x)))%>%
                                                                            filter(term==a)%>%
                                                                            filter(comparisonFinder(comparison)==TRUE)%>%
                                                                            separate(comparison,into=c("group1","group2"),sep="-")%>%
                                                                            mutate(p.signif=ifelse(adj.p.value<=0.0001, '****',
                                                                                                   ifelse(adj.p.value<=0.001, "***", 
                                                                                                          ifelse(adj.p.value<=0.01,"**",
                                                                                                                 ifelse(adj.p.value<=0.05,"*",
                                                                                                                        ifelse(adj.p.value<=0.1,"'",
                                                                                                                               "ns"))))))%>%
                                                                            filter(adj.p.value<=0.1)%>%
                                                                            mutate(range=plotData%>%
                                                                                     filter(featureName==y&normMethod==z)%>%
                                                                                     ungroup()%>%
                                                                                     summarise(range=max(foldChange)-min(foldChange))%>%
                                                                                     pull(),
                                                                                   big.max.y=plotData%>%
                                                                                     filter(featureName==y&normMethod==z)%>%
                                                                                     ungroup()%>%
                                                                                     summarise(max(foldChange))%>%
                                                                                     pull())%>%
                                                                            ungroup()%>%
                                                                            add_tally()%>%
                                                                            mutate(y.pos=big.max.y+range * 0.07 + range * 0.1* c(0:(n-1))),
                                                                          aes(xmin=group1, xmax=group2, 
                                                                              annotations=p.signif, 
                                                                              y_position=y.pos),
                                                                          manual=TRUE,vjust=0.6, textsize=10,size=0.7))))
      )
      
      
      plots <- rbind(plots1,anti_join(plots,sumMets,
                                      by = c("normMethod", "featureName")))%>%
        mutate(savedPlot=pmap(list(x=data,
                                   y=featureName,z=normMethod,
                                   a=factors1,b=form,
                                   p=savedPlot),
                              function(x,y,z,a,b,p)
                                (p+
                                   scale_x_discrete(breaks=x$factors2,
                                                    labels=x$factors4)+
                                   scale_y_continuous(
                                     labels = scales::number_format(accuracy = 1,
                                                                    decimal.mark = ','))+
                                   scale_colour_discrete(name=fillPlot,
                                                         breaks=waiver(),
                                                         labels=fillVector,
                                                         aesthetics=c("colour","fill"))+
                                   scale_shape(name=fillPlot,
                                               breaks=waiver(),
                                               labels=fillVector)
                                )
        )
        )
    }else if (plotPlot&!plotP){
      plots <- rbind(plots,anti_join(plots,sumMets,
                                      by = c("normMethod", "featureName")))%>%
        mutate(savedPlot=pmap(list(x=data,
                                   y=featureName,z=normMethod,
                                   a=factors1,b=form,
                                   p=savedPlot),
                              function(x,y,z,a,b,p)
                                (p+
                                   scale_y_continuous(
                                     labels = scales::number_format(accuracy = 1,
                                                                    decimal.mark = ',')))))
    }
  }else if(!plotPeak){
    plots <- plotData%>%
      group_by(featureName,saveName)%>%
      group_nest()
  }
  if (plotPeak){
    plots <- plots%>%
      mutate(savedPeak=map(featureName,
                            possibly(~plotEIC(integratedset,
                                              featureName = .,featureClass="met",
                                              displayIntegrationRange=TRUE,integrationRangeAlphaScalar=0.8),
                                     NA)
      )
      )
  }
  return(plots%>%
           left_join(sampleData%>%
                       select(featureName,KEGG_ID)%>%
                       distinct(),by=c("featureName")))
}

savePeaks <- function(plots,final.cmpd.df,
                      met_list=c(),pathway_list=c(),
                      norm_list=c(),
                      saveName,
                      widthPlot=7, heightPlot=5,
                      widthPeak=7, heightPeak=5,
                      savePeak=TRUE,savePlot=TRUE,
                      path="",
                      fileType=".png"){
  if (!(is.null(pathway_list))){
    met_list=plots%>%
      inner_join(final.cmpd.df%>%
                   filter(Pathway%in%pathway_list),by="KEGG_ID")%>%
      distinct(featureName)
  }else if(!(is.data.frame(met_list))){
    met_list=as.data.frame(met_list)%>%
      rename(featureName="met_list")
  }else{
    met_list=plots%>%
      distinct(featureName)
  }
  if(!(is.null(norm_list))){
    norm_list=as.data.frame(norm_list)%>%
      rename(normMethod="norm_list")
  }else{
    norm_list=plots%>%
      distinct(normMethod)
  }
  
  filteredPlots <- plots%>%
    semi_join(met_list)%>%
    semi_join(norm_list)
  
  if (savePlot){
    map2(paste0(path,filteredPlots$saveName,fileType),
         filteredPlots$savedPlot,
         base_height=heightPlot,base_width=widthPlot,
         save_plot) 
  }
  
  if (savePeak){
    map2(paste0(path,filteredPlots$saveName,"peak_",fileType),
         filteredPlots$savedPeak,
         base_width=widthPeak,base_height=heightPeak,
         base_aspect_ratio=0.8,
         save_plot)   
  }
  dev.off()
}

comparisonFinder <- function(x){
  
  sapply(x, function(x){
    x1 <- strsplit(strsplit(x,"-")[[1]][1],":")[[1]]
    x2 <- strsplit(strsplit(x,"-")[[1]][2],":")[[1]]
                
    if(length(setdiff(x1,x2))==1){
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  
}
