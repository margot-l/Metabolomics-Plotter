require("gdtools")
require("dplyr")
require("tidyr")
require("broom")
require("cowplot")
require("ggsignif")
require("rlang")
require("purrr")

signANOVA <- function(normData,sampleData,
                      pVal=0.05,q=c(pVal,(1-pVal)),
                      comp=character(),
                      saveName){
  mets <<- normData%>%
    group_by(normMethod,featureName)%>%
    nest()%>%
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
  
  save_plot(paste0(saveName,"volcano_01.png"),volcanoPlot,
            ncol=dim(data1.cor%>%ungroup()%>%distinct(term))[1],
            nrow=dim(data1.cor%>%ungroup()%>%distinct(normMethod))[1],
            base_height = 2,base_aspect_ratio=1)
  
  sign_mets <- mets%>%
    unnest(mets,.drop=FALSE)%>%
    #filter(!(comparison%in% comp))%>%
    left_join(sampleData%>%
                select(featureName,rtCenter,mzCenter,Kegg.ID)%>%
                distinct(),by="featureName")%>%
    select(normMethod,featureName,term,,Kegg.ID,mzCenter,rtCenter,comparison,estimate,adj.p.value)%>%
    ungroup()%>%
    group_by(term)%>%
    filter(estimate<quantile(estimate,q[1])|
             estimate>quantile(estimate,q[2]))%>%
    filter(adj.p.value<=pVal)
  
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

plotPeaks <- function(plotData, 
                      met_list=c(),integratedset,
                      yPlot="foldChange", xPlot="Contents", fillPlot="Drug",
                      sampleData,cBlacklist=c(),
                      plotP=FALSE,plotPeak=TRUE,plotPlot=TRUE,
                      ...){
  
  if (plotPlot&plotP){
    plotData <- plotData%>%
      filter(ifelse(is.null(met_list),!(featureName%in%met_list),featureName%in%met_list))%>%
      ungroup()%>%
      mutate(factors1 = factors %>% map(~as.character(interaction(as.list(.), sep=':', lex.order = TRUE))) %>% unlist(),
             factors2 = factors %>% imap(~interaction(plotData[.y, match(.x, names(plotData))], sep=":", lex.order = TRUE)) %>% unlist())%>%
      mutate(Contents2=interaction(!!!syms(c(xPlot,fillPlot)),sep="\n",lex.order=TRUE))%>%
      #mutate(KCN2=ifelse(KCN=="0","-",ifelse(KCN=="1","+","++")))
      mutate(KCN2=ifelse(KCN=="0KCN","-","+"))
      #mutate(Contents3=ifelse(Contents=="M9","-\n-",ifelse(Contents=="2DG", "-\n+",ifelse(Contents=="KCN","+\n-","+\n+"))))
    
    plots <- plotData%>%
      ungroup()%>%
      group_by(normMethod,featureName,form,factors1,saveName)%>%
      nest()%>%
      mutate(savedPlot=
               pmap(list(x=data,
                         y=featureName,z=normMethod,
                         a=factors1,b=form),
                    function(x,y,z,a,b)
                      ggplot(data=x,
                             aes_string(y=yPlot,x="factors2"))+
                      geom_boxplot(aes_string(fill=fillPlot,colour=fillPlot),
                                   outlier.shape = NA,alpha=0.5)+
                      #geom_violin(aes_string(fill="Mutant",colour="Mutant"),
                      #            alpha=0.4)+
                      #stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, 
                      #             geom = "crossbar",width=0.5) + 
                      geom_jitter(aes_string(shape=fillPlot),width=0.2)+
                      theme_few(base_size=17)+
                      theme(plot.title = element_text(size = 30,
                                                      hjust=0.5),
                            legend.title = element_text(size=20),
                            legend.text = element_text(size = 20),
                            aspect.ratio=1.0,
                            axis.text =  element_text(size=25),
                            axis.title=element_text(size=27))+
                      #geom_text(aes(label = "KCN", x = -Inf, y = -Inf),size=7,vjust=1.18,hjust=1)+
                      coord_cartesian(clip="off")+
                      labs(title=gsub("_","-",gsub("L_","",gsub("ic_Acid","ate",gsub("ic_acid", "ate",y)))),
                           #subtitle=paste0("m/z:",sampleData%>%
                           #                   select(featureName,rtCenter,mzCenter)%>%
                           #                  distinct()%>%
                           #                 filter(featureName==y)%>%
                           #                pull(mzCenter)," RT:",sampleData%>%
                           #               select(featureName,rtCenter,mzCenter)%>%
                           #              distinct()%>%
                           #             filter(featureName==y)%>%
                           #            pull(rtCenter),"s"),
                           x="KCN",
                           y="log2(FoldChange)"
                           #y=paste0(z,"\nlog2(Fold Change)")
                      )
               )
      )
    
    anovaMets <<- plotData%>%
      ungroup()%>%
      group_by(normMethod,featureName,factors1,saveName)%>%
      nest()%>%
      mutate(mets=map(data,~tidy(TukeyHSD(aov(as.formula(.$form),.)))))%>%
      unnest(mets,.drop=F)
    
    sumMets <- anovaMets%>%
      filter(adj.p.value<=0.1&term==factors1)%>%
      filter(!(comparison%in%cBlacklist))%>%
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
                                                                          filter(!(comparison%in%cBlacklist))%>%
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
                                 #scale_colour_discrete(name="Strain",
                                #                        breaks=c("N2","kynu1"),
                                #                        labels=c("N2","kynu-1(e1003)"),
                                #                        aesthetics=c("colour","fill"))+
                                #  scale_shape(name="Strain",
                                #                        breaks=c("N2","kynu1"),
                                #              labels=c("N2","kynu-1(e1003)"))+
                                 scale_colour_discrete(name="Rotenone",
                                                       breaks=c("0Rot","12.5Rot"),
                                                       labels=c("0.8% DMSO","12.5uM"),
                                                       aesthetics=c("colour","fill"))+
                                 scale_shape(name="Rotenone",
                                             breaks=c("0Rot","12.5Rot"),
                                             labels=c("0.8% DMSO","12.5uM"))+
                                 scale_x_discrete(breaks=x$factors2,labels=x$KCN2)+
                                 #scale_x_discrete(breaks=x$factors2,labels=x$Contents3)+
                                 scale_y_continuous(
                                   labels = scales::number_format(accuracy = 1,
                                                                  decimal.mark = ','))
                              )
      )
      )
  }else if (plotPlot&!plotP){
    plotData <- plotData%>%
      filter(ifelse(is.null(met_list),!(featureName%in%met_list),featureName%in%met_list))%>%
      group_by(featureName,normMethod,saveName)%>%
      nest()
    plots <- plotData%>%
      mutate(savedPlot=
               pmap(list(x=data,y=featureName,z=normMethod),
                    function(x,y,z)ggplot(data=x,
                                          aes_string(y=yPlot,x=xPlot,fill=fillPlot))+
                      geom_boxplot(aes_string(fill=fillPlot),
                                   position = position_dodge(width=0.9))+
                      #geom_point(size=3,
                      #           position=position_jitterdodge(dodge.width=0.9),
                      #           show.legend=FALSE) +
                      theme_few(base_size=17)+
                      theme(plot.title = element_text(size = 18),
                            #legend.key.size = 20,
                            legend.title = element_text(size=20,face="bold"),
                            legend.text = element_text(size = 20),
                            aspect.ratio=0.9,
                            axis.title = element_text(face="bold"))+
                      ggtitle(paste0(y,"\nm/z:",sampleData%>%
                                       select(featureName,rtCenter,mzCenter)%>%
                                       distinct()%>%
                                       filter(featureName==y)%>%
                                       pull(mzCenter)," RT:",sampleData%>%
                                       select(featureName,rtCenter,mzCenter)%>%
                                       distinct()%>%
                                       filter(featureName==y)%>%
                                       pull(rtCenter),"s"))+
                      ylab(paste0(z,"\nlog2(Fold Change)"))
               )
      )
  }
  if (plotPeak){
    plotData <- plotData%>%
      filter(ifelse(is.null(met_list),!(featureName%in%met_list),featureName%in%met_list))%>%
      group_by(featureName,normMethod,saveName)%>%
      nest()
    plots <- plots%>%
      mutate(savedPeak=map2(data,featureName,
                            possibly(~plotEIC(integratedset,sampleName=sampleName,
                                              featureName = .y,featureClass="met",
                                              displayIntegrationRange=TRUE,integrationRangeAlphaScalar=0.8),
                                     NA)
      )
      )
  }
  return(plots)
}

savePeaks <- function(plots,met_list,
                      saveName,
                      widthPlot=7, heightPlot=5,
                      widthPeak=7, heightPeak=5,
                      savePeak=TRUE,savePlot=TRUE,fileType=".png"){
  if(!(is.data.frame(met_list))){
    met_list=as.data.frame(met_list)%>%
      rename(featureName="met_list")
  }
  
  filteredPlots <- plots%>%
    semi_join(met_list)
  
  if (savePlot){
    map2(paste0(filteredPlots$saveName,filteredPlots$featureName,fileType),
         filteredPlots$savedPlot,
         base_height=heightPlot,base_width=widthPlot,
         save_plot) 
  }
  
  if (savePeak){
    map2(paste0(filteredPlots$saveName,filteredPlots$featureName,"_peak",fileType),
         filteredPlots$savedPeak,
         base_width=widthPeak,base_height=heightPeak,
         base_aspect_ratio=0.8,
         save_plot)   
  }
}