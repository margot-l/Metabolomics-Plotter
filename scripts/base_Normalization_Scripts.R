### Title: Normalization Script
### Summary: A list of normalization functions you must modify to suit your data.

## Guide:

## Begin by naming your function using the '[name] <- function (data){}' format.
##
## Each normalization must have its own function.
##
## The normalization code should be written in the brackets using 'data' as the base dataframe.
##
## Use dplyr and tidyr to make each to read code, see 
## https://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf for guide.
##
## Begin by grouping the factors you want to include in your normalization ie normalizing for 
## each feature, rep, factor. Use 'group_by(x,y)'.
##
## Chose the mathemetical normalization you'll be imploying such as division, log2. After
## the last step, you'll be use 'norm'as your base peak. If you want to normalize by the WT you
## would perhaps write 'norm[Mutant=="WT"]'. If you want multiple specifications it might be
## 'norm[Mutant=="WT"&Contents=="Buffer"]'. See different examples throughout this base script.
## 
## Next create a new column called 'factors' using the 'mutate' function. You want to
## have the following syntax 'list(c("factor1","factor2"))'. These are the factors which you 
## are interested in analyzing their effect on each peak.
## 
## Finally write out the formula for the linear model. ie 'form=foldChange~factor1*factor2'. An
## additive model will not assess synergy while a multiplicative one will. All the variables 
## in your function must be factors you listed in 'factors'.
##
## Be sure to include every normalization function you create in 'list_model' at the end.

require("tidyverse")

rationorm <- function(data){
  data%>%
    # This generates a ratio for each drug condition of kynu-1 mutant relative to the N2.
    group_by(featureName,Contents,Replicate)%>% 
    mutate(foldChange=log2((norm))-log2(mean(norm[Mutant=="N2"])))%>%
    ungroup()%>%
    filter(Mutant=="kynu1")%>%
    mutate(factors=list(c("KCN")))%>%
    # This is a fast and easy way to generalte your formula.
    mutate(form=paste0('foldChange~"',paste0(factors[[1]],collapse = "*")))
  }
mutants <- function(data){
  data%>%
    # This normalization normalizes everything to the WT+buffer control.
    group_by(featureName)%>% 
    mutate(foldChange=log2(norm)-log2(mean(norm[Mutant=="N2"&Contents=="M9"])))%>%
    ungroup()%>%
    mutate(factors=list(c("KCN","Mutant")))%>%
    mutate(form=paste0('foldChange~',paste0(factors[[1]],collapse = "*")))%>%
    mutate(Mutant=factor(Mutant,levels=c("N2","kynu1")))
  
}
mutantnorm <- function(data){
  data%>%
    # This us an extra normalization step which normalizes all the peaks in a single sample
    # to conrol for any sample size variation. Other than this step, this normalization method
    # is the same as 'mutants'.
    filter(Rot=="0Rot")%>%
    filter(KCN!="1KCN")%>%
    group_by(sampleName)%>% 
    mutate(norm1=norm/median(norm))%>%
    ungroup()%>%
    group_by(featureName)%>% 
    mutate(foldChange=log2(norm1)-log2(mean(norm1[Mutant=="N2"&Contents=="M9"])))%>%
    ungroup()%>%
    mutate(factors=list(c("KCN","Mutant")))%>%
    mutate(form=paste0('foldChange~',paste0(factors[[1]],collapse = "*")))%>%
    mutate(Mutant=factor(Mutant,levels=c("N2","kynu1")))
  
}

mutantnorm2 <- function(data){
  data%>%
    # This normalizes to each mutant for each replicate so normalizes for strain density
    # variation.
    group_by(Replicate,Mutant, Drug)%>% 
    mutate(norm1=norm/median(norm))%>%
    ungroup()%>%
    group_by(featureName)%>% 
    mutate(foldChange=log2(norm1)-log2(mean(norm1[Mutant=="N2"&Contents=="M9"])))%>%
    ungroup()%>%
    mutate(factors=list(c("KCN","Mutant")))%>%
    mutate(form=paste0('foldChange~',paste0(factors[[1]],collapse = "*")))%>%
  mutate(Mutant=factor(Mutant,levels=c("N2","kynu1")))
  
}
mutantsmean <- function(data){
  data%>%
    group_by(sampleName)%>% 
    mutate(foldChange=norm/median(norm))%>%
    ungroup()%>%
    # This normalizes the the mean buffer only between all the mutants.
    group_by(featureName,Mutant)%>% 
    mutate(foldChange=log2(foldChange)-log2(mean(foldChange[Contents=="M9"])))%>%
    ungroup()%>%
    mutate(factors=list(c("KCN","Mutant")))%>%
    mutate(form=paste0('foldChange~',paste0(factors[[1]],collapse = "*")))%>%
    mutate(Mutant=factor(Mutant,levels=c("N2","kynu1")))
  
}
repmean <- function(data){
  data%>%
    group_by(sampleName)%>% 
    mutate(foldChange=norm/median(norm))%>%
    ungroup()%>%
    # This normalizes to the mean buffer only for every mutant and replicate.
    group_by(featureName,Mutant,Replicate)%>% 
    mutate(foldChange=log2(foldChange)-log2(mean(foldChange[Contents=="M9"])))%>%
    ungroup()%>%
    mutate(factors=list(c("KCN","Mutant")))%>%
    mutate(form=paste0('foldChange~',paste0(factors[[1]],collapse = "*")))%>%
    mutate(Mutant=factor(Mutant,levels=c("N2","kynu1")))
  
}

## All the models must be listed here. Hash out models if you don't want them to be used.
list_model <- 
  list(
    #mutants=mutants,
    mutantnorm=mutantnorm
    #mutantnorm2=mutantnorm2,
    #mutantsmean=mutantsmean,
    #repmean=repmean
  )

fn_model <- function(.normalizedData,df){
  df%>%
    mutate(normData=map(data,.normalizedData))
  }