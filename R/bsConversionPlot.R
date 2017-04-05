#' Mehtylation bias plot
#'
#'Plot the bisulfite conversion rate for each sample
#'based on the pheno data in the bs object
#'
#'@param Takes the bs object as the input
#'@return Plot showing bisulfite conversion rate for each sample
#'
#'@examples
#'bsConversionPlott(bs)
#'@export


bsConversionPlot<-function(bs){
  phenoData<-pData(bs)
  if ('bsconversion' %in% colnames(phenoData)) {
    bscDf<-data.frame(sample=rownames(phenoData),bsc=phenoData$bsconversion)

  g<-ggplot2::ggplot(meltedDf,aes(sample,bsc))+geom_point()+ylim(min(meltedDf$bsc)-0.05,max(meltedDf$bsc)+0.05)+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
    xlab('samples')+ylab('bisulfite conversion rate')+ggtitle('Bisulfite conversion rate across samples')
  return(g)
  }else
    message("Provide a bs object with bisufite conversion to produce the plot")
}
