#' Bisulfite conversion rate visualization
#'
#'Plot the bisulfite conversion rate for each sample
#'based on the pheno data in the bs object
#'
#'@param bs bsseq object
#'@return Plot showing bisulfite conversion rate for each sample
#'
#'@examples
#'bsConversionPlot(bsseqObject)
#'@export


bsConversionPlot<-function(bs){
  phenoData<-bsseq::pData(bs)
  if ('bsconversion' %in% colnames(phenoData)) {
    bscDf<-data.frame(sample=rownames(phenoData),bsc=phenoData$bsconversion)

    if (requireNamespace("ggplot2",quietly = TRUE)){
    g<-ggplot(bscDf,aes_string('sample','bsc'))+geom_point()+ylim(max(min(bscDf$bsc)-0.05,0),min(max(bscDf$bsc)+0.05,1))+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))+xlab('samples')+ylab('bisulfite conversion rate')+
    ggtitle('Bisulfite conversion rate across samples')
    return(g)
    }else{
      warning("ggplot package needed for plot rendering")
    }
  }else
    message("Provide a bs object with bisufite conversion to produce the plot")
}
