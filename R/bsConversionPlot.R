#' Bisulfite conversion rate visualization
#'
#'Plot the bisulfite conversion rate for each sample
#'based on the pheno data in the bs object
#'@param bs bsseq object
#'@return Plot showing bisulfite conversion rate for each sample
#'@examples
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'bsConversionPlot(bs)
#'@export


bsConversionPlot<-function(bs){
    phenoData<-bsseq::pData(bs)
    phenoData$bsconversion <- 1 - (phenoData$CHH_meth+phenoData$CHG_meth)/
                                  (phenoData$CHH_meth+phenoData$CHH_unmeth+
                                   phenoData$CHG_meth+phenoData$CHG_unmeth)
    bscDf<-data.frame(sample=phenoData$cell_id,bsc=phenoData$bsconversion)

    g<-ggplot2::ggplot(bscDf, aes(x="", y=bsc))
    g<-g+ggplot2::geom_boxplot()
    g<-g+ggplot2::ylim(max(min(bscDf$bsc)-0.02,0),min(max(bscDf$bsc)+0.02,1))
    g<-g+ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                            hjust = 1))
    g<-g+ggplot2::xlab('')+ggplot2::ylab('bisulfite conversion rate')
    g<-g+ggplot2::ggtitle('Bisulfite conversion rate across samples')
    return(g)

}
