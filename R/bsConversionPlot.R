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
    if ('bsconversion' %in% colnames(phenoData)) {
        bscDf<-data.frame(sample=rownames(phenoData),bsc=phenoData$bsconversion)

        g<-ggplot2::ggplot(bscDf,ggplot2::aes_string('sample','bsc'))
        g<-g+ggplot2::geom_point()
        g<-g+ggplot2::ylim(max(min(bscDf$bsc)-0.05,0),min(max(bscDf$bsc)+0.05,1))
        g<-g+ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                            hjust = 1))
        g<-g+ggplot2::xlab('samples')+ggplot2::ylab('bisulfite conversion rate')
        g<-g+ggplot2::ggtitle('Bisulfite conversion rate across samples')
        return(g)
    }else
        warning("Provide a bs object with BS conversion info to produce the plot")
}
