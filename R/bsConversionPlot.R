#' Bisulfite conversion rate visualization
#'
#'Plot the bisulfite conversion rate for each sample
#'based on the pheno data in the bs object
#'@param bs bsseq object
#'@return Plot showing bisulfite conversion rate for each sample
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-HDF5Array::loadHDF5SummarizedExperiment(directory)
#'bsConversionPlot(bs)
#'@export


bsConversionPlot<-function(bs){
    phenoData<-bsseq::pData(bs)
    phenoData$bsconversion <- 1 - (phenoData$CHH_meth+phenoData$CHG_meth)/
                                    (phenoData$CHH_meth+phenoData$CHH_unmeth+
                                    phenoData$CHG_meth+phenoData$CHG_unmeth)
    bscDf<-data.frame(sample=rownames(phenoData),bsc=phenoData$bsconversion)

    g<-ggplot2::ggplot(bscDf, ggplot2::aes_string(x="''", y='bsc'))
    g<-g+ggplot2::geom_boxplot()
    g<-g+ggplot2::ylim(max(min(bscDf$bsc)-0.02,0),min(max(bscDf$bsc)+0.02,1))
    g<-g+ggplot2::theme_bw()
    g<-g+ggplot2::geom_jitter()
    g<-g+ggplot2::xlab('')+ggplot2::ylab('bisulfite conversion rate')
    g<-g+ggplot2::ggtitle('Bisulfite conversion rate across samples')
    return(g)

}
