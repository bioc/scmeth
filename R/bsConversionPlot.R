#' Bisulfite conversion rate visualization
#'
#'Plot the bisulfite conversion rate for each sample
#'based on the pheno data in the bs object
#'@param bs bsseq object
#'@return Plot showing bisulfite conversion rate for each sample
#'@examples
#'directory <- system.file("extdata/bismark_data",package='scmeth')
#'bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
#'bsConversionPlot(bs)
#'@export


bsConversionPlot <- function(bs){
    phenoData <- bsseq::pData(bs)
    phenoData$bsconversion <- 1 - (phenoData$CHH_meth+phenoData$CHG_meth)/
                                    (phenoData$CHH_meth+phenoData$CHH_unmeth+
                                    phenoData$CHG_meth+phenoData$CHG_unmeth)
    bscDf <- data.frame(sample=rownames(phenoData),bsc=phenoData$bsconversion)
    return(bscDf)
}
