#' Provide graphics for read information
#'
#'Plot the mapped and unmapped reads
#'@param bs bsseq object
#'@return Plot showing the mapped and unmapped read information for each cell
#'@examples
#'directory <- system.file("extdata/bismark_data",package='scmeth')
#'bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
#'readmetrics(bs)
#'@export

readmetrics <- function(bs){

    phenotypicData <- Biobase::pData(bs)
    dat <- data.frame(sample=rownames(phenotypicData),
                    total=as.vector(phenotypicData$total_reads),
                    mapped=as.vector(phenotypicData$uniquely_aligned_reads+
                                phenotypicData$non_uniquely_aligned_reads))
    dat$unmapped <- dat$total-dat$mapped
    return(dat)
}
