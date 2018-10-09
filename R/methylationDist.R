#' Provide graphics for methylation distribution
#'
#'Plot the methylation distribution for the cells in bsseq object
#'@param bs bsseq object
#'@return mean methylation for each sample
#'@examples
#'directory <- system.file("extdata/bismark_data", package='scmeth')
#'bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
#'methylationDist(bs)
#'@importFrom bsseq getCoverage
#'@importFrom DelayedArray colMeans
#'@export


methylationDist <- function(bs){
    covMatrix <- bsseq::getCoverage(bs)
    methMatrix <- bsseq::getCoverage(bs, type='M')/covMatrix
    methVec <- DelayedArray::colMeans(methMatrix, na.rm=TRUE)
    return(methVec)
}
