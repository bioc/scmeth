#' Coverage for single cells
#'
#'Provides Coverage for each cell in a library pool
#'@param bs bsseq object
#'@param subSample number of CpGs to subsample.
#'Default value is 1000000.
#'@param offset how many CpGs to offset when subsampling
#'Default value is set to be 50000, i.e. first 50000 CpGs will
#'be ignored in subsampling.
#'@return vector of coverage for the cells in bs object
#'@examples
#'directory <- system.file("extdata/bismark_data",package='scmeth')
#'bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
#'coverage(bs)
#'@importFrom DelayedArray colSums
#'@importFrom bsseq getCoverage
#'@export


coverage <- function(bs,subSample=1e6,offset=50000) {
    nCpGs <- nrow(bs)

    if (subSample == 'all'){
        bs <- bs
        ratio <- 1
    }else{
        if (nCpGs < (subSample + offset)){
            bs <- bs
            subSample <- nCpGs
            ratio <- 1
        }else{
            bs <- bs[offset:(subSample+offset)]
            ratio <- nCpGs/subSample
        }
    }
    covMatrix <- bsseq::getCoverage(bs)
    covVec <- DelayedArray::colSums(covMatrix>0,na.rm=TRUE)
    covVec <- covVec*ratio
    return(covVec)
}
