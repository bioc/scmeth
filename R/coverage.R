#' Coverage for single cells
#'
#'Provides Coverage for each cell in a library pool
#'@param bs bsseq object
#'@param subSample number of CpGs to subsample
#'@param offset how many CpGs to offset when subsampling
#'@return vector of coverage for the cells in bs object
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-HDF5Array::loadHDF5SummarizedExperiment(directory)
#'coverage(bs)
#'@importFrom DelayedArray colSums
#'@importFrom bsseq getCoverage
#'@export


coverage <- function(bs,subSample=1e6,offset=50000) {
    nCpGs<-nrow(bs)

    if (nCpGs<(subSample+offset)){
        bs<-bs
    }else{
        bs<-bs[offset:(subSample+offset)]
    }

    covMatrix<-bsseq::getCoverage(bs)
    covVec<- DelayedArray::colSums(covMatrix>0,na.rm=TRUE)

    covVec<-covVec*(nCpGs/subSample)
    return(covVec)
}
