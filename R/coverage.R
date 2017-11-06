#' Coverage for single cells
#'
#'Provides Coverage for each cell in a library pool
#'@param bs bsseq object
#'@param subSample number of CpGs to subsample
#'@return vector of coverage for the cells in bs object
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'coverage(bs)
#'@importFrom DelayedArray colSums
#'@importFrom bsseq getCoverage
#'@export


coverage <- function(bs,subSample=1e6) {
    nCpGs<-nrow(bs)
    subSampleCpGs<-min(nCpGs,subSample)
    covMatrix<-bsseq::getCoverage(bs[1:subSampleCpGs,])
    covVec<- DelayedArray::colSums(covMatrix>0,na.rm=TRUE)

    covVec<-covVec*(nCpGs/subSampleCpGs)
    return(covVec)
}
