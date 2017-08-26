#' Coverage for single cells
#'
#'Provides Coverage for each cell in a library pool
#'@param bs bsseq object
#'@return vector of coverage for the cells in bs object
#'@examples
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'coverage(bs)
#'@importFrom DelayedArray colSums
#'@importFrom bsseq getCoverage
#'@export


coverage <- function(bs) {
    covMatrix<-bsseq::getCoverage(bs)
    covVec<- DelayedArray::colSums(covMatrix>0,na.rm=TRUE)
    return(covVec)
}
