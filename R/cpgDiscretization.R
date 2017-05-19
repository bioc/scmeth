#' Discretize the CpG methylation values
#' to align with single cell analysis
#'
#'In single cell analysis overwhelmingly large number of CpGs have binary methylation
#'Due to errors in sequencing and amplification many CpGs tend to have non-binary methylation. Hence
#'this function catergorizes the non-binary CpGs as methylated if the methyation is above 0.8 and
#'unmethylated if the methylation is below 0.2
#'@param bs bsseq object
#'@return meth discretized methylation matrix
#'@return discard total number of removed CpGs from each sample
#'@return Percentage of CpGs discarded compared to the total number of CpGs
#'@examples
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'cpgDiscretization(bs)
#'@export
#

cpgDiscretization<-function(bs){
    covMatrix<-bsseq::getCoverage(bs)
    methMatrix<-bsseq::getCoverage(bs,type='M')
    methMatrix<-methMatrix/covMatrix
    tempMethylationMatrix<-methMatrix
    methMatrix[methMatrix<=0.2]<-0
    methMatrix[methMatrix>=0.8]<-1
    removedCpGs<-colSums(methMatrix>0.2 & methMatrix<0.8, na.rm=TRUE)
    removedCpGFrac<-(removedCpGs/(scmeth::coverage(bs)))*100
    returnList<-list('meth' = methMatrix, 'discard' = removedCpGs,
                   'discard-perc' = removedCpGFrac)
    return(returnList)
}
