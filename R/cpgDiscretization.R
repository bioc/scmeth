#' Discretize the CpG methylation values
#' to align with single cell analysis
#'
#'In single cell analysis overwhelmingly large number of CpGs have binary
#'methylation
#'Due to errors in sequencing and amplification many CpGs tend to have
#'non-binary methylation. Hence
#'this function catergorizes the non-binary CpGs as methylated if the
#'methyation is above 0.8 and
#'unmethylated if the methylation is below 0.2
#'@param bs bsseq object
#'@param subSample number of CpGs to subsample
#'@param coverageVec If coverage vector is already calculated provide it to
#'speed up the process
#'@return meth discretized methylation matrix
#'@return discard total number of removed CpGs from each sample
#'@return Percentage of CpGs discarded compared to the total number of CpGs
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'cpgDiscretization(bs)
#'@importFrom DelayedArray colSums
#'@importFrom bsseq getCoverage
#'@export


cpgDiscretization<-function(bs,subSample=1e6, coverageVec=NULL){
    # subsampling
    nCpGs<-nrow(bs)
    subSampleCpGs<-min(nCpGs,subSample)
    bs<-bs[1:subSampleCpGs,]

    covMatrix<-bsseq::getCoverage(bs)
    methMatrix<-bsseq::getCoverage(bs,type='M')
    nSamples<-ncol(methMatrix)
    methMatrix<-methMatrix/covMatrix
    if (is.null(coverageVec)){
      covVec<- DelayedArray::colSums(covMatrix>0,na.rm=TRUE)
    }else{
      covVec<-coverageVec
    }
    #methMatrix[methMatrix>=0.8]<-1
    #methMatrix[methMatrix<=0.2]<-0

    # Only consider methylation between 0.2 and 0.8
    methCutOff<-c(0.01,0.19,0.79,1.0)

    methylationDistMatrix<-sapply(1:nSamples, function(i) {
      mv = as.vector(methMatrix[,i])
      mv<-mv[!is.na(mv)]
      mvBin<-cut(mv,methCutOff)
      tab <- table(mvBin)
      x<- tab
      x
    })


    removedCpGs<-methylationDistMatrix[2,]
    removedCpGFrac<-(removedCpGs/(covVec))*100
    returnList<-list('discard' = removedCpGs,
                     'discardPerc' = removedCpGFrac)
    return(returnList)

}
