#' Downsample analysis
#'
#'Downsample the CpG covergae matrix for saturation analysis
#'@param bs bsseq object
#'@param dsRates downsampling rate. i.e. the probabaility of sampling a single CpG
#'default is list of probabilities ranging from 0.01 to 1
#'This can be changed by the user, for more continuous saturation curve
#'dsRates can be changed to add more sampling rates
#'@return Data frame with the CpG coverage for each sample at each sampling rate
#'@examples
#'downsample(bsObject)
#'downsample(bsObject,seq(0,1,length.out=20))
#'
#'@export

downsample <- function(bs,dsRates = c(0.01, 0.02, 0.05, seq(0.1, 0.9, 0.1),0.99,1)) {
  covMatrix<-getCoverage(bs)
  Samples<-sampleNames(bs)
  downSampleMatrix<-matrix(nrow=length(dsRates),ncol=length(Samples))

  for (i in 1:length(dsRates)){
    #covSubMatrix<-as.data.frame((unlist(lapply(covMatrix,rbinom,n=1,prob=dsRates[i]))))
    for (j in 1:ncol(covMatrix)){

      cellCoverage<-covMatrix[,j]
      cellNonZeroCoverage<-cellCoverage[cellCoverage>0]
      covSubList<-lapply(cellNonZeroCoverage,rbinom,n=1,prob=dsRates[i])

      #covSubMatrix<-as.data.frame(as.vector(covMatrix))
      #colnames(covSubMatrix)<-c("CpG")
      #covSubMatrix$sample<-Samples
      #covSubMatrix$IndCpG<-as.numeric(covSubMatrix$CpG>0)
      #CpGcount<-table(covSubMatrix$IndCpG,covSubMatrix$sample)
      #idx <- match(colnames(cov_cll),ReadData$samplename)
      downSampleMatrix[i,j]<- sum(covSubList>0)
    }
  }
  rownames(downSampleMatrix)<-dsRates
  return(downSampleMatrix)

}

