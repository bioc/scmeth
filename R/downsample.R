#' Downsample metrics
#'
#'Downsample the CpG covergae matrix for saturation analysis
#'@param takes bs object as input
#'@return Data frame with the saturation Data
#'@example
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
  return(downSampleMatrix)

}
