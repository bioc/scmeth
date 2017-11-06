#' Downsample analysis
#'
#'Downsample the CpG covergae matrix for saturation analysis
#'@param bs bsseq object
#'@param subSample number of CpGs to subsample
#'@param dsRates downsampling rate. i.e. the probabaility of sampling
#'a single CpG
#'default is list of probabilities ranging from 0.01 to 1
#'This can be changed by the user, for more continuous saturation curve
#'dsRates can be changed to add more sampling rates
#'@return Data frame with the CpG coverage for each sample at each
#'sampling rate
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'downsample(bs)
#'@importFrom stats dbinom
#'@importFrom bsseq getCoverage
#'
#'@export




downsample <-function(bs,subSample=1e6,dsRates = c(0.01,0.02,0.05, seq(0.1,0.9,0.1))){
    nCpGs<-nrow(bs)
    subSampleCpGs<-min(nCpGs,subSample)
    bs<-bs[1:subSampleCpGs,]


    covMatrix<-bsseq::getCoverage(bs)
    nSamples<-dim(covMatrix)[2]
    downSampleMatrix<-matrix(nrow=length(dsRates)+1,ncol=nSamples)
    maxCov<-20

    nonZeroProbMatrix<-matrix(nrow=(length(dsRates)+1),ncol=maxCov)

    for (i in 1:length(dsRates)){
      nonZeroProbMatrix[i,]<-1 - dbinom(0,1:maxCov,dsRates[i])
    }

    nonZeroProbMatrix[(length(dsRates)+1),]<-1


    countMatrix<-sapply(1:ncol(covMatrix), function(i) {
      cv = as.vector(covMatrix[,i])
      cv[cv>maxCov ] <- maxCov
      tab <- table(cv)
      tab <- tab[names(tab)!="0"]
      x <- rep(0, maxCov)
      x[as.numeric(names(tab))] <- tab
      x
    })

    downSampleMatrix<-round(nonZeroProbMatrix %*% countMatrix)
    downSampleMatrix<-downSampleMatrix*(nCpGs/subSampleCpGs)

    #for (i in 1:length(dsRates)){
    #    for (j in 1:nSamples){
    #    cellCoverage<-as.vector(covMatrix[,j])
    #    cellNonZeroCoverage<-cellCoverage[cellCoverage>0]
    #    covSubList<-lapply(cellNonZeroCoverage,rbinom,n=1,prob=dsRates[i])
    #    downSampleMatrix[i,j]<- sum(covSubList>0)

    #    }
    #}
    #downSampleMatrix[length(dsRates)+1,]<-DelayedArray::colSums(covMatrix>0)
    rownames(downSampleMatrix)<-c(dsRates,1)
    return(downSampleMatrix)
}
