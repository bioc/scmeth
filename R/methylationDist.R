#' Provide graphics for methylation distribution
#'
#'Plot the methylation distribution for the cells in bsseq object
#'@param bs bsseq object
#'@param subSample number of CpGs to subsample
#'Default value is 1000000.
#'@param offset how many CpGs to offset when subsampling
#'Default value is set to be 50000, i.e. first 50000 CpGs will
#'be ignored in subsampling.
#'@param coverageVec If coverage vector is already calculated provide it to
#'speed up the process
#'@return plot of the methylation distribution
#'@examples
#'directory <- system.file("extdata/bismark_data",package='scmeth')
#'bs <- HDF5Array::loadHDF5SummarizedExperiment(directory)
#'methylationDist(bs)
#'@importFrom bsseq getCoverage
#'@importFrom ggthemes theme_tufte
#'@export


methylationDist <- function(bs,subSample=1e6, offset=50000,coverageVec=NULL){
    # subsampling
    nCpGs <- nrow(bs)

    if (nCpGs<(subSample+offset)){
        bs <- bs
        subSample <- nCpGs
    }else{
        bs <- bs[offset:(subSample+offset)]
    }

    covMatrix <- bsseq::getCoverage(bs)
    methMatrix <- bsseq::getCoverage(bs,type='M')/covMatrix
    nSamples <- ncol(methMatrix)

    methCutOff <- c(-0.01,0.2,0.4,0.6,0.8,1.0)
    methylIntervals <- length(methCutOff)-1
    if (is.null(coverageVec)){
        totCpGs <- DelayedArray::colSums(covMatrix>0,na.rm=TRUE)
    }else{
        totCpGs <- coverageVec
    }

    methylationDistMatrix <- sapply(seq_len(nSamples), function(i) {
        mv = as.vector(methMatrix[,i])
        mv <- mv[!is.na(mv)]
        mvBin <- cut(mv,methCutOff)
        tab <- table(mvBin)
        x <- tab
        x
    })
    #
    methylationDistMatrix <- t(methylationDistMatrix)
    methylationDistMatrix <- methylationDistMatrix*(nCpGs/subSample)

    methylationDistMatrix <- apply(methylationDistMatrix,2,function(x) x/totCpGs)
    orderdMeth <- order(methylationDistMatrix[,1])
    methylationDistMatrix <- methylationDistMatrix[orderdMeth,]

    colnames(methylationDistMatrix) <- c('[0,0.2]','(0.2,0.4]','(0.4,0.6]','(0.6,0.8]','(0.8,1]')
    meltedMDistMatrix <- reshape2::melt(methylationDistMatrix)
    return(meltedMDistMatrix)
}
