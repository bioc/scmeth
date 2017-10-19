#'Provides Coverage metrics in the repeat masker region
#'@param bs bsseq object
#'@param organism scientific name of the organism of interest,
#'i.e. Mus musculus or Homo sapiens
#'@param windowLength Length of the window to calculate the density
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'library(BSgenome.Mmusculus.UCSC.mm10)
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'cpgDensity(bs,Mmusculus,1000)
#'@import BSgenome
#'@importFrom bsseq getCoverage
#'@export



cpgDensity<-function(bs,organism,windowLength=1000){
    cov<-bsseq::getCoverage(bs)
    cpgd<-Repitools::cpgDensityCalc(GenomicRanges::granges(bs),
                                organism,window = windowLength)
    cpgdBin<-cut(cpgd,seq(0,max(cpgd,5)))

    cpgdCov <- by(cov>0, cpgd, colSums)
    cpgdCov <- do.call("rbind", cpgdCov)
    return(cpgdCov)
}
