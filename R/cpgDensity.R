#'Provides Coverage metrics in the repeat masker region
#'@param bs bsseq object
#'@param organism scientific name of the organism of interest, i.e. Mus musculus or Homo sapiens
#'@param windowLength Length of the window to calculate the density
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'library(BSgenome.Mmusculus.UCSC.mm10)
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'cpgDensity(bs,Mmusculus,1000)
#'@import BSgenome
#'@export


cpgDensity<-function(bs,organism,windowLength=1000){
    cov<-bsseq::getCoverage(bs)
    cpgd<-Repitools::cpgDensityCalc(GenomicRanges::granges(bs),
                                  organism,window = windowLength)
    if (max(cpgd)>50){
        cpgdBin<-cut(cpgd,c(seq(0,50,5),max(cpgd)))
    }else{
        cpgdBin<-cut(cpgd,c(seq(0,20),max(cpgd)))
    }
    cpgdCov <- by(cov>0, cpgdBin, colSums)
    cpgdCov <- do.call("rbind", cpgdCov)
    return(cpgdCov)
}
