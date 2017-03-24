#'Provides Coverage metrics in the repeat masker region
#'@param bs bsseq object
#'@param organism scientific name of the organism of interest, i.e. Mus musculus or Homo sapiens
#'@param genome reference alignment, i.e. mm10 or hg38
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'repMask(bsseqObject,'Mus musculus','mm10')
#'@export


cpgDensity<-function(bs,organism,windowLength=100){
  if (organism == 'Mus musculus'){

  library(BSgenome.Mmusculus.UCSC.mm10)
  cov<-getCoverage(bs)
  cpgd<-Repitools::cpgDensityCalc(granges(bs),Mmusculus,window = windowLength)
  cpgdBin<-cut(cpgd,c(seq(0,20),max(cpgd)))
  cpgdCov <- by(cov>0, cpgdBin, colSums)
  cpgdCov <- do.call("rbind", cpgdCov)

  a}else if (organism == 'Homo sapiens'){
    library(BSgenome.Hsapiens.UCSC.hg19)
    cpgd<-Repitools::cpgDensityCalc(granges(bs),Hsapeins,window = windowLength)
    cpgdBin<-cut(cpgd,c(seq(0,20),max(cpgd)))
    cpgdCov <- by(cov>0, cpgdBin, colSums)
    cpgdCov <- do.call("rbind", cpgdCov)

  }
  return(cpgdCov)
}
