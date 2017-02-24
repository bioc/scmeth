#' Coverage statistics based on chromosome
#'
#'Provides Coverage metrics for the sample by each chromosome
#'@param takes bs object as input
#'@return matrix of numbers as output with
#'column and rows indicating the samples and the chromosome
#'@examples
#'chromosomeCoverage(bsObject)
#'@export


chromosomeCoverage <- function(bs) {
  covMatrix<-getCoverage(bs)
  Granges<-granges(bs)
  allChr <- as.vector(seqnames(Granges))
  Granges<-keepStandardChromosomes(Granges)
  standardChr <- as.vector(seqnames(Granges))
  standardChrInd<- allChr %in% standardChr
  covMatrix<-covMatrix[standardChrInd,]
  chrCov <- by(covMatrix>0, standardChr, colSums)
  chrCov <- do.call("rbind", chrCov)
  return(chrCov)
}
