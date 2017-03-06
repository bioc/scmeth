#' CpG covergae in each chromosome
#'
#'Provides Coverage metrics for each sample by the chromosome
#'@param bs bsseq object
#'@return matrix of chromsome covergae with
#'column and rows indicating the samples and the chromosome respectively
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
