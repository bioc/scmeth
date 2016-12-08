#' square a number
#'
#'Take in a numeric value and squares it
#'@param x A numeric input
#'@return The sqaure of the output
#'@example
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
