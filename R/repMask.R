#' square a number
#'
#'Take in a numeric value and squares it
#'@param x A numeric input
#'@return The sqaure of the output
#'@example
#'@export


repMask<-function(bs,organism,genome){
  hub <- AnnotationHub()
  repeatGr <- hub[[names(query(hub, c("rmsk", organism, genome)))]]
  rep <- countOverlaps(bs, repeatGr)>0
  cov=getCoverage(bs)
  covDf <- data.frame(sample=sampleNames(bs), coveredCpgs=colSums(cov[!rep,]>=1))
  return(covDf)

}
