#' Repeat masker coverage statistics
#'
#'Provides Coverage metrics based on the repeat masker region
#'@param Takes bs object, name of the organism and reference genome
#'@return Data frame with sample name and coverage in repeat masker regions
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
