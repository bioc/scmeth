#' Combine multiple bsseq rda object
#'
#'
#'Combine multiple bsseq rda object into one bs object
#'Preferably for a library pool of cells
#'
#'@param Takes bs object, name of the organism and reference genome
#'@return Data frame with sample name and coverage in repeat masker regions
#'@examples
#'combineRDA(c(bsObject1, bsObject2))
#'@export


combineRDA<-function(rdaList){
  message("Reading bsseq objects")


  message("Combining bs object")

  bs<-Reduce(combine,rdaList)
  bs@parameters<-rdaList[[1]]@parameters

  #pData<-pData(bs)
  #pData$col<-as.numeric(factor(pData$type))
  #pData(bs)<-pData

  #message("Saving bs to ")
  #save(bs, file=outfile)

  return(bs)

}
