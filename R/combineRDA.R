#' Combine multiple bsseq rda object
#'
#'
#'Combine multiple bsseq rda object into one bs object
#'Preferably for a library pool of cells
#'
#'@param Takes bs object, name of the organism and reference genome
#'@return Data frame with sample name and coverage in repeat masker regions
#'@example
#'@export


combineRDA<-function(rdaList){
  message("Reading bsseq objects")
  #bsList<-lapply(rdaList,function(file){
  #  load(file)
  #  bs
  #})

  message("Combining bs object")

  bs<-Reduce(combine,bsList)
  bs@parameters<-bsList[[1]]@parameters

  #pData<-pData(bs)
  #pData$col<-as.numeric(factor(pData$type))
  #pData(bs)<-pData

  #message("Saving bs to ")
  #save(bs, file=outfile)

  return(bs)

}
