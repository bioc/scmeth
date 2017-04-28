#' Combine multiple bsseq rda object
#'
#'
#'Combine multiple bsseq rda object into one bs object
#'Preferably for a library pool of cells
#'
#'@param rdaList list of bsseq rda objects
#'@return Combined one rda object for all the cells
#'@examples
#'bsObject1<-createRDA(system.file("extdata",'sc-RRBS_zyg_01_chr1_CpG.bedGraph',package='scmeth'))
#'bsObject2<-createRDA(system.file("extdata",'sc-RRBS_zyg_02_chr1_CpG.bedGraph',package='scmeth'))
#'combineRDA(c(bsObject1, bsObject2))
#'@export


combineRDA<-function(rdaList){
  message("Reading bsseq objects")

  message("Combining bs object")
  bs<-Reduce(BiocGenerics::combine,rdaList)
  bs@parameters<-rdaList[[1]]@parameters

  return(bs)

}
