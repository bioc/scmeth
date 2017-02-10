#' Generate bsseq rda object
#'
#'
#'Generate bsseq rda object from the pileometh output
#'
#'@param Takes bs object, name of the organism and reference genome
#'@return Data frame with sample name and coverage in repeat masker regions
#'@example
#'@export


createRDA<-function(file){
  tab<-read.delim(file,sep='\t')
  message('Creating BSseq object with', nrow(tab),'loci.')
  sample<-sub("_CpG.bedGraph","",basename(file))
  m<-tab[,5] # methylated count
  um<-tab[,6] # unmethylated count
  cov<-m+um

  message("generating bs object")
  bs<-BSseq(chr=tab[,1], pos=tab[,3], M= matrix(m), Cov=matrix(cov),sampleNames=sample)
  message("Done.")

  return(bs)

}
