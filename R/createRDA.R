#' Generate bsseq rda object
#'
#'
#'Generate bsseq rda object from the methyldackel/pileometh output
#'
#'@param file CpG.bedGraph file output from methyldackel/pileometh
#'@return bsseq object with coverga and methylation information for the sample
#'@examples
#'createRDA(system.file("extdata",'sc-RRBS_zyg_01_chr1_CpG.bedGraph',package='scmeth'))
#'@export


createRDA<-function(file){
  tab<-read.delim(file,sep='\t')
  message('Creating BSseq object with', nrow(tab),'loci.')
  sample<-sub("_CpG.bedGraph","",basename(file))
  m<-tab[,5] # methylated count
  um<-tab[,6] # unmethylated count
  cov<-m+um

  message("generating bs object")
  bs<-bsseq::BSseq(chr=tab[,1], pos=tab[,3], M= matrix(m), Cov=matrix(cov),sampleNames=sample)
  message("Done.")

  return(bs)

}
