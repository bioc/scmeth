#' Generates an inclusive report on methylation analysis
#'
#' This function uses most of the functions in this package to generate a report for the user
#'@param bsObj bsseq object
#'@param outdirectory name of the output directory where the report will be saved
#'@param organism scientific name of the organism of interest, i.e. Mus musculus or Homo sapiens
#'@param genome reference alignment, i.e. mm10 or hg38
#'@param readData read information file on the samples. This is an optional argument and if given
#'the report will have graphics on read information
#'@return Report will be an html file
#'@examples
#'report(bsseqObject,outputDirectory,'Mus musculus','mm10')
#
#'@export
#

report <- function(bsObj,outdirectory,organism,genome,readData = NULL) {
  #bsseqObject<-Sys.getenv(rdaFile)
  #if (is.null(outdirectory)){
  #  outdirectory=getwd()
  #}

  #RmdFile<-file.path(PROJHOME,'scmeth/R','qcReport.Rmd')
  RmdFile<-system.file(".",'qcReport.Rmd',package="scmeth")

  rmarkdown::render(RmdFile,params=list(outdir=outdirectory,samples=bsObj,organism=organism,genome=genome,reads=readData),output_file=paste0(outdirectory,"/qcReport.html"))

  #return(23)
}


