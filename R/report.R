#' Generates an inclusive report on methylation analysis
#'
#'Takes bs object, output directory, name of the organism, reference genome
#'@param list of inputs
#'@return A report will be the output based on rmarkdown file
#'@example
#'@export
#

report <- function(bsObj,outdirectory,organism,genome,meta = NULL, cacheable = NA) {
  #bsseqObject<-Sys.getenv(rdaFile)
  #if (is.null(outdirectory)){
  #  outdirectory=getwd()
  #}

  rmarkdown::render("./qcReport.Rmd",params=list(outdir=outdirectory,samples=bsObj,organism=organism,genome=genome),output_file=paste0(outdirectory,"/qcReport.html"))

  #return(23)
}


