#' square a number
#'
#'Take in a numeric value and squares it
#'@param x A numeric input
#'@return The sqaure of the output
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

