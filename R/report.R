#' Generates an inclusive report on methylation analysis
#'
#' This function uses most of the functions in this package to
#' generate a report for the user
#'@param bsObj bsseq object
#'@param outdirectory name of the output directory where the
#'report will be saved
#'@param organism scientific name of the organism of interest,
#'i.e. Mus musculus or Homo sapiens
#'@param genome reference alignment, i.e. mm10 or hg38
#'the report will have graphics on read information
#'@return Report will be an html file
#'@examples
#'library(BSgenome.Mmusculus.UCSC.mm10)
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'report(bs,'~',Mmusculus,'mm10')
#'@import knitr
#'@export
#

report <- function(bsObj,outdirectory,organism,genome) {
    RmdFile<-system.file(".",'qcReport.Rmd',package="scmeth")
    rmarkdown::render(RmdFile,params=list(outdir=outdirectory,samples=bsObj
                        ,organism=organism,genome=genome)
                        ,output_file=paste0(outdirectory,"/qcReport.html"))

}


