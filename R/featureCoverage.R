#' Coverage based on the genomic feature
#'
#'Provides Coverage metrics for the sample by each genomic features provided by the user
#'@param bs bsseq object
#'@param features list of genomic features, i.e. genes_exons, genes_introns, cpg_islands, cpg_shelves
#'Names are based on the annotatr packages, so all the features provided by the annotatr
#'package will be supported in this function
#'@param genome reference alignment, i.e. mm10 or hg38
#'@return a data frame with genomic feature names and the number of CpG covered in each feature
#'
#'@examples
#'featureCoverage(bsObject,c('cpg_islands','genes_exons'),'mm10')
#'@export


featureCoverage <-function(bs,features,genome){
  if (!requireNamespace('annotatr',quietly=TRUE)){
    stop("Pkg needed for this function to work. Please install it")
  }
  library(annotatr)
  annotationFeatures<-c()
  for (i in features){
    annotationFeatures<-c(paste0(genome,'_',i),annotationFeatures)

  }
  annots_gr = annotatr::build_annotations(genome = genome, annotations = annotationFeatures)


  # Intersect the regions with the reference annotations
  dm_annotated = annotatr::annotate_regions(
    regions = granges(bs),
    annotations = annots_gr,
    ignore.strand = TRUE,
    quiet = TRUE)
  sumAnnot<-annotatr::summarize_annotations(dm_annotated,quiet=TRUE)
  sumAnnotDf<-as.data.frame(sumAnnot)
  return(sumAnnotDf)

}
