#' Coverage statistics based on the feature
#'
#'Provides Coverage metrics for the sample by each features provided by the user
#'@param takes bs object and list of features as input
#'@return matrix a data frame with feature names and th number of CpG covered
#'
#'@examples
#'featureCoverage(bsObject,c('cpg_island','exon'),'mm10')
#'@export


featureCoverage <-function(bs,features,genome){
  if (!requireNamespace('annotatr',quietly=TRUE)){
    stop("Pkg needed for this function to work. Please install it")
  }
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
