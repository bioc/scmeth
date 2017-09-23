#' Provide graphics for read information
#'
#'Plot the mapped and unmapped reads
#'@param dir directory that contains bismark output
#'@return Plot showing the mapped and unmapped read information for each cell
#'@examples
#'directory<-system.file("extdata/bismark_data",package='scmeth')
#'bs<-SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
#'readmetrics(bs)
#'@export
#'@importFrom utils read.delim
#'@importFrom stats relevel

readmetrics<-function(bs){
    phenotypicData<-Biobase::pData(bs)
    dat<-data.frame(sample=phenotypicData$cell_id,
                    total=as.vector(phenotypicData$total_reads),
                    mapped=as.vector(phenotypicData$uniquely_aligned_reads+
                                phenotypicData$non_uniquely_aligned_reads))
    dat$unmapped<-dat$total-dat$mapped
    o<-order(dat$total,dat$sample)
    sampleOrder<-dat$sample[o]
    dat$sample<-factor(dat$sample,levels=sampleOrder)

    m<-reshape2::melt(dat[,c("sample","mapped","unmapped")],
                      id.vars="sample",variable.name="Mapping_status")
    m$sample<-factor(m$sample,levels=sampleOrder)
    m$Mapping_status <- relevel(m$Mapping_status, ref="unmapped")
    g<-ggplot2::ggplot(m,ggplot2::aes_string('sample','value',
                                             fill='Mapping_status'))
    g<-g+ggplot2::geom_bar(stat="identity")+ggplot2::coord_flip()
    g<-g+ggplot2::scale_y_continuous(name="Number of reads")
    g<-g+ggplot2::xlab("samples")
    g<-g+ggplot2::ggtitle("Read mapping stats")
    g<-g+ggplot2::theme_bw()
    return(g)

}


