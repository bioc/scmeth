#' Provide graphics for read information
#'
#'Plot the mapped and unmapped reads
#'@param bs bsseq object
#'@return Plot showing the mapped and unmapped read information for each cell
#'@examples
#'load(system.file("extdata",'bsObject.rda',package='scmeth'))
#'readmetrics(bs)
#'@export
#'@importFrom utils read.delim


readmetrics<-function(bs){
    phenotypicData<-Biobase::pData(bs)
    pDcolNames<-colnames(phenotypicData)
    if (('total_reads' %in% pDcolNames & 'uniquely_aligned_reads' %in% pDcolNames)){
        dat<-data.frame(sample=rownames(phenotypicData),
                        total=phenotypicData$total_reads,
                        mapped=phenotypicData$uniquely_aligned_reads)
        # Attempt to convert to numeric (from factor or string) if necessary
        if (!is.numeric(dat$mapped)) dat$mapped <- as.numeric(as.character(dat$mapped))
        if (!is.numeric(dat$total)) dat$total <- as.numeric(as.character(dat$total))
        dat$unmapped<-dat$total-dat$mapped
        o<-order(dat$total,dat$sample)
        sampleOrder<-dat$sample[o]
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

    }else{
        warning('Read information not provided in the pData slot')
    }
}


