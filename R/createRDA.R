#' Generate bsseq rda object
#'
#'
#'Generate bsseq rda object from the methyldackel/pileometh output
#'
#'@param CpG_file CpG.bedGraph file output from methyldackel/pileometh
#'@param readmetric_file file containing read information
#'@param bsConv_file file containing bisulfite conversion information
#'@return bsseq object with coverga and methylation information for the sample
#'@examples
#'CpGBedGraphFile<-system.file("extdata",'sc-RRBS_zyg_01_chr1_CpG.bedGraph',
#'package='scmeth')
#'readMetricsFile<-system.file("extdata",'sc-RRBS-zygote_01.read_metrics.txt',
#'package='scmeth')
#'bsConversionFile<-system.file("extdata",'sc-RRBS-zygote_01.bsConv.txt',
#'package='scmeth')
#'createRDA(CpGBedGraphFile,readMetricsFile,bsConversionFile)
#'@export


createRDA<-function(CpG_file, readmetric_file,bsConv_file){
    tab<-read.delim(CpG_file,sep='\t')
    message('Creating BSseq object with', nrow(tab),'loci.')
    sample<-sub("_CpG.bedGraph","",basename(CpG_file))
    m<-tab[,5] # methylated count
    um<-tab[,6] # unmethylated count
    cov<-m+um

    # Get the phenotypic data
    readInfo<-read.table(readmetric_file,sep=' ')
    bsconvInfo<-read.table(bsConv_file,sep=' ', header=TRUE)
    pd<-data.frame(row.names=sample,totalReads=as.numeric(readInfo[2]),
                    mappedReads=as.numeric(readInfo[3]),
                    bsconversion=as.numeric(bsconvInfo[1,1]),
                    stringsAsFactors=FALSE)

    message("generating bs object")
    bs<-bsseq::BSseq(chr=tab[,1], pos=tab[,3], M= matrix(m),
                    Cov=matrix(cov),sampleNames=sample,pData=pd)
    message("Done.")

    return(bs)
}
