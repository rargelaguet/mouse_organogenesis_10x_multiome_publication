here::i_am("rna/processing/2_QC.R")

source(here::here("settings.R"))

#####################
## Define arguments ##
#####################

p <- ArgumentParser(description='')
p$add_argument('--metadata',       type="character",                    help='Metadata')
p$add_argument('--outputdir',       type="character",                    help='Output directory')
p$add_argument('--min_nFeature_RNA',       type="integer",                    help='Minimum number of expressed genes')
p$add_argument('--max_nFeature_RNA',       type="integer",                    help='Maximum number of expressed genes')
p$add_argument('--mitochondrial_percent_RNA',       type="integer",                    help='Maximum percentage of mitochondrial reads')
p$add_argument('--ribosomal_percent_RNA',       type="integer",                    help='Maximum percentage of ribosomal reads')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$samples <- opts$samples
# # args$metadata <- paste0(io$basedir,"/processed_new/rna_new/metadata.txt.gz")
# args$metadata <- io$metadata
# args$min_nFeature_RNA <- 2000
# args$max_nFeature_RNA <- 10000
# args$mitochondrial_percent_RNA <- 40
# args$ribosomal_percent_RNA <- 20
# args$outputdir <- paste0(io$basedir,"/results_new/rna/qc")
## END TEST ##

# Sanity checks
stopifnot(args$samples%in%opts$samples)

###############
## Load data ##
###############

metadata <- fread(args$metadata) %>% 
    .[sample%in%args$samples] %>%
    # .[,pass_rnaQC:=nFeature_RNA>args$min_nFeature_RNA & nCount_RNA>2**args$log_nCount_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA]
    .[,pass_rnaQC:=nFeature_RNA<=args$max_nFeature_RNA & nFeature_RNA>=args$min_nFeature_RNA & mitochondrial_percent_RNA<args$mitochondrial_percent_RNA & ribosomal_percent_RNA<args$ribosomal_percent_RNA]

# temporary
metadata[stage=="E8.55",stage:="E8.5"]
table(metadata$pass_rnaQC)

# Ad hoc editing because of outlier cluster of high rRNA in the E7.75 sample
metadata[sample=="E7.75_rep1" & ribosomal_percent_RNA>=8,pass_rnaQC:=FALSE]

#####################
## Plot QC metrics ##
#####################

to.plot <- metadata %>% .[pass_rnaQC==TRUE] %>%
    .[nFeature_RNA<=8000 & mitochondrial_percent_RNA<=60 & ribosomal_percent_RNA<=16] %>% # remove massive outliers for plotting
    # .[,log_nFeature_RNA:=log10(nFeature_RNA)] %>%
    melt(id.vars=c("sample","cell","stage"), measure.vars=c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))
    # melt(id.vars=c("sample","cell"), measure.vars=c("nFeature_RNA"))

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mitochondrial_percent_RNA" = "Mitochondrial %", "ribosomal_percent_RNA" = "Ribosomal %")
    
## Box plot 

p <- ggplot(to.plot, aes_string(x="sample", y="value", fill="stage")) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    scale_fill_manual(values=opts$stage.colors) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir), width=9, height=5)
# pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outputdir))
print(p)
dev.off()

## histogram 

tmp <- data.table(
    variable = c("nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA"),
    value = c(args$min_nFeature_RNA, args$mitochondrial_percent_RNA, args$ribosomal_percent_RNA)
)
# tmp <- data.table(
#     variable = c("nFeature_RNA"),
#     value = c(args$min_nFeature_RNA)
# )

p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
# p <- ggdensity(to.plot, x="value", fill="sample") +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free", nrow=1) +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.75))
    )
    
pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outputdir), width=13, height=6)
print(p)
dev.off()


#########################################################
## Plot fraction of cells that pass QC for each sample ##
#########################################################

to.plot <- metadata %>%
    .[,mean(pass_rnaQC,na.rm=T),by=c("sample","stage")]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="stage") +
    scale_fill_manual(values=opts$stage.colors) +
    labs(x="", y="Fraction of cells that pass QC (RNA)") +
    # facet_wrap(~stage)
    theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
    )

pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outputdir), width=6, height=5)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(metadata, paste0(args$outputdir,"/sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")

