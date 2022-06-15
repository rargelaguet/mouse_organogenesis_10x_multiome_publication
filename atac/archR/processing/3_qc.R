here::i_am("atac/archR/processing/3_qc.R")

source(here::here("settings.R"))

suppressPackageStartupMessages(library(ArchR))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',    type="character",    help='ArchR directory')
p$add_argument('--metadata',    type="character",    help='metadata file')
p$add_argument('--outdir',     type="character",    help='Output directory')
p$add_argument('--min_tss_enrichment',     type="integer",    default=8,   help='Minimum TSS enrichment')
p$add_argument('--max_tss_enrichment',     type="integer",    default=40,   help='Minimum TSS enrichment')
p$add_argument('--min_number_fragments',     type="integer",    default=3000,    help='Maximum number of ATAC fragments')
p$add_argument('--max_number_fragments',     type="integer",    default=300000,    help='Maximum number of ATAC fragments')
p$add_argument('--max_blacklist_ratio',     type="double",    default=0.05,    help='Maximum Blacklist Ratio')
p$add_argument('--threads',     type="integer",    default=1,    help='Number of threads')

args <- p$parse_args(commandArgs(TRUE))


#####################
## Define settings ##
#####################

# START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$metadata <- file.path(io$basedir,"processed/atac/archR/sample_metadata_after_archR.txt.gz")
# args$min_tss_enrichment <- 9
# args$max_tss_enrichment <- 35
# args$min_number_fragments <- 3500
# args$max_number_fragments <- 250000
# args$max_blacklist_ratio <- 0.05
# args$threads <- 2
# args$outdir <- file.path(io$basedir,"results/atac/archR/qc")
## END TEST ##

# Options
opts$chr <- paste0("chr",1:3)
opts$test <- TRUE

########################
## Load cell metadata ##
########################

sample_metadata <- fread(args$metadata)

# temporary
sample_metadata[is.na(stage),stage:=strsplit(sample,"_") %>% map_chr(1)]

########################
## Load ArchR project ##
########################

# source(here::here("atac/archR/load_archR_project.R"))

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = args$threads)

ArchRProject <- loadArchRProject(args$archr_directory)[sample_metadata$cell]

##################
## Subset ArchR ##
##################

if (opts$test) {
  cells.to.use <- split(ArchRProject$cellNames,ArchRProject$sample) %>% map(~ head(.,n=100)) %>% unlist
  ArchRProject <- ArchRProject[cells.to.use,]
}

# Subset chr for faster computations
tss.granges <- getTSS(ArchRProject)
tss.granges <- tss.granges[seqnames(tss.granges)%in%opts$chr]

#########################
## Plot TSS Enrichment ##
#########################

data_tss.dt <- opts$samples %>% map(function(i) {
  plotTSSEnrichment(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE,
    TSS = tss.granges
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample") %>% 
  melt(id.vars=c("sample","x"))
fwrite(data_tss.dt, sprintf("%s/qc_TSSenrichment.txt.gz",args$outdir))

# data_tss.dt <- fread(sprintf("%s/qc_TSSenrichment.txt.gz",args$outdir)) %>%
#   .[,.(value=mean(value)), by = c("sample","x","variable")]

to_plot_tss.dt <- data_tss.dt %>% 
  merge(unique(sample_metadata[,c("sample","stage")])) %>%
  .[,.(value=mean(value)), by = c("sample","stage","x","variable")] %>%
  .[variable=="normValue"]

p <- ggline(to_plot_tss.dt, x="x", y="value", plot_type="l") +
  facet_wrap(~stage, scales="fixed") +
  # scale_colour_manual(values=opts$stage.colors) +
  # scale_x_continuous(breaks=seq(-2000,2000,1000)) +
  labs(x="Distance from TSS (bp)", y="TSS enrichment (normalised)") +
  theme(
    axis.text.y = element_text(size=rel(0.65), color="black"),
    axis.text.x = element_text(size=rel(0.6), color="black"),
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title = element_text(size=rel(0.75), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  )

pdf(file.path(args$outdir,"qc_TSSenrichment.pdf"), width=6, height=5)
print(p)
dev.off()

#####################################
## Plot Fragment size distribution ##
#####################################

# to.plot.fragmentsize <- plotFragmentSizes(ArchRProject, groupBy = "Sample", returnDF=T) %>% 
#   as.data.table %>% setnames("group","sample")
data_fragmentsize.dt <- opts$samples %>% map(function(i) {
  plotFragmentSizes(
    ArchRProj = ArchRProject[ArchRProject$Sample==i,], 
    groupBy = "Sample", 
    returnDF = TRUE
  ) %>% as.data.table %>% return
}) %>% rbindlist %>% 
  setnames("group","sample")

fwrite(data_fragmentsize.dt, sprintf("%s/qc_FragmentSizeDistribution.txt.gz",args$outdir))

# data_fragmentsize.dt <- fread(sprintf("%s/qc_FragmentSizeDistribution.txt.gz",args$outdir)) %>%
#   .[,.(fragmentPercent=mean(fragmentPercent)), by = c("sample","fragmentSize")]

to_plot_fragmentsize.dt <- data_fragmentsize.dt %>% 
  merge(unique(sample_metadata[,c("sample","stage")])) %>%
  .[,.(fragmentPercent=mean(fragmentPercent)), by = c("sample","stage","fragmentSize")] %>%
  .[fragmentSize<=400]

p <- ggline(to_plot_fragmentsize.dt, x="fragmentSize", y="fragmentPercent", plot_type="l") +
  facet_wrap(~stage, scales="fixed") +
  # scale_x_continuous(breaks=seq(125,750,125)) +
  # scale_colour_manual(values=opts$sample.colors) +
  labs(x="Fragment Size (bp)", y="Percentage of fragments (%)") +
  theme(
    axis.text = element_text(size=rel(0.55)),
    axis.title = element_text(size=rel(0.75)),
    legend.position = "none",
    legend.title = element_blank()
  )

pdf(file.path(args$outdir,"qc_FragmentSizeDistribution.pdf"), width=8, height=4)
print(p)
dev.off()

##################################
## Plot histogram of QC metrics ##
##################################

to.plot <- sample_metadata %>%
  .[!is.na(nFrags_atac)] %>%
  .[,log_nFrags_atac:=log10(nFrags_atac)] %>%
  melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags_atac","BlacklistRatio_atac"))
  # melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags"))

# tmp <- data.table(
#   variable = c("TSSEnrichment_atac", "nFrags_atac", "BlacklistRatio_atac"),
#   value = c(args$min_tss_enrichment, args$min_number_fragments, args$max_blacklist_ratio)
# )
tmp <- data.table(
  variable = c("TSSEnrichment_atac", "TSSEnrichment_atac", "log_nFrags_atac", "log_nFrags_atac","BlacklistRatio_atac"),
  value = c(args$min_tss_enrichment, args$max_tss_enrichment, log10(args$min_number_fragments), log10(args$max_number_fragments),args$max_blacklist_ratio)
)
# p <- gghistogram(to.plot, x="value", fill="sample", bins=50) +
p <- gghistogram(to.plot, x="value", y="..density..", bins=70, fill="sample") +
  geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
  facet_wrap(~variable, scales="free") +
  theme(
    axis.text =  element_text(size=rel(0.8)),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size=rel(0.75))
  )
pdf(sprintf("%s/qc_metrics_histogram.pdf",args$outdir), width=8, height=5)
print(p)
dev.off()


#############
## Call QC ##
#############

sample_metadata %>%
  .[,pass_atacQC:=TSSEnrichment_atac>=args$min_tss_enrichment & TSSEnrichment_atac<=args$max_tss_enrichment & 
      nFrags_atac>=args$min_number_fragments & nFrags_atac<=args$max_number_fragments & 
      BlacklistRatio_atac<=args$max_blacklist_ratio] %>%
  .[is.na(pass_atacQC),pass_atacQC:=FALSE]

print(sample_metadata[,mean(pass_atacQC,na.rm=T),by="sample"])
# print(sample_metadata[,mean(is.na(nFrags_atac)),by="sample"])

# Filter low quality cells that did not pass RNA QC
sample_metadata[pass_rnaQC==FALSE & nFrags_atac<=1e4,pass_atacQC:=FALSE]

# Save
fwrite(sample_metadata, file.path(args$outdir,"sample_metadata_after_qc.txt.gz"), quote=F, na="NA", sep="\t")


###########################################
## Plot QC statistics after QC filtering ##
###########################################

# Barplot of the fraction of cells that pass QC for each sample

to.plot <- sample_metadata %>%
  .[,mean(pass_atacQC,na.rm=T),by=c("sample","stage")]

p <- ggbarplot(to.plot, x="sample", y="V1", fill="stage") +
  scale_fill_manual(values=opts$stage.colors) +
  labs(x="", y="Fraction of cells that pass ATAC QC") +
  coord_cartesian(ylim=c(0,1)) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(colour="black",size=rel(0.8)),
    axis.text.x = element_text(colour="black",size=rel(0.65), angle=20, hjust=1, vjust=1),
  )

pdf(sprintf("%s/qc_metrics_barplot.pdf",args$outdir), width=6, height=5)
print(p)
dev.off()


# Boxplots of QC metrics
to.plot <- sample_metadata %>%
  .[pass_atacQC==TRUE] %>%
  .[nFrags_atac<=150000 & TSSEnrichment_atac<=27] %>% # remove massive outliers for plotting
  .[,log_nFrags_atac:=log10(nFrags_atac)] %>%
  # melt(id.vars=c("sample","cell"), measure.vars=c("TSSEnrichment_atac","log_nFrags","BlacklistRatio_atac"))
  melt(id.vars=c("cell","sample","stage"), measure.vars=c("TSSEnrichment_atac","log_nFrags_atac"))

facet.labels <- c("log_nFrags_atac" = "Num. of fragments (log10)", "TSSEnrichment_atac" = "TSS enrichment")

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

pdf(sprintf("%s/qc_metrics_boxplot.pdf",args$outdir), width=9, height=5)
print(p)
dev.off()
