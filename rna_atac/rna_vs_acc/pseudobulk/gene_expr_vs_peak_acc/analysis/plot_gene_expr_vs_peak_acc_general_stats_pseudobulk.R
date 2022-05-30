here::i_am("rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_peak_acc/analysis/plot_gene_expr_vs_peak_acc_general_stats_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$gene_expr_vs_peak_acc <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_peak_acc/cor_gene_expr_vs_peak_acc_pseudobulk.txt.gz")
io$atac_peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
io$peak2gene <- file.path(io$basedir,"results/atac/archR/peak_calling/peaks2genes/peaks2genes_all.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna_atac/rna_vs_acc/pseudobulk/gene_expr_vs_peak_acc/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$max_genomic_distance <- 5e4

###############
## Load data ##
###############

gene_expr_vs_peak_acc.dt <- fread(io$gene_expr_vs_peak_acc) %>%
  .[is.na(cor), c("cor","pvalue"):=list(0,1)] %>% .[,dist:=NULL]

# Filter genes
genes.to.use <- grep("*Rik|^Gm|^Mt-|^Rps|^Rpl|^Olfr",unique(gene_expr_vs_peak_acc.dt$gene), invert=T,value=T)
gene_expr_vs_peak_acc.dt <- gene_expr_vs_peak_acc.dt[gene%in%genes.to.use]

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$atac_peak_metadata) %>%
  .[,peak:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("peak","peakType","distToGeneStart","nearestGene")]

##########################
## Load peak2gene links ##
##########################

peak2gene.dt <- fread(io$peak2gene) %>% 
  .[gene%in%genes.to.use] %>%
  .[,c("peak","gene","dist")] %>%
  merge(peak_metadata.dt[,c("peak","peakType")])

# Note that this does not contain peaks with no links to genes?
length(unique(peak2gene.dt$peak))

stopifnot(unique(peak2gene.dt$peak)%in%unique(peak_metadata.dt$peak))

###############################################
## Calculate number of peaks linked to genes ##
###############################################

peak_metadata.dt[distToGeneStart<=opts$max_genomic_distance,.N] / nrow(peak_metadata.dt)

# tmp <- peak2gene.dt[,sum(dist<=opts$max_genomic_distance), by="peak"]
# mean(tmp$V1>0)

###########
## Merge ##
###########

gene_expr_vs_peak_acc_filt.dt <- gene_expr_vs_peak_acc.dt %>% 
  merge(peak2gene.dt,by=c("peak","gene")) %>%
  .[,sign:=c("-","+")[as.numeric(cor>0)+1]]

length(unique(gene_expr_vs_peak_acc_filt.dt$peak))

###########################################################################################
## Calculate number of peaks that are correlated with RNA expression of the nearest gene ##
###########################################################################################

opts$max.pval <- 0.05
opts$min.cor <- 0.05

tmp <- gene_expr_vs_peak_acc_filt.dt %>% .[,sum(pvalue<=opts$max.pval & abs(cor)>=opts$min.cor), by=c("peak")] %>% .[,mean(V1>=1)]

to.plot <- gene_expr_vs_peak_acc_filt.dt %>%
  .[,sum(pvalue<=opts$max.pval & abs(cor)>=opts$min.cor), by=c("peak","sign","peakType")] %>%
  .[,100*mean(V1>=1),by=c("sign","peakType")]

p <- ggbarplot(to.plot, x="peakType", fill="sign", y="V1", position=position_dodge(width = 0.75)) +
  labs(x="", y="% of peaks correlated with proximal gene expr.") +
  # scale_fill_brewer(palette="Dark2") +
  # guides(color=F) + 
  theme(
    axis.title = element_text(size=rel(0.80)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.80), color="black"),
    axis.text.y = element_text(size=rel(0.80), color="black")
  )

pdf(file.path(io$outdir,"barplots_number_gene_expr_vs_peak_acc_correlations.pdf"), width=7, height=5)
print(p)
dev.off()

####################################################
## Histogram of ATAC peak acc vs RNA correlations ##
####################################################

to.plot <- gene_expr_vs_peak_acc_filt.dt# %>% .[pvalue<=0.05]

p <- gghistogram(to.plot, x="cor", y="..density..", bins=50, fill="gray70") +
  labs(x="Pearson correlation coefficient\n(ATAC peak vs gene expr.)", y="Density") +
  theme(
    axis.title.x = element_text(size=rel(0.9)),
    axis.title.y = element_text(size=rel(1.0)),
    legend.title = element_blank(),
    axis.text.x = element_text(size=rel(0.8), color="black"),
    axis.text.y = element_text(size=rel(0.8), color="black")
  )


pdf(file.path(io$outdir,"histogram_gene_expr_vs_peak_acc_corr.pdf"), width=7, height=5)
print(p)
dev.off()

##############################
## Number of peaks per gene ##
##############################

tmp <- peak2gene.dt %>% .[dist<=opts$max_genomic_distance] %>% .[,.(npeaks=.N),by=c("gene")]

to.plot <- peak2gene.dt %>% .[dist<=opts$max_genomic_distance] %>% 
  .[,.(npeaks=.N),by=c("gene")] %>%
  .[,.N,by="npeaks"] %>%
  .[npeaks>=2] %>%
  .[npeaks>=200,npeaks:=200] # for viz purposes

p <- ggline(to.plot, x="npeaks", y="N") +
  labs(x="Number of peaks per gene", y="Number of genes") +
  theme(
    axis.title = element_text(size=rel(1.0), color="black"),
    axis.text = element_text(size=rel(0.80), color="black")
  )


pdf(file.path(io$outdir,"lineplot_npeaks_per_gene.pdf"), width=7, height=5)
print(p)
dev.off()
  

to.plot <- peak2gene.dt %>% .[dist<=opts$max_genomic_distance] %>% .[,.(npeaks=.N),by=c("gene")]
p <- gghistogram(to.plot, x="npeaks", y="..density..", bins=50, fill="gray70") +
  scale_x_continuous(limits=c(0,100)) +
  labs(x="Number of peaks per gene", y="Density") +
  theme(
    axis.title = element_text(size=rel(1.0), color="black"),
    axis.text = element_text(size=rel(1.0), color="black")
  )

pdf(file.path(io$outdir,"histogram_npeaks_per_gene.pdf"), width=7, height=5)
print(p)
dev.off()


##############################
## Number of genes per peak ##
##############################

tmp <- peak2gene.dt %>% .[dist<=opts$max_genomic_distance] %>% .[,.(ngenes=.N),by=c("peak")]
mean(tmp$ngenes)

to.plot <- peak2gene.dt %>% .[dist<=opts$max_genomic_distance] %>% 
  .[,.(ngenes=.N),by=c("peak")] %>%
  .[,.N,by="ngenes"] %>%
  .[ngenes<=30]

p <- ggline(to.plot, x="ngenes", y="N") +
  labs(x="Number of genes per peak", y="Number of peaks") +
  theme(
    axis.title = element_text(size=rel(1.0), color="black"),
    axis.text = element_text(size=rel(0.80), color="black")
  )


pdf(file.path(io$outdir,"lineplot_ngenes_per_peak.pdf"), width=7, height=5)
print(p)
dev.off()
