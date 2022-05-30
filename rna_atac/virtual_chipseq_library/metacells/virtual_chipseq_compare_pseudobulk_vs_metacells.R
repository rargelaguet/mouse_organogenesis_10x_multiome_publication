# here::i_am("rna_atac/virtual_chipseq_library/virtual_chipseq_plot_stats.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

opts$TFs <- c("TAL1", "GATA1", "RUNX1", "FOXA2", "GATA4", "CDX2","NKX2-5","TBX5", "SOX10")
# opts$TFs <- c("TAL1", "GATA1", "FOXA2")

# I/O
io$virtual_chip_pseudobulk.dir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s",opts$motif_annotation))
io$virtual_chip_metacells.dir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/%s",opts$motif_annotation))
# io$motifmatcher_chip.se <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/motifmatchr_virtual_chip.rds")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/metacells/%s/pdf",opts$motif_annotation)); dir.create(io$outdir, showWarnings = F)

###################################
## Load virtual ChIP-seq library ##
###################################

# pseudobulk
virtual_chip_pseudobulk.mtx <- readRDS(file.path(io$virtual_chip_pseudobulk.dir,"virtual_chip_matrix.rds"))

# metacells
virtual_chip_metacells.mtx <- readRDS(file.path(io$virtual_chip_metacells.dir,"virtual_chip_matrix.rds"))

################
## Parse data ##
################

stopifnot(opts$TFs%in%colnames(virtual_chip_pseudobulk.mtx))
stopifnot(opts$TFs%in%colnames(virtual_chip_metacells.mtx))

virtual_chip_pseudobulk.mtx <- virtual_chip_pseudobulk.mtx[,opts$TFs]
virtual_chip_metacells.mtx <- virtual_chip_metacells.mtx[,opts$TFs]

stopifnot(rownames(virtual_chip_metacells.mtx)==rownames(virtual_chip_pseudobulk.mtx))

########################################
## Scatterplot of distribution scores ##
########################################

i <- opts$TFs[1]

to.plot <- data.table(metacells = virtual_chip_metacells.mtx[,i], pseudobulk = virtual_chip_pseudobulk.mtx[,i]) %>%
  .[abs(pseudobulk)>0.01]

ggscatter(to.plot, x="pseudobulk", y="metacells", size=0.5) +
  geom_abline(slope=1, intercept=0) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  coord_cartesian(ylim=c(-0.75,0.75), xlim=c(-0.75,0.75)) +
  labs(x="Pseudoblk score", y="Metacells score") +
  theme(
    axis.text = element_text(size=rel(0.75))
  )

####################################
## Compare distribution of scores ##
####################################

to.plot <- rbind(
  data.table(value = virtual_chip_metacells.mtx[,i], class="metacells"),
  data.table(value = virtual_chip_pseudobulk.mtx[,i], class="pseudobulk")
) %>% .[value!=0]

gghistogram(to.plot, x="value", fill="class", bins=100) +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size=rel(0.75), color="black")
  )

###################################################
## Plot score thresholds vs number binding sites ##
###################################################

# seq.ranges <- seq(0.25,1,by=0.01)
# 
# to.plot <- seq.ranges %>% map(function(j) {
#   data.table(
#     tf = opts$TFs,
#     min_score = j,
#     log2_N = log2(colSums(virtual_chip_filt.mtx>=j)+1)
#   )
# }) %>% rbindlist
# 
# p.lineplot <- ggplot(to.plot[log2_N>0], aes_string(x="min_score", y="log2_N", color="tf")) +
#   geom_line(size=1) +
#   labs(y="Number of predicted binding sites (log2)", x="Minimum in silico binding score") +
#   scale_color_brewer(palette="Dark2") +
#   # coord_cartesian(xlim=c(0,0.9)) +
#   theme_classic() +
#   theme(
#     axis.text = element_text(color="black"),
#     # legend.position = "right",
#     legend.position = c(.9,.65),
#     legend.title = element_blank()
#   )
# 
# pdf(file.path(io$outdir,"lineplot_number_binding_sites_per_tf.pdf"), width=6, height=5)
# print(p.lineplot)
# dev.off()
