here::i_am("rna_atac/virtual_chipseq_library/pseudobulk/analysis/stats/virtual_chipseq_plot_stats.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# opts$TFs <- c("TAL1", "GATA1", "RUNX1", "FOXA2", "GATA4", "CDX2","NKX2-5","TBX5", "SOX10")
opts$TFs <- c("CDX2", "TAL1", "GATA1", "FOXA2", "GATA4","TBX5")

# I/O
io$virtual_chip.dir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s",opts$motif_annotation))
io$virtual_chip.mtx <- file.path(io$virtual_chip.dir,"virtual_chip_matrix.rds")
# io$motifmatcher_chip.se <- file.path(io$basedir,"results_new/rna_atac/virtual_chipseq/CISBP/motifmatchr_virtual_chip.rds")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/stats",opts$motif_annotation)); dir.create(io$outdir, showWarnings = F)

######################
## Load motifmatchR ##
######################

# io$motifmatcher.se <- sprintf("%s/Annotations/%s-Scores.rds",io$archR.directory,opts$motif_annotation)
# motifmatcher.se <- readRDS(io$motifmatcher.se)

###################################
## Load virtual ChIP-seq library ##
###################################

virtual_chip.dt <- opts$TFs %>% map(function(i) {
# virtual_chip.dt <- list.files(io$virtual_chip.dir,".txt.gz") %>% strsplit("\\.") %>% map_chr(1) %>% map(function(i) {
    fread(sprintf("%s/%s.txt.gz",io$virtual_chip.dir,i)) %>%
    # .[,c("chr","start","end"):=NULL] %>%
    .[,tf:=i] %>%
    return
}) %>% rbindlist

# Load matrix
virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)[,opts$TFs]

sum(virtual_chip.dt[tf=="GATA1",correlation_score]>=0.30,na.rm=T)

###################################################
## Plot score thresholds vs number binding sites ##
###################################################

seq.ranges <- seq(0.25,0.8,by=0.01)

to.plot <- seq.ranges %>% map(function(j) {
  data.table(
    tf = opts$TFs,
    min_score = j,
    log2_N = log2(colSums(virtual_chip.mtx>=j)+1)
  )
}) %>% rbindlist

p <- ggplot(to.plot[log2_N>0], aes_string(x="min_score", y="log2_N", color="tf")) +
  geom_line(size=1) +
  labs(y="Number of predicted binding sites (log2)", x="Minimum in silico binding score") +
  scale_color_brewer(palette="Dark2") +
  geom_vline(xintercept=0.3, linetype="dashed") +
  # coord_cartesian(xlim=c(0,0.9)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    # legend.position = "right",
    legend.position = c(.9,.65),
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"lineplot_number_binding_sites_per_tf.pdf"), width=4.5, height=5)
print(p)
dev.off()

##############################################################
## Compare the number of putative binding sites for each TF ##
##############################################################

opts$min_acc_threshold <- 0.15
opts$min_cor_threshold <- 0.15

# Model 1: DNA
# THIS IS NOT CORRECT, WE NEED TO RUN FOR ALL GENOME!
# model1.dt <- data.table(
#   tf = colnames(motifmatcher.se),
#   V1 = colSums(assay(motifmatcher.se,"motifMatches")),
#   model = "DNA"
# )

# Model 2: DNA + ATAC
model2.dt <- virtual_chip.dt %>%
  .[,sum(max_accessibility_score>=opts$min_acc_threshold & !is.na(motif_score),na.rm=T),by="tf"] %>% 
  .[,model:="DNA + ATAC"]

# Model 3: DNA + ATAC + RNA
model3.dt <- virtual_chip.dt %>%
  .[,sum(correlation_score>=opts$min_cor_threshold & max_accessibility_score>=opts$min_acc_threshold & !is.na(motif_score),na.rm=T), by="tf"] %>%
  .[,model:="DNA + ATAC + RNA"]

# max.score <- 0.75
seq.ranges <- seq(0,1,by=0.1); names(seq.ranges) <- as.character(1:length(seq.ranges))

to.plot <- rbindlist(list(model2.dt,model3.dt)) %>%
  .[,model:=factor(model,levels=c("DNA + ATAC + RNA","DNA + ATAC","DNA"))]

# p.barplot <- ggbarplot(to.plot, x="V1", y="tf", fill="model", color="black", position = position_dodge(width=0.75)) +
#   labs(y="Number of predicted binding sites", x="") +
#   theme(
#     axis.text = element_text(size=rel(0.75)),
#     legend.title = element_blank()
#   )
# 
# pdf(sprintf("%s/barplot_compare_models.pdf",io$outdir), width=8, height=6.5)
# print(p.barplot)
# dev.off()

##  Box plot ##

to.plot <- to.plot %>%
  # .[V1<=3e4] %>%
  .[,log10_N:=log10(V1)]

# p <- ggboxplot(to.plot, x="tf", fill="model", y="log10_N", outlier.shape=NA) +
#   labs(x="", y="Number of predicted binding sites per TF (log2)") +
#   scale_fill_brewer(palette="Dark2") +
#   # coord_cartesian(ylim=c(0,35000)) +
#   stat_compare_means(aes(label = paste0("p = ", ..p.format..)), method="t.test") +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "none",
#     axis.text = element_text(size=rel(0.75), color="black")
#     # axis.text.y = element_text(size=rel(0.70), color="black")
#   )
# 
# pdf(file.path(io$outdir,"boxplot_number_predicted_sites_per_tf.pdf"), width=5, height=4)
# print(p)
# dev.off()

## Barplot ##


to.plot[,tf:=factor(tf,levels=model2.dt %>% setorder(-V1) %>% .$tf)]

p <- ggbarplot(to.plot, x="tf", fill="model", y="V1", position = position_dodge(width = 0.8)) +
  labs(x="", y="Number of predicted binding sites") +
  scale_fill_brewer(palette="Dark2") +
  # coord_cartesian(ylim=c(0,35000)) +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    axis.text = element_text(size=rel(0.75), color="black")
    # axis.text.y = element_text(size=rel(0.70), color="black")
  )

pdf(file.path(io$outdir,"barplot_number_predicted_sites_per_tf.pdf"), width=5, height=4)
print(p)
dev.off()

##  Density plot ##

# to.plot <- to.plot[V1<=1e5 & V1>=100] %>% .[,log10_N:=log10(V1)]
# 
# p <- ggdensity(to.plot, x="log10_N", fill="model", y="..density..") +
#   # geom_vline(aes(xintercept=V1), color="black", data=to.plot[,median(value),by="variable"], linetype="dashed") +
#   labs(x="Number of predicted binding sites (log10)", y="Density") +
#   # scale_x_continuous(limits = c(0,165)) +
#   # scale_fill_brewer(palette="Dark2") +
#   # guides(color=F) +
#   theme(
#     axis.title = element_text(size=rel(0.75)),
#     legend.title = element_blank(),
#     axis.text.x = element_text(size=rel(0.70), color="black"),
#     axis.text.y = element_text(size=rel(0.70), color="black")
#   )
# # 
# # 
# pdf(file.path(io$outdir,"density_number_predicted_sites_per_tf.pdf"), width=5, height=4)
# print(p)
# dev.off()

################################################################
## Compare the number of putative binding sites for each peak ##
################################################################

# Model 1: DNA
# model1.dt <- data.table(
#   peak = rownames(motifmatcher.se),
#   V1 = rowSums(assay(motifmatcher.se,"motifMatches")),
#   model = "DNA"
# )

# Model 2: DNA + ATAC
model2.dt <- virtual_chip.dt %>%
  .[,sum(max_accessibility_score>=opts$min_acc_threshold & !is.na(motif_score),na.rm=T),by="peak"] %>% 
  .[,model:="DNA + ATAC"]

# Model 3: DNA + ATAC + RNA
model3.dt <- virtual_chip.dt %>%
  .[,sum(correlation_score>=opts$min_cor_threshold & max_accessibility_score>=opts$min_acc_threshold & !is.na(motif_score),na.rm=T), by="peak"] %>%
  .[,model:="DNA + ATAC + RNA"]

# max.score <- 0.75
seq.ranges <- seq(0,1,by=0.1); names(seq.ranges) <- as.character(1:length(seq.ranges))

to.plot <- rbindlist(list(model2.dt,model3.dt)) %>%
  .[,model:=factor(model,levels=c("DNA + ATAC + RNA","DNA + ATAC","DNA"))]

##  Box plot ##

p <- ggboxplot(to.plot, x="model", fill="model", y="V1", outlier.shape=NA) +
  labs(x="", y="Number of predicted binding sites per peak") +
  scale_fill_brewer(palette="Dark2") +
  coord_cartesian(ylim=c(0,180)) +
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), method="t.test") +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    # axis.text.x = element_text(size=rel(0.70), color="black"),
    # axis.text.y = element_text(size=rel(0.70), color="black")
  )

pdf(file.path(io$outdir,"boxplot_number_predicted_sites_per_peak.pdf"), width=4, height=5)
print(p)
dev.off()

## Density plot ##
# to.plot <- to.plot %>%
#   .[V1>1] %>%
#   .[,log10_N:=log10(V1+1)]

# p <- ggdensity(to.plot, x="log10_N", fill="model", y="..density..") +
#   # geom_vline(aes(xintercept=V1), color="black", data=to.plot[,median(value),by="variable"], linetype="dashed") +
#   labs(x="Number of predicted binding sites (log10)", y="Density") +
#   # scale_x_continuous(limits = c(0,165)) +
#   # scale_fill_brewer(palette="Dark2") +
#   # guides(color=F) +
#   theme(
#     axis.title = element_text(size=rel(0.85)),
#     legend.title = element_blank(),
#     axis.text.x = element_text(size=rel(0.70), color="black"),
#     axis.text.y = element_text(size=rel(0.70), color="black")
#   )
# 
# pdf(file.path(io$outdir,"density_number_predicted_sites_per_peak.pdf"), width=5, height=4)
# print(p)
# dev.off()

##################################################################
## Plot score thresholds vs number of (+) and (-) binding sites ##
##################################################################

# sum(virtual_chip.mtx[,"SOX10"]>=0.15)
# virtual_chip.dt[tf=="SOX10" & score>=0.15]
## Lineplot for TF binding sites ##

seq.ranges <- seq(0.15,0.70,by=0.05)

to.plot <- seq.ranges %>% map(function(j) {
  data.table(
    tf = opts$TFs,
    min_score = j,
    N_negative = colSums(virtual_chip.mtx<=(-j)),
    N_positive = colSums(virtual_chip.mtx>=j)
  )
}) %>% rbindlist %>% melt(id.vars=c("tf","min_score"), variable.name="cor_sign", value.name="N")

to.plot %>% .[,log10_N:=log10(N+1)]
# active.tfs <- celltype_virtual_chip_to_plot.dt[min_score>=0.40,max(log2_N),by="tf"] %>% .[V1>=6] %>% .$tf

p <- ggplot(to.plot[N>0], aes_string(x="min_score", y="log10_N", color="cor_sign")) +
  geom_line(size=1) +
  facet_wrap(~tf, scales="free_y") +
  labs(y="Number of predicted binding sites (log10)", x="In silico TF binding score") +
  scale_color_brewer(palette="Set1") +
  # coord_cartesian(xlim=c(0,0.9)) +
  theme_classic() +
  theme(
    axis.text = element_text(color="black"),
    # legend.position = "right",
    legend.position = "top",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"lineplot_number_binding_sites_per_tf_positive_vs_negative.pdf"), width=7, height=5)
print(p)
dev.off()

