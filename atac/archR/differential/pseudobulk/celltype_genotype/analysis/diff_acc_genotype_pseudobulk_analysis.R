# here::i_am("atac/archR/differential/wt_vs_ko/analysis/analysis.R")

source(here::here("settings.R"))
source(here::here("utils.R"))
# source(here::here("atac/archR/differential/pseudobulk/utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$min.logFC <- 1.5        # Minimum logFC
opts$min.acc <- 2          # Minimum accessibility per group
opts$threshold_fdr <- 0.001
opts$atac.matrix <- "PeakMatrix"

opts$wt.class <- "WT"
opts$ko.class <- "KO"

# I/O
io$atac_marker_peaks <- file.path(io$basedir,"results/atac/archR/differential/pseudobulk/celltype/PeakMatrix/markers/PeakMatrix_markers_upregulated_filt.txt.gz")
io$diff <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/celltype_genotype/%s/parsed/diff_results.txt.gz",opts$atac.matrix))
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/pseudobulk/celltype_genotype/%s/pdf",opts$atac.matrix))

dir.create(io$outdir, showWarnings = F)
dir.create(file.path(io$outdir,"volcano_plots"), showWarnings = F)
dir.create(file.path(io$outdir,"polar_plots"), showWarnings = F)

##################
## Load results ##
##################

diff.dt <- fread(io$diff) %>%
  .[, sign := ifelse(logFC>0,"Upregulated in KO","Downregulated in KO")] %>%
  .[, sig := (padj_fdr<=opts$threshold_fdr & abs(logFC)>=opts$min.logFC & (mean_groupA>=opts$min.acc | mean_groupB>=opts$min.acc))]

##################
## Filter results ##
##################

# Only consider cell types with enough number of cells
stats.dt <- fread("/Users/argelagr/data/gastrulation_multiome_10x/test/results/atac/archR/pseudobulk/celltype_genotype/PeakMatrix/stats.txt") %>%
  .[,celltype:=strsplit(group,"-") %>% map_chr(1)] %>% .[,genotype:=strsplit(group,"-") %>% map_chr(2)] %>%
  dcast(celltype~genotype, value.var="N", fill=0)
celltypes.to.use <- stats.dt[T_KO>=50 & WT>=50,celltype]

celltypes.to.use <- celltypes.to.use[!celltypes.to.use%in%c("Parietal_endoderm","ExE_endoderm","Intermediate_mesoderm")]

diff.dt <- diff.dt[celltype%in%celltypes.to.use]

# Subset to marker peaks
marker_peaks.dt <- fread(io$atac.markers_peaks.pseudobulk.filt)
diff_markers.dt <- diff.dt[feature%in%unique(marker_peaks.dt$feature)]

################
## Polar plot ##
################

to.plot <- diff_markers.dt %>% copy %>%
  .[,.(N=sum(sig,na.rm=T)) ,by=c("celltype")]

# to.plot[N>=100,N:=100]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  coord_polar() +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.text = element_text(size=rel(0.75)),
    legend.title = element_blank(),
    axis.title=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.line=element_blank(),
    axis.text.x = element_blank()
    # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
  )

pdf(file.path(io$outdir,"DA_polar_plots_wt_vs_ko.pdf"), width=11, height=8)
print(p)
dev.off()

##############
## Bar plot ##
##############

to.plot <- diff_markers.dt[sig==T] %>% .[,.N, by=c("celltype")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = celltype), color="black", stat = 'identity') + 
  scale_fill_manual(values=opts$celltype.colors, drop=F) +
  labs(x="", y="Number of DA peaks") +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    # axis.text.x = element_text(color="black", size=rel(0.75))
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


pdf(file.path(io$outdir,"DA_barplots_wt_vs_ko.pdf"), width=12, height=4)
print(p)
dev.off()

#############################
## Bar plot, split by sign ##
#############################

to.plot <- diff_markers.dt[sig==T] %>% .[,.N, by=c("celltype","sign")]

p <- ggplot(to.plot, aes(x=celltype, y=N)) +
  geom_bar(aes(fill = sign), color="black", stat="identity") + 
  labs(x="", y="Number of DA peaks") +
  guides(x = guide_axis(angle = 90)) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.75))
  )

pdf(file.path(io$outdir,"DA_barplots_sign_wt_vs_ko.pdf"), width=11, height=8)
print(p)
dev.off()


###################
## Volcano plots ##
###################

label_groups <- c("WT", "T_KO")

# i <- "NMP"
for (i in unique(diff_markers.dt$celltype)) {
  # to.plot <- diff.dt[celltype==i] %>% .[!is.na(sig)] %>% .[abs(logFC)>=0.5]
  to.plot <- diff_markers.dt[celltype==i] %>% .[!is.na(sig)] %>% .[abs(logFC)>=0.25]
  to.plot %>% .[,log_padj_fdr:=-log10(padj_fdr)] %>% .[log_padj_fdr>=40,log_padj_fdr:=40]
  # p <- gg_volcano_plot(to.plot, top_genes=0, label_groups = c("WT","KO"), xlim = 75, ylim = 20)
  
  negative_hits <- to.plot[sig==TRUE & logFC<0,feature]
  positive_hits <- to.plot[sig==TRUE & logFC>0,feature]
  all <- nrow(to.plot)
  
  # filter small peaks for better viz
  to.plot <- to.plot %>% .[mean_groupA>=1 | mean_groupB>=1]
  to.plot <- rbind(to.plot[sig==T], to.plot[sig==F][sample(.N/2.5)])
  xlim <- max(abs(to.plot$logFC), na.rm=T)
  ylim <- max(-log10(to.plot$padj_fdr), na.rm=T)
  
  p <- ggplot(to.plot, aes(x=logFC, y=log_padj_fdr)) +
    geom_point(aes(color=sig, size=log_padj_fdr)) +
    # ggrastr::geom_point_rast(aes(color=sig), size=0.75) +
    geom_segment(x=0, xend=0, y=0, yend=ylim-3, color="orange") +
    scale_size_continuous(range = c(0.01,1.25)) + 
    scale_color_manual(values=c("black","red")) +
    scale_x_continuous(limits=c(-xlim-0.5,xlim+0.5)) +
    # scale_x_continuous(limits=c(-xlim,xlim)) +
    # scale_y_continuous(limits=c(0,ylim+2.5)) +
    labs(x="Accessibility difference (%)", y=expression(paste("-log"[10],"(FDR)"))) +
    annotate("text", x=0, y=ylim+0.5, size=5, label=sprintf("(%d)", all)) +
    annotate("text", x=-7, y=ylim+0.5, size=5, label=sprintf("%d (-)",length(negative_hits))) +
    annotate("text", x=7, y=ylim+0.5, size=5, label=sprintf("%d (+)",length(positive_hits))) +
    # ggrepel::geom_text_repel(data=head(to.plot[sig==T],n=top_genes), aes(x=logFC, y=-log10(FDR), label=symbol), size=5) +
    guides(size="none") +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(1), color='black'),
      axis.title = element_text(size=rel(1), color='black'),
      # axis.title = element_text(),
      legend.position="none"
    )
  
  if (length(label_groups)>0) {
    p <- p +
      annotate("text", x=-6, y=0, size=4.5, label=sprintf("Up in %s",label_groups[2])) +
      annotate("text", x=7, y=0, size=4.5, label=sprintf("Up in %s",label_groups[1]))
  }
  

  pdf(file.path(io$outdir,sprintf("volcano_plots/%s_WT_vs_KO_volcano.pdf",i)), width=7, height=5)
  # png(sprintf("%s/%s_WT_vs_KO_DA_peaks_volcano.png",io$outdir,i,j), width=900, height=500)
  print(p)
  dev.off()
}

################################
## Cell fate bias polar plots ##
################################

to.plot <- diff_markers.dt %>%
  merge(
    marker_peaks.dt[,c("feature","celltype")] %>% setnames("celltype","celltype_marker"), by = "feature", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","sign")]

to.plot <- to.plot[!celltype_marker%in%c("PGC")]
to.plot[V1>=100,V1:=100] # for viz

# i <- "NMP"
for (i in unique(to.plot$celltype)) {
  
  p <- ggplot(to.plot[celltype==i], aes(x=celltype_marker, y=V1)) +
    geom_bar(aes(fill = celltype_marker), color="black", stat = 'identity') + 
    facet_wrap(~sign, nrow=1) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    scale_y_continuous(limits = c(0,100)) +
    coord_polar() +
    theme_bw() +
    theme(
      legend.position = "none",
      # legend.text = element_text(size=rel(0.75)),
      # legend.title = element_blank(),
      axis.title=element_blank(),
      # axis.text.y=element_blank(),
      axis.text.y=element_text(),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
      # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
    )
  
  pdf(file.path(io$outdir,sprintf("polar_plots/%s_DA_polar_plot_fate_bias.pdf",i)), width=6, height=3)
  print(p)
  dev.off()
}

###################################
## Calculate TF motif enrichment ##
###################################

# Load Motif annotations 
peak_annotation.se <- readRDS(io$motifmatcher_scores.se)[unique(diff_markers.dt$feature),]

motifs.to.use <- colnames(peak_annotation.se)# %>% head(n=10)
# motifs.to.use <- c("TFAP2D_2","MESP1_69","TBXT")
celltypes.to.use <- unique(diff_markers.dt$celltype)
# celltypes.to.use <- c("NMP","Spinal_cord")

motif_enrichment.dt <- motifs.to.use %>% map(function(i) { 
  print(i)
  
  celltypes.to.use %>% map(function(j) {
    
    foreground.peaks <- diff_markers.dt[celltype==j & sign=="Downregulated in KO" & sig==TRUE,feature]
    background.peaks <- diff_markers.dt[celltype==j & sig==FALSE,feature]
    foreground.nmatches <- assay(peak_annotation.se,"motifMatches")[,i][foreground.peaks] %>% sum
    background.nmatches <- assay(peak_annotation.se,"motifMatches")[,i][background.peaks] %>% sum
    
    p.value <- phyper(foreground.nmatches-1, background.nmatches, length(background.peaks)-background.nmatches,length(foreground.peaks), lower.tail = F)
    # fisher.test(matrix(c(Overlap, group2-Overlap, group1-Overlap, Total-group2-group1 +Overlap), 2, 2), alternative='greater')$p.value
    # background.total = length(background.peaks)
    # p.value <- phyper(
    #   q = foreground.nmatches-1, 
    #   m = length(foreground.peaks), 
    #   n = length(background.peaks), 
    #   k = foreground.nmatches+background.nmatches,
    #   lower.tail = F
    # )
    # data.table(motif=i, celltype=j, pval=format(p.value,digits=3))
    data.table(motif=i, celltype=j, pval=p.value)
}) %>% rbindlist }) %>% rbindlist

fwrite(motif_enrichment.dt, file.path(io$outdir,"motif_enrichment_DA_peaks_genotype.txt.gz"), sep="\t", na="NA", quote=F)

# Load precomputed
motif_enrichment.dt <- fread(file.path(io$outdir,"motif_enrichment_DA_peaks_genotype.txt.gz")) %>%
  .[celltype%in%celltypes.to.use]

##############################
## Plot TF motif enrichment ##
##############################

to.plot <- motif_enrichment.dt[pval<=0.10] %>% 
  .[,pval:=as.numeric(pval)] %>%
  .[,log_pval:=-log10(pval)] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes[opts$celltypes%in%celltypes.to.use])]

to.plot[,dot_size:=minmax.normalisation(abs(log_pval))]
  
# cap p-values for better viz
to.plot <- to.plot[log_pval>=40,log_pval:=40]
to.plot.text <- to.plot[log_pval>=20]
  
p <- ggplot(to.plot, aes_string(x="celltype", y="log_pval", size="dot_size", fill="celltype")) +
  geom_point(shape=21) +
  ggrepel::geom_text_repel(data=to.plot.text, aes(x=celltype, y=log_pval, label=motif), size=3,  max.overlaps=100, segment.color = NA) +
  # scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  scale_fill_manual(values=opts$celltype.colors) +
  scale_x_discrete(drop=F) +
  scale_size_continuous(range = c(0.01,7)) +
  guides(x = guide_axis(angle = 90), size="none") +
  labs(x="", y=expression(paste("-log"[10],"(p.value)"))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black")
  )

pdf(file.path(io$outdir,"DA_peaks_motif_enrichment_per_celltype.pdf"), width=3.5, height=5)
print(p)
dev.off()

####################################
## Plot TF motif enrichment vs DE ##
####################################

motif2gene.dt <- fread(io$archR.motif2gene) %>% .[,c("motif","gene")]
io$diff.expr <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype_genotype/parsed/diff_expr_results.txt.gz")

diff_expr.dt <- fread(io$diff.expr) %>% 
  .[,gene:=toupper(gene)] %>%
  .[gene%in%motif2gene.dt$gene] %>%
  .[celltype%in%celltypes.to.use] %>%
  .[is.na(logFC),logFC:=0]

to.plot <- motif_enrichment.dt %>% 
  .[,log_pval:=-log10(as.numeric(pval))] %>%
  .[,celltype:=factor(celltype,levels=celltypes.to.use)] %>%
  merge(motif2gene.dt,by="motif") %>%
  merge(diff_expr.dt[,c("logFC","gene","celltype")],by=c("gene","celltype"))

ggplot(to.plot, aes_string(x="log_pval", y="logFC")) +
  geom_point(shape=21) +
  facet_wrap(~celltype) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black")
  )


#####################################
## Overlay with in silico chip-seq ##
#####################################

# Load virtual ChIP-seq
io$virtual_chip <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/JASPAR/virtual_chip_matrix.rds")
virtual_chip.mtx <- readRDS(io$virtual_chip)
virtual_chip.mtx[virtual_chip.mtx<=0] <- 0

# convert to binary matrix
opts$min.chip.score <- 0.25
virtual_chip_logical.mtx <- virtual_chip.mtx
virtual_chip_logical.mtx[virtual_chip_logical.mtx>=opts$min.chip.score] <- 1
virtual_chip_logical.mtx[virtual_chip_logical.mtx<opts$min.chip.score] <- 0

# Select DA peaks
diff_acc_peaks <- diff_markers.dt[celltype=="NMP" & sign=="Downregulated in KO" & sig==TRUE,feature]

####################################################
## Brachyury in silico ChIP-seq score in DA peaks ##
####################################################

motifmatcher.se <- readRDS(io$motifmatcher_scores.se)

peaks.with.motif <- which(assay(motifmatcher.se[,"TBXT"],"motifMatches")[,1]) %>% names

peaks <- Reduce("intersect",list(diff_markers.dt[celltype=="NMP",feature], rownames(virtual_chip.mtx),peaks.with.motif))
da_results.dt <- diff_markers.dt[feature%in%peaks &celltype=="NMP"] 
da_results.dt$insilico_chip_score <- virtual_chip.mtx[da_results.dt$feature,"T"]

to.plot <- da_results.dt[insilico_chip_score>=0 & sign=="Downregulated in KO"]

p <- ggboxplot(to.plot, x="sig", y="insilico_chip_score", fill="sig", outlier.shape=NA) +
  labs(x="Differential acc. of peaks with T motif", y="Virtual ChIP-seq score") +
  # geom_boxplot(outlier.shape=NA) + 
  stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.y=0.40, method="t.test") +
  coord_cartesian(ylim=c(0,0.45)) +
  scale_fill_brewer(palette="Dark2") +
  # geom_violin() +
  theme_classic() +
  theme(
    axis.text.y = element_text(color="black"),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top"
  )

pdf(file.path(io$outdir,"DA_peaks_virtual_chip_score.pdf"), width=3.5, height=5)
print(p)
dev.off()



##########
## TEST ##
##########

# opts$min.logFC <- 1.5        # Minimum logFC
# opts$min.acc <- 2          # Minimum accessibility per group
# 
# tmp <- sort(virtual_chip.mtx[,"T"]) %>% tail(n=100)
# tmp <- tmp[names(tmp)%in%unique(marker_peaks.dt$feature)]
# 
# names(tmp)%in%diff_markers.dt[celltype=="NMP" & sig==T,feature]
# 
# diff_markers.dt[feature=="chr19:37874126-37874726" & celltype=="NMP"]
