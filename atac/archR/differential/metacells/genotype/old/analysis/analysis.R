here::i_am("atac/archR/differential/wt_vs_ko/analysis/analysis.R")

source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("atac/archR/differential/utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$min.MeanDiff <- 0.10
opts$threshold_fdr <- 0.01
opts$atac.matrix <- "PeakMatrix"

opts$wt.class <- "WT"
opts$ko.class <- "KO"

# I/O
io$atac_marker_peaks <- file.path(io$basedir,"results/atac/archR/differential/celltype.mapped/PeakMatrix/markers/marker_peaks_upregulated_filtered.txt.gz")
io$indir <- file.path(io$basedir,sprintf("results/atac/archR/differential/wt_vs_ko/%s",opts$atac.matrix))
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/differential/wt_vs_ko/%s/pdf",opts$atac.matrix))

dir.create(io$outdir, showWarnings = F)
dir.create(file.path(io$outdir,"volcano_plots"), showWarnings = F)
dir.create(file.path(io$outdir,"polar_plots"), showWarnings = F)

##################
## Load results ##
##################

source(here::here("atac/archR/differential/wt_vs_ko/analysis/load_data.R"))

####################
## Filter results ##
####################

# Filter by number of cells
# diff.dt <- diff.dt[groupA_N>=opts$min.cells & groupB_N>=opts$min.cells]

# Subset to marker peaks
marker_peaks.dt <- fread(io$atac_marker_peaks)
diff_markers.dt <- diff.dt[idx%in%unique(marker_peaks.dt$idx)]

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

for (i in unique(diff_markers.dt$celltype)) {
  to.plot <- diff.dt[celltype==i] %>% .[!is.na(sig)] %>% .[,MeanDiff:=100*MeanDiff] %>% .[abs(MeanDiff)>=2]

  p <- gg_volcano_plot(to.plot, top_genes=0, label_groups = c("WT","KO"), xlim = 75, ylim = 20)

  pdf(file.path(io$outdir,sprintf("volcano_plots/%s_WT_vs_KO_volcano.pdf",i)), width=9, height=5)
  # png(sprintf("%s/%s_WT_vs_KO_DA_peaks_volcano.png",io$outdir,i,j), width=900, height=500)
  print(p)
  dev.off()
}

################################
## Cell fate bias polar plots ##
################################

to.plot <- diff_markers.dt %>%
  merge(
    marker_peaks.dt[,c("idx","celltype")] %>% setnames("celltype","celltype_marker"), by = "idx", allow.cartesian=TRUE
  ) %>% .[,sum(sig), by=c("celltype","celltype_marker","sign")]

to.plot[V1>=50,V1:=50] # for viz

# i <- "NMP"
for (i in unique(to.plot$celltype)) {
  
  p <- ggplot(to.plot[celltype==i], aes(x=celltype_marker, y=V1)) +
    geom_bar(aes(fill = celltype_marker), color="black", stat = 'identity') + 
    facet_wrap(~sign, nrow=1) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    scale_y_continuous(limits = c(0,50)) +
    coord_polar() +
    theme_bw() +
    theme(
      legend.position = "none",
      # legend.text = element_text(size=rel(0.75)),
      # legend.title = element_blank(),
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
      # axis.text.x = element_text(angle= -76 - 360 / length(unique(to.plot$celltype)) * seq_along(to.plot$celltype))
    )
  
  pdf(file.path(io$outdir,sprintf("polar_plots/%s_DA_polar_plot_fate_bias.pdf",i)), width=6, height=3)
  print(p)
  dev.off()
}

#########################
## TF motif enrichment ##
#########################

# Load Motif annotations 
peak_annotation.se <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/processed/atac/archR/Annotations/CISBP-Scores.rds")
peak_annotation_filt.se <- peak_annotation.se[unique(diff_markers.dt$idx),]

motifs.to.use <- colnames(peak_annotation_filt.se)# %>% head(n=10)
celltypes.to.use <- unique(diff_markers.dt$celltype)

# library(furrr)
# opts$ncores <- 4
# if (opts$ncores>1){
#   plan(multicore, workers=opts$ncores)
#   motifs_split <- split(motifs.to.use, cut(seq_along(motifs.to.use), opts$ncores, labels = FALSE)) 
# } else {
#   plan(sequential)
#   motifs_split <- list(motifs.to.use)
# }

# i="T_789"; j <- "NMP"
motif_enrichment.dt <- motifs.to.use %>% map(function(i) { 
  print(i)
  
  celltypes.to.use %>% map(function(j) {
    
    foreground.peaks <- diff_markers.dt[celltype==j & sign=="Downregulated in KO" & sig==TRUE,idx]
    background.peaks <- diff_markers.dt[celltype==j & sig==FALSE,idx]
    foreground.nmatches <- assay(peak_annotation_filt.se,"motifMatches")[,i][foreground.peaks] %>% sum
    background.nmatches <- assay(peak_annotation_filt.se,"motifMatches")[,i][background.peaks] %>% sum
    
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

fwrite(motif_enrichment.dt, file.path(io$outdir,"motif_enrichment_DA_peaks_wt_vs_ko.txt.gz"), sep="\t", na="NA", quote=F)

to.plot <- motif_enrichment.dt[pval<=0.10] %>% 
  .[,pval:=as.numeric(pval)] %>%
  .[,log_pval:=-log10(pval)] %>%
  .[,celltype:=factor(celltype,levels=celltypes.to.use)]

to.plot[,dot_size:=minmax.normalisation(abs(log_pval))]
  
to.plot.text <- to.plot[pval<=0.01]
  
ggplot(to.plot, aes_string(x="celltype", y="log_pval", group="motif", size="dot_size")) +
  geom_point(shape=21) +
  scale_x_discrete(drop=F) +
  ggrepel::geom_text_repel(data=to.plot.text, aes(x=diff_CG, y=diff_GC, label=tf), size=3,  max.overlaps=100, segment.color = NA) +s
  # scale_fill_gradient2(low = "gray50", mid="gray90", high = "red") +
  scale_size_continuous(range = c(0.25,2)) +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black")
  )

#####################################
## Overlay with in silico chip-seq ##
#####################################
  
  
# io$virtual_chip_mtx <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/pseudobulk/CISBP/virtual_chip_matrix.rds")
io$virtual_chip_mtx <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/CISBP/virtual_chip_matrix.rds")
virtual_chip.mtx <- readRDS(io$virtual_chip_mtx)

virtual_chip_logical.mtx <- virtual_chip.mtx
virtual_chip_logical.mtx[virtual_chip_logical.mtx>=0.25] <- 1
virtual_chip_logical.mtx[virtual_chip_logical.mtx<0.25] <- 0

diff_acc_peaks <- diff_markers.dt[celltype=="NMP" & sign=="Downregulated in KO" & sig==TRUE,idx]
# to.plot <- data.table(
#   tf = colnames(virtual_chip.mtx),
#   N = colSums(virtual_chip.mtx[diff_acc_peaks,]>=0.25)
# )

# i="T_789"; j <- "NMP"
tfs.to.use <- colnames(virtual_chip.mtx)
celltypes.to.use <- "NMP"
# i <- "T"
tf_enrichment_insilico_chip.dt <- tfs.to.use %>% map(function(i) { 
  print(i)
  
  celltypes.to.use %>% map(function(j) {
    
    foreground.peaks <- diff.dt[celltype==j & sign=="Downregulated in KO" & sig==TRUE,idx]
    background.peaks <- diff.dt[celltype==j & sig==FALSE,idx]
    foreground.nmatches <- virtual_chip_logical.mtx[foreground.peaks,i] %>% sum
    background.nmatches <- virtual_chip_logical.mtx[background.peaks,i] %>% sum
    
    p.value <- phyper(foreground.nmatches-1, background.nmatches, length(background.peaks)-background.nmatches,length(foreground.peaks), lower.tail = F)

    data.table(tf=i, celltype=j, pval=p.value)
  }) %>% rbindlist }) %>% rbindlist


to.plot <- tf_enrichment_insilico_chip.dt[pval<=0.10] %>% 
  .[,pval:=as.numeric(pval)] %>%
  .[,log_pval:=-log10(pval)] %>%
  .[,celltype:=factor(celltype,levels=celltypes.to.use)]

to.plot[,dot_size:=minmax.normalisation(abs(log_pval))]

to.plot.text <- to.plot[pval<=0.01]

ggplot(to.plot[pval<=0.01], aes_string(x="log_pval", y="tf", size="dot_size")) +
  geom_point(shape=21) +
  # scale_x_discrete(drop=F) +
  # scale_size_continuous(range = c(0.25,2)) +
  # guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black")
  )

##############################
## Network of TF2peak links ##
##############################

# TO-DO: SUBSET PEAKS WITH THE T MOTIF
motifmatcher.se <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/processed/atac/archR/Annotations/JASPAR-Scores.rds")

peaks.with.motif <- which(assay(motifmatcher.se[,"TBXT_121"],"motifMatches")[,1]) %>% names
peaks <- Reduce("intersect",list(diff_markers.dt[celltype=="NMP",idx], rownames(virtual_chip.mtx),peaks.with.motif))
da_results.dt <- diff_markers.dt[idx%in%peaks &celltype=="NMP"] 
da_results.dt$insilico_chip_score <- virtual_chip.mtx[da_results.dt$idx,"T"]

to.plot <- da_results.dt[insilico_chip_score>=0 & sign=="Downregulated in KO"]

ggscatter(to.plot, x="MeanDiff", y="insilico_chip_score") +
  stat_cor(method = "pearson") +
  stat_smooth(method="lm")

# to.plot <- 
# ggscatter(edge_list.dt, x="weight", y="da")
tmp <- virtual_chip.mtx[da_results.dt$idx,"T"]  

edge_list.dt <- data.table(
  from = "T",
  to = names(tmp),
  weight = tmp,
  da = da_results.dt$MeanDiff
)


# Create igraph object
igraph.net <- graph_from_data_frame(d = edge_list.dt)

igraph::degree(igraph.net)
