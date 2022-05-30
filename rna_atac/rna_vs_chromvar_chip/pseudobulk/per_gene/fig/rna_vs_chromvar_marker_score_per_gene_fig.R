source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

opts$motif_annotation <- "CISBP"

io$tf_markers_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_tfs_all.txt.gz")
io$tf_markers_atac <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/CISBP/parsed/markers_all.txt.gz")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_gene/%s/pdf/marker_scores",opts$motif_annotation))

dir.create(io$outdir, showWarnings = F)

################################################
## Load pseudobulk RNA and chromVAR estimates ##
################################################

source(here::here("rna_atac/load_rna_atac_pseudobulk.R"))

############################
## Load TF marker scores  ##
############################

tf_markers_rna.dt <- fread(io$tf_markers_rna) %>% setnames("score","tf_marker_score_rna")
tf_markers_atac.dt <- fread(io$tf_markers_atac) %>% setnames("score","tf_marker_score_atac")

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

##########
## Plot ##
##########

# genes.to.plot <- c("FOXA2","RFX4","TAL1")
genes.to.plot <- c("ETS1","TFAP2C","SOX9")
# genes.to.plot <- unique(rna_chromvar.dt$gene)

for (i in genes.to.plot) {
  
  to.plot <- rna_chromvar.dt[gene==i]  %>%
      merge(tf_markers_rna.dt[gene==i,c("celltype","gene","tf_marker_score_rna")],by=c("gene","celltype")) %>%
      merge(tf_markers_atac.dt[gene==i,c("celltype","gene","tf_marker_score_atac")],by=c("gene","celltype"))
  
  to.plot %>% 
    .[,tf_marker_score:=tf_marker_score_rna*tf_marker_score_atac] %>% 
    .[,tf_marker_score:=tf_marker_score/max(tf_marker_score)]
  

  # to.plot %>% .[,dot_size:=minmax.normalisation(tf_marker_score)]
  # to.plot.text <- to.plot[tf_marker_score>=0.75]
  # pos.low <- position_jitter(width = 1, height=0.12, seed = 42)
  # pos.high <- position_jitter(width = 0.5, height=0, seed = 42)
  
  p <- ggplot(to.plot, aes(x=celltype, y=tf_marker_score)) +
    # geom_jitter(aes(size=dot_size, fill=celltype), stat = 'identity', shape=21, color="black", alpha=0.9, position=pos.low, data=to.plot[tf_marker_score<0.5]) +
    # geom_jitter(aes(size=dot_size, fill=celltype), stat = 'identity', shape=21, color="black", alpha=0.9, data=to.plot[tf_marker_score>0.5], position=pos.high) +
    geom_bar(aes(fill=celltype), stat = 'identity', color="black") +
    ggrepel::geom_text_repel(aes(label=gsub("_"," ",celltype)), size=4, max.overlaps=Inf, data=to.plot[tf_marker_score>0.5]) +
    scale_fill_manual(values=opts$celltype.colors) +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by=0.25)) +
    coord_polar() + 
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(size=rel(1.5)),
      axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line.x=element_blank()
    )
  
  pdf(file.path(io$outdir,sprintf("%s_tf_marker_score_polar_plot.pdf",i)), width=5, height=4)
  print(p)
  dev.off()
  
}


#######################################
## Stacked barplots per TF (all TFs) ##
#######################################

genes.to.plot <- unique(tf_markers.dt[tf_marker_score>=0.75,gene])

to.plot <- tf_markers.dt %>%
  .[gene%in%genes.to.plot] %>%
  # .[,var:=var(chromvar_zscore),by="gene"] %>% .[var>0] %>% .[,var:=NULL] %>%
  .[,value:=minmax.normalisation(tf_marker_score), by="gene"]

# Sort by euclidean distance
tmp <- to.plot %>% dcast(celltype~gene,value.var="tf_marker_score",fill=0) %>% matrix.please
TF.order <- colnames(tmp)[order(irlba::prcomp_irlba(t(tmp), n=1)$x[,1])]
TF.order <- TF.order[TF.order%in%to.plot$gene]
to.plot[,gene:=factor(gene,levels=TF.order)]

# to.plot <- to.plot[value>0.25] 

p <- ggplot(to.plot, aes(x=gene, y=value)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black", position="fill", size=0.15) +
  # facet_wrap(~cor_sign, nrow=2, scales="free_x") +
  scale_fill_manual(values=opts$celltype.colors[names(opts$celltype.colors)%in%unique(to.plot$celltype)]) +
  theme_classic() +
  labs(x="", y="TF marker score (scaled)") +
  # guides(x = guide_axis(angle = 90)) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    # axis.text.x = element_text(color="black", size=rel(0.85)),
    axis.text.y = element_text(color="black", size=rel(0.8)),
    # axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(file.path(io$outdir,"TF_marker_score_stacked_barplots.pdf"), width=9, height=3.5)
print(p)
dev.off()

