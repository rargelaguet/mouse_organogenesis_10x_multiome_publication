here::i_am("rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype/rna_vs_chromvar_pseudobulk_per_celltype.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

opts$motif_annotation <- "CISBP"
io$rna_sce_pseudobulk <- file.path(io$basedir,"results/rna/pseudobulk/celltype/SingleCellExperiment_TFs_pseudobulk.rds")
io$atac_chromvar_chip_pseudobulk <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk/chromVAR_chip_%s_archr.rds",opts$motif_annotation))
io$tf_markers_rna <- file.path(io$basedir,"results/rna/differential/pseudobulk/celltype/parsed/marker_tfs_all.txt.gz")
io$tf_markers_atac <- file.path(io$basedir,"results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltypes/CISBP/parsed/markers_all.txt.gz")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype/%s/fig",opts$motif_annotation))

# I/O
dir.create(io$outdir, showWarnings=F, recursive = T)

#######################################
## Load pseudobulk RNA and ATAC data ##
#######################################

# Load pseudobulk TF RNA expression
rna_pseudobulk.sce <- readRDS(io$rna_sce_pseudobulk)

# Load pseudobulk ATAC chromVAR-ChIP matrix
atac_chromvar_pseudobulk.se <- readRDS(io$atac_chromvar_chip_pseudobulk)

# Fetch common TFs
TFs <- intersect(rownames(rna_pseudobulk.sce),rownames(atac_chromvar_pseudobulk.se))
rna_pseudobulk.sce <- rna_pseudobulk.sce[TFs,]
atac_chromvar_pseudobulk.se <- atac_chromvar_pseudobulk.se[TFs,]

########################
## Prepare data table ##
########################

atac_chromvar_pseudobulk.dt <- assay(atac_chromvar_pseudobulk.se,"z") %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","celltype") %>%
  melt(id.vars=c("celltype"), variable.name="gene", value.name="chromvar_zscore")

rna_tf_pseudobulk.dt <- logcounts(rna_pseudobulk.sce) %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  data.table::melt(id.vars="gene", variable.name="celltype", value.name="expr")

###########
## Merge ##
###########

rna_chromvar.dt <- merge(
  rna_tf_pseudobulk.dt,
  atac_chromvar_pseudobulk.dt,
  by = c("celltype","gene")
)

length(unique(rna_chromvar.dt$gene))

######################
## Load TF evidence ##
######################

tf_evidence.dt <- fread("/Users/argelagr/data/gastrulation_multiome_10x/test/results/rna_atac/rna_vs_chromvar_chip/pseudobulk/per_celltype/tf_evidence.txt", sep=",") %>%
  .[,c("gene","celltype","species")] %>% setnames("species","evidence") %>% .[is.na(evidence),evidence:="None"]

table(tf_evidence.dt$evidence)

tf_evidence.colors <- c(
  "zebrafish" = "#1B9E77", 
  "mouse" = "#D95F02", 
  "human" = "#7570B3", 
  # "in silico" = "#7570B3", 
  "chick" = "#E7298A", 
  "frog" = "#66A61E",
  "None" = "gray80"
)

# library(RColorBrewer)
# tmp <- brewer.pal(length(unique(tf_evidence.dt$evidence)), "Dark2")
# names(tmp) <- unique(tf_evidence.dt$evidence)
# dput(tmp)

###########################
## Load TF marker scores ##
###########################

tf_markers_rna.dt <- fread(io$tf_markers_rna) %>% setnames("score","tf_marker_score_rna")
tf_markers_atac.dt <- fread(io$tf_markers_atac) %>% setnames("score","tf_marker_score_atac")

########################################
## Sort TFs based on motif similarity ##
########################################

# Load motif similarity 
io$motif_similarity <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/motif_similarity/%s/motif_similarity.rds",opts$motif_annotation))
motif_similarity.mtx <- readRDS(io$motif_similarity)

# Multidimensional scaling
motif_similarity.mds <- stats::cmdscale(motif_similarity.mtx, k=1)[,1]
motif_similarity.mds[names(motif_similarity.mds)%in%"T_789"]

# Rename from motif id to TF name
io$archR.motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",opts$motif_annotation))
motif2gene.dt <- fread(io$archR.motif2gene)
tmp <- motif2gene.dt$gene; names(tmp) <- motif2gene.dt$motif
names(motif_similarity.mds) <- tmp[names(motif_similarity.mds)]

# Remove duplicated TFs
motif_similarity.mds <- motif_similarity.mds[!duplicated(names(motif_similarity.mds))]

# Define TF order
TF.order <- names(motif_similarity.mds)

# Sanity checks
stopifnot(sum(is.na(TF.order))==0)
stopifnot(sum(duplicated(TF.order))==0)

#####################################
## Dot plot coloured per cell type ##
#####################################

celltypes.to.plot <- unique(rna_chromvar.dt$celltype)

# i <- "Neural_crest"
for (i in celltypes.to.plot) {
  
  to.plot <- rna_chromvar.dt[celltype==i]  %>%
    merge(tf_markers_rna.dt[celltype==i,c("celltype","gene","tf_marker_score_rna")],by=c("gene","celltype")) %>%
    merge(tf_markers_atac.dt[celltype==i,c("celltype","gene","tf_marker_score_atac")],by=c("gene","celltype"))
  
  to.plot %>% .[,tf_marker_score:=tf_marker_score_rna*tf_marker_score_rna] %>% .[,tf_marker_score:=tf_marker_score/max(tf_marker_score)]
  to.plot %>% .[,dot_size:=minmax.normalisation(tf_marker_score)]
    
  TF.order.i <- TF.order[TF.order%in%to.plot$gene]
  to.plot %>% .[,gene:=factor(gene,levels=TF.order.i)]
  
  # to.plot.text <- to.plot[tf_marker_score_rna>=0.85 & tf_marker_score_atac>=0.85]
  to.plot.text <- to.plot[tf_marker_score>=0.75]
  
  # p <- ggplot(to.plot, aes(x=gene, y=tf_marker_score, group=1)) +
  #   geom_point(aes(size=dot_size), stat = 'identity', shape=21, color="black", fill="gray95", alpha=0.9) +
  #   ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
  #   # geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line[celltype%in%to.plot$celltype]) +
  #   scale_size_continuous(range = c(0.05,3)) +
  #   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by=0.15)) +
  #   guides(size="none") +
  #   coord_polar() + 
  #   theme_bw() +
  #   theme(
  #     panel.background=element_rect(fill = c("white")),
  #     panel.grid.major.x=element_blank(),
  #     panel.grid.major.y=element_line(size=rel(1.5)),
  #     axis.title=element_blank(),
  #     axis.text.x=element_blank(),
  #     axis.text.y=element_text(),
  #     axis.ticks.y=element_blank(),
  #     axis.line.x=element_blank()
  #   )
  
  p <- ggplot(to.plot, aes(x=gene, y=tf_marker_score)) +
    geom_point(aes(size=dot_size), stat = 'identity', shape=21, color="black", fill="gray95", alpha=0.9) +
    # ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
    geom_text(data=to.plot.text, aes(label=gene)) +
    scale_size_continuous(range = c(0.05,3.5)) +
    labs(x="TF", y="TF marker score") +
    guides(size="none") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black")
    )
  
  # outfile <- file.path(io$outdir,sprintf("%s_rna_vs_chromvar_chip_%s_polar_plots_pseudobulk.pdf",i,opts$motif_annotation))
  outfile <- file.path(io$outdir,sprintf("%s_rna_vs_chromvar_chip_%s_dot_plots_pseudobulk.pdf",i,opts$motif_annotation))
  # outfile <- file.path(io$outdir,"evidence_legend.pdf")
  pdf(outfile, width = 8, height = 6)
  print(p)
  dev.off()
  
}

#######################################################
## Dot plot per cell type, coloured by TF evidence ##
#######################################################

# celltypes.to.plot <- unique(rna_chromvar.dt$celltype)
celltypes.to.plot <- c("PGC","Neural_crest")

# i <- "Neural_crest"
# i <- "PGC"
for (i in celltypes.to.plot) {
  
  to.plot <- rna_chromvar.dt[celltype==i]  %>%
    merge(tf_markers_rna.dt[celltype==i,c("celltype","gene","tf_marker_score_rna")],by=c("gene","celltype")) %>%
    merge(tf_markers_atac.dt[celltype==i,c("celltype","gene","tf_marker_score_atac")],by=c("gene","celltype")) %>% 
    .[,tf_marker_score:=tf_marker_score_rna*tf_marker_score_rna] %>% .[,tf_marker_score:=tf_marker_score/max(tf_marker_score)] %>% 
    .[,dot_size:=minmax.normalisation(tf_marker_score)]
  
  to.plot <- to.plot %>% 
    merge(tf_evidence.dt,by=c("gene","celltype"), all.x=T) %>%
    .[is.na(evidence),evidence:="None"]
    
  to.plot[,foo:=ifelse(tf_marker_score<=0.75,"foo","bar")]
  to.plot.text <- to.plot[tf_marker_score>=0.75] %>% .[,gene:=str_to_title(gene)]
  
  # p <- ggplot(to.plot, aes(x=gene, y=tf_marker_score, group=1)) +
  #   geom_point(aes(size=dot_size, fill=evidence, alpha=foo), stat = 'identity', shape=21, color="black") +
  #   ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
  #   scale_fill_manual(values=tf_evidence.colors) +
  #   scale_size_continuous(range = c(0.05,6)) +
  #   scale_alpha_manual(values=c(1,0.70)) +
  #   scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by=0.15)) +
  #   coord_polar() + 
  #   theme_bw() +
  #   guides(size="none", alpha="none") +
  #   theme(
  #     legend.position = "right",
  #     panel.background=element_rect(fill = c("white")),
  #     panel.grid.major.x=element_blank(),
  #     panel.grid.major.y=element_line(size=rel(1.5)),
  #     axis.title=element_blank(),
  #     axis.text.x=element_blank(),
  #     axis.text.y=element_text(),
  #     axis.ticks.y=element_blank(),
  #     axis.line.x=element_blank()
  #   )
  
  p <- ggplot(to.plot, aes(x=gene, y=tf_marker_score)) +
    geom_point(aes(size=dot_size, fill=evidence, alpha=foo), stat = 'identity', shape=21, color="black") +
    # ggrepel::geom_text_repel(data=to.plot.text, aes(label=gene), size=3, max.overlaps=Inf) +
    geom_text(data=to.plot.text, aes(label=gene)) +
    scale_fill_manual(values=tf_evidence.colors) +
    scale_size_continuous(range = c(0.05,4)) +
    geom_hline(yintercept=0.75, linetype="dashed") +
    scale_alpha_manual(values=c(1,0.70)) +
    labs(x="TF", y="TF marker score") +
    guides(size="none", alpha="none") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color="black")
    )
  
  outfile <- file.path(io$outdir,sprintf("%s_rna_vs_chromvar_chip_%s_dot_plots_tf_evidence.pdf",i,opts$motif_annotation))
  # outfile <- file.path(io$outdir,"evidence_legend.pdf")
  pdf(outfile, width = 4, height = 4)
  print(p)
  dev.off()
  
}
