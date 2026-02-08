source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$cell_cycle_g2m_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/cell_cycle_g2m.txt"
io$cell_cycle_s_genes <- "/Users/rargelaguet/data/hg38_regulation/gene_lists/cell_cycle_s.txt"
io$outdir <- file.path(io$basedir,"rna/results/cell_cycle"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata)

#########################
## Load RNA expression ##
#########################

sce <- readRDS(io$sce)[,sample_metadata.dt$sample]

colData(sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame
rownames(sce) <- rowData(sce)$gene_name

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

cell_cycle_s_genes <- fread(io$cell_cycle_s_genes, header=F)[[1]]
cell_cycle_g2m_genes <- fread(io$cell_cycle_g2m_genes, header=F)[[1]]

cell_cycle_s_genes <- cell_cycle_s_genes[cell_cycle_s_genes%in%rownames(sce)]
cell_cycle_g2m_genes <- cell_cycle_g2m_genes[cell_cycle_g2m_genes%in%rownames(sce)]

###################
## Normalisation ##
###################

assay(sce,"logcounts") <- log2(1e6*(sweep(assay(sce,"counts"),2,colSums(counts(sce)),"/"))+1)

########################
## Cell cycle scoring ##
########################

cell_cycle_scores.dt <- data.table(
  sample = colnames(sce),
  s_score = scale(colMeans(assay(sce,"logcounts")[cell_cycle_s_genes,]))[,1] %>% round(3),
  g2m_score = scale(colMeans(assay(sce,"logcounts")[cell_cycle_g2m_genes,]))[,1] %>% round(3)
)

to.plot <- sample_metadata.dt %>% 
  merge(cell_cycle_scores.dt, by="sample")

to.plot$expr <- assay(sce,"logcounts")["MKI67",][to.plot$sample]

p <- ggplot(to.plot, aes(x=g2m_score, y=expr)) +
  geom_point(aes(color=expr), size=3, alpha=0.9) +
  facet_wrap(~condition, nrow=2, scales="free_x") +
  scale_color_viridis() +
  labs(x="G2M score", y="MKI67 expression") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text = element_text(colour="black",size=rel(0.90)),
    axis.title = element_text(colour="black",size=rel(1.0)),
    legend.position = "right"
  )


pdf(file.path(io$outdir,"cell_cycle_scores_vs_mki67_hdf_reik.pdf"), width=8, height=6)
print(p)
dev.off()

#########
## PCA ##
#########

# sce.filt <- sce[cell_cycle_g2m_genes,]
# sce.filt <- runPCA(sce.filt, ncomponents=2)
# sce.filt$g2m_score <- cell_cycle_scores.dt %>% setkey(sample) %>% .[colnames(sce)] %>% .$g2m_score
# plotPCA(sce.filt, colour_by="MKI67", ncomponents=c(1,2), point_alpha=1, point_size=5)
# plotPCA(sce.filt, colour_by="g2m_score", ncomponents=c(1,2), point_alpha=1, point_size=5)

###################################
## Plot cell cycle score vs age ##
###################################

to.plot <- sample_metadata.dt %>% 
  merge(cell_cycle_scores.dt, by="sample")

to.plot <- to.plot[condition%in%c("hTERT","Replicative_Senescence","Radiation_Senescence")]
  
to.plot.means <- to.plot[,.(g2m_score=mean(g2m_score), sd=sd(g2m_score)), by=c("condition","time")]

ylim.max <- max(to.plot$g2m_score)


p <- ggplot(to.plot, aes(x=time, y=g2m_score)) +
  geom_line(aes(group=condition), color="gray40", linewidth=1, alpha=0.75, data=to.plot.means) +
  geom_point(color="gray40", size=2.50) +
  geom_errorbar(aes(ymin=g2m_score-sd, ymax=g2m_score+sd), width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
  facet_wrap(~condition, nrow=1, scales="free_x") +
  labs(x="Time", y="G2/M score", title="") +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(1.10)),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour="black",size=rel(1.0)),
    axis.title = element_text(colour="black",size=rel(1.10)),
    axis.text.y = element_text(colour="black",size=rel(1.10)),
    legend.position = "none",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,"cell_cycle_scores_vs_age_chan2022_senescence.pdf"), width=10, height=5)
print(p)
dev.off()

##########
## Save ##
##########

to.save <- sample_metadata.dt %>% merge(cell_cycle_scores.dt[,c("sample","g2m_score")] %>% setnames("g2m_score","cell_cycle_score"), by="sample")
fwrite(to.save, file.path(io$outdir,"sample_metadata_with_cell_cycle_scores.txt"), sep="\t", quote=F, na="NA")
