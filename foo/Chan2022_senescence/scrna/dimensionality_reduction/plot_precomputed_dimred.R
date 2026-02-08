source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"cell_metadata.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction"); dir.create(io$outdir, showWarnings = F)

# Options
opts$samples <- c( 
  "D0",
  "D2",
  "D8",
  "D12",
  "D16_TRA-1-60_negative",
  "D16_TRA-1-60_positive"
)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,sample:=NULL] %>% setnames("alias","sample") %>%
  .[sample%in%opts$samples] %>% .[,sample:=factor(sample,levels=opts$samples)] %>%
  .[,stage:=factor(stage,levels=opts$stages)]

table(cell_metadata.dt$sample,cell_metadata.dt$celltype_liu2020)

###############################################
## Load precomputed dimensionality reduction ##
###############################################

io$dimred.dir <- paste0(io$basedir,"/results/rna/dimensionality_reduction")

pca.dt <- fread(file.path(io$dimred.dir,"precomputed_pca.txt.gz"))# %>% .[,cell:=gsub("scRNA#","",cell)]
umap.dt <- fread(file.path(io$dimred.dir,"precomputed_umap.txt.gz"))# %>% .[,cell:=gsub("scRNA#","",cell)]

umap.dt$V2 <- -umap.dt$V2
umap.dt$V1 <- -umap.dt$V1

##########################################################
## Load precomputed dimensionality reduction per sample ##
##########################################################

# io$dimred.dir <- paste0(io$basedir,"/results/rna/dimensionality_reduction/per_stage")
# 
# pca_per_stage.dt <- opts$samples %>% map(function(i) fread(sprintf("%s/%s/X_pca.txt",io$dimred.dir,i))) %>% rbindlist
# umap_per_stage.dt <- opts$samples %>% map(function(i) fread(sprintf("%s/%s/X_umap.txt",io$dimred.dir,i))) %>% rbindlist
# fa_per_stage.dt <- opts$samples %>% map(function(i) fread(sprintf("%s/%s/X_fa.txt",io$dimred.dir,i))) %>% rbindlist

##########
## Plot ##
##########

to.plot <- umap.dt %>% setnames(c("cell","V1","V2")) %>%
# to.plot <- umap_per_stage.dt %>% setnames(c("cell","V1","V2")) %>%
  # to.plot <- pca_per_stage.dt %>% setnames(c("cell",sprintf("V%d",1:(ncol(pca_per_stage.dt)-1)))) %>%
  # to.plot <- fa_per_stage.dt %>% setnames(c("cell","V1","V2")) %>%
  merge(cell_metadata.dt, by="cell")

p <- ggplot(to.plot, aes(x=V1, y=V2)) +
  ggrastr::geom_point_rast(aes(fill=stage), shape=21, size=0.90, stroke=0.03, color="black", alpha=0.90, raster.dpi=300) +
  guides(fill = guide_legend(override.aes = list(size=2))) +
  # ggrepel::geom_text_repel(aes(label=stage), data=to.plot[,.(V1=median(V1), V2=median(V2)), by=c("stage","dataset")]) +
  # facet_wrap(~stage, nrow=2) +
  # scale_fill_manual(values=opts$stage.colors) +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    legend.position  ="right",
    legend.direction = "vertical",
    legend.title=element_blank()
  )

# pdf(file.path(io$outdir,"Nair2023_UMAP_split_stage_cluster.pdf"), width=9, height=5.5)
pdf(file.path(io$outdir,"Xing2020_UMAP_stage.pdf"), width=7, height=5.5)
print(p)
dev.off()

###############################
## Plot one sample at a time ##
###############################

# for (i in unique(to.plot$sample)) {
#   print(i)
#   
#   to.plot_i <- to.plot %>% copy %>%
#     .[,alpha:=1.0] %>%
#     .[sample!=i,c("sample","alpha"):=list("None",0.25)]
#   
#   p <- ggplot() +
#     ggrastr::geom_point_rast(aes(x=x, y=y), size=1, color="grey", alpha=0.25, data=to.plot_i[sample!=i]) +
#     ggrastr::geom_point_rast(aes(x=x, y=y, fill=sample), size=3, shape=21, stroke=0.15, alpha=0.9, data=to.plot_i[sample==i]) +
#     scale_fill_manual(values=opts$sample.colors) +
#     # viridis::scale_fill_viridis(discrete=TRUE,option="C") +
#     # guides(fill = guide_legend(override.aes = list(size=4))) +
#     labs(x="", y="", title=i) +
#     theme_classic() +
#     ggplot_theme_NoAxes() +
#     theme(
#       plot.title = element_text(size=rel(2.0), hjust = 0.5),
#       legend.title = element_blank(),
#       legend.position="none"
#     )
#   
#   # pdf(sprintf("%s/umap_stage_%s.pdf",io$outdir,i), width=6.5, height=5)
#   png(sprintf("%s/fdl_sample_%s.png",io$outdir,i), width=700, height=500)
#   print(p)
#   dev.off()
# }

###########################
## Cell type annotations ##
###########################

# cell_metadata.dt[celltype=="nic",celltype:=NA]

to.plot <- cell_metadata.dt[,.N,by=c("celltype_liu2020","sample")] %>% 
  .[!is.na(celltype_liu2020)]

ggplot(to.plot, aes(x=sample, y=N, fill=celltype_liu2020)) +
  geom_bar(stat="identity", color="black") +
  # facet_wrap(~cell_line, nrow=1, scales = "free_x") +
  guides(x = guide_axis(angle = 90)) +
  labs(x="", y="Number of cells") +
  theme_classic() +
  # ggplot_theme_NoAxes() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(color="black", size=rel(1)),
    axis.text.x = element_text(color="black", size=rel(0.85)),
    axis.title = element_text(color="black", size=rel(1.10)),
    legend.title = element_blank(),
    legend.position = "right"
  )

#########
## GIF ##
#########

# library(magick)
# 
# io$outdir <- paste0(io$basedir,"/results/dimensionality_reduction/gif"); dir.create(io$outdir, showWarnings = F)
# 
# samplees.to.plot <- c(
#   "D0-fibroblast", "D4-fibroblast", "D8-fibroblast", 
#   "D12-primed", "D16-primed", "D20-primed", "D24-primed", "P20-primed", 
#   "D12-naive", "D16-naive", "D20-naive", "D24-naive", "P3-naive", "P20-naive"
# )
# 
# # 2) Create gif
# files <- sprintf("%s/fdl_sample_%s.png",io$outdir,samplees.to.plot)
# stopifnot(file.exists(files))
# 
# files %>%
#   map(image_read) %>%
#   image_join %>%
#   image_animate(fps=1) %>%
#   image_write(quality=100, path = sprintf("%s/fdl_per_stage.gif",io$outdir))

