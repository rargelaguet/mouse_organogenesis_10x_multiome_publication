here::i_am("rna_atac/gene_regulatory_networks/metacells/build_global_GRN_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

library(cowplot)
library(furrr)

#####################
## Define settings ##
#####################

# I/O
io$basedir <- file.path(io$basedir,"test")
io$trajectory <- file.path("results/rna/trajectories/nmp/nmp_trajectory.txt.gz")
io$sce <- file.path(io$basedir, 'results/rna/metacells/trajectories/nmp/SingleCellExperiment_metacells.rds')
io$tf2gene_virtual_chip <- file.path(io$basedir,"results/rna_atac/virtual_chipseq/metacells/trajectories/nmp/JASPAR/TF2gene_after_virtual_chip.txt.gz")
io$outdir <-  file.path(io$basedir,"results/rna_atac/gene_regulatory_networks/metacells/trajectories/nmp"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
# opts$celltypes <- setdiff(opts$celltypes, c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm"))
opts$trajectory_name <- "nmp"
opts$min_chip_score <- 0.15
opts$max_distance <- 5e4
opts$ncores <- 4

####################################################
## Load TF2gene links based on in silico ChIP-seq ##
####################################################

tf2gene_chip.dt <- fread(io$tf2gene_virtual_chip) %>%
  .[chip_score>=opts$min_chip_score & dist<=opts$max_distance] %>% 
  .[,c("tf","gene")] %>% unique # Only keep TF-gene links

# tf2gene_chip.dt[tf=="T"] %>% View

##############################
## Load RNA expression data ##
##############################

sce <- readRDS(io$sce)

# (Optional) restrict to marker genes
# marker_genes.dt <- fread(io$rna.atlas.marker_genes)
# sce <- sce[rownames(sce)%in%unique(marker_genes.dt$gene),]

##########################
## Filter TFs and genes ##
##########################

TFs <- intersect(unique(tf2gene_chip.dt$tf),toupper(rownames(sce)))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(sce))

tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

# Fetch RNA expression matrices
rna_tf.mtx <- logcounts(sce)[str_to_title(unique(tf2gene_chip.dt$tf)),]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
rna_targets.mtx <- logcounts(sce)[unique(tf2gene_chip.dt$gene),]

# Filter out lowly variable genes and TFs
rna_tf.mtx <- rna_tf.mtx[apply(rna_tf.mtx,1,var)>=1,]
rna_targets.mtx <- rna_targets.mtx[apply(rna_targets.mtx,1,var)>=0.1,]

TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(rna_tf.mtx))
genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna_targets.mtx))
tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

print(sprintf("Number of TFs: %s",length(TFs)))
print(sprintf("Number of genes: %s",length(genes)))

####################
## run regression ##
####################

# cvfit <- cv.glmnet(x, y, alpha=0)
# opt_ridge <- glmnet(x, y, alpha = 0, lambda  = cvfit$lambda.min)
# df <- data.frame(tf=rownames(opt_ridge$beta), gene=i, beta=as.matrix(opt_ridge$beta)[,1])

if (opts$ncores>1){
  plan(multicore, workers=opts$ncores)
  genes_split <- split(genes, cut(seq_along(genes), opts$ncores, labels = FALSE)) 
} else {
  plan(sequential)
  genes_split <- list(genes)
}

# length(genes_split); sapply(genes_split,length)

GRN_coef.dt <- genes_split %>% future_map(function(genes) {
  tmp <- tf2gene_chip.dt[gene%in%genes]
  genes %>% map(function(i) {
    # print(sprintf("%s (%d/%d)",i,match(i,genes),length(genes)))
    tfs <- tmp[gene==i,tf]
    tfs %>% map(function(j) {
      x <- rna_tf.mtx[j,]
      y <- rna_targets.mtx[i,]
      lm.fit <- lm(y~x)
      data.frame(tf=j, gene=i, beta=round(coef(lm.fit)[[2]],3), pvalue=format(summary(lm.fit)$coefficients[2,4], digits=3))
    }) %>% rbindlist
  }) %>% rbindlist %>% return
}) %>% rbindlist

# save
fwrite(GRN_coef.dt, file.path(io$outdir,'global_chip_GRN_coef.txt.gz'), sep="\t")

