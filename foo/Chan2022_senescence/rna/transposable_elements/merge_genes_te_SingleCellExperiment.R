source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")

# I/O (SalmonTE)
# io$te.sce <- file.path(io$basedir,"processed/salmonTE/SalmonTE_SingleCellExperiment.rds")
# io$gene.sce <- file.path(io$basedir,"processed/salmon/salmon_genes_SingleCellExperiment.rds")
# io$outfile <- file.path(io$basedir,"processed/salmonTE/salmonTE_TE_genes_SingleCellExperiment.rds")

# I/O (TEtranscripts)
io$te.sce <- file.path(io$basedir,"processed/rna/TEtranscripts/TEtranscripts_TE_SingleCellExperiment.rds")
io$gene.sce <- file.path(io$basedir,"processed/rna/SingleCellExperiment.rds")
io$outfile <- file.path(io$basedir,"processed/rna/TEtranscripts/TEtranscripts_TE_genes_SingleCellExperiment.rds")

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata)

###############################
## Load SingleCellExperiment ##
###############################

# Genes
genes.sce <-  readRDS(io$gene.sce)[,sample_metadata.dt$sample]
# rownames(genes.sce) <- rowData(genes.sce)$gene_name
# colnames(genes.sce) <- sample_metadata.dt$sample
# genes.sce <- genes.sce %>% as("SingleCellExperiment")

# Rename to gene symbols
tmp <- rowData(genes.sce)$symbol==""
rowData(genes.sce)$symbol[rowData(genes.sce)$symbol==""] <- rownames(genes.sce)[tmp]
genes.sce <- genes.sce[!duplicated(rowData(genes.sce)$symbol),]
rownames(genes.sce) <- rowData(genes.sce)$symbol

# TE
te.sce <- readRDS(io$te.sce)[,sample_metadata.dt$sample]# %>% as("SingleCellExperiment")
rowData(te.sce)$name <- rownames(te.sce)

############
## Filter ##
############

# filter samples
samples <- intersect(colnames(genes.sce),colnames(te.sce))
genes.sce <- genes.sce[,samples]; te.sce <- te.sce[,samples] 
sample_metadata.dt <- sample_metadata.dt[sample%in%samples] %>% setkey(sample) %>% .[samples]

# add colData
colData(genes.sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame
colData(te.sce) <- sample_metadata.dt %>% tibble::column_to_rownames("sample") %>% DataFrame
stopifnot(colnames(genes.sce)==colnames(te.sce))

stopifnot(sum(duplicated(rownames(genes.sce)))==0)

###########
## Merge ##
###########

# prepare rowData
foo <- rowData(genes.sce)[,c("ID"),drop=F]; colnames(foo) <- c("gene_name")
bar <- rowData(te.sce)[,c("name"),drop=F]; colnames(bar) <- c("gene_name")
df <- rbind(foo,bar) %>% DataFrame
df$transcript_type <- c(rep("gene",nrow(genes.sce)),rep("te",nrow(te.sce)))

# ensembl ID
# df$ens_id <- rownames(df)
# df$ens_id[df$transcript_type%in%c("te")] <-  NA

# merge
sce <- SingleCellExperiment::rbind(genes.sce,te.sce)

# add rowData
# rownames(sce) <- df$gene_name
rowData(sce) <- df

# Filter duplicated gene names
sce <- sce[!duplicated(rownames(sce)),]

# save
saveRDS(sce, io$outfile)
