# here::i_am("rna/processing/multiome/1_create_seurat_rna.R")

source(here::here("settings.R"))

# suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(scran))

#####################
## Define settings ##
#####################

# I/O
io$inputdir <- paste0(io$basedir,"/cellranger")
io$outdir.data <- paste0(io$basedir,"/processed/rna"); dir.create(io$outdir.data, recursive = T, showWarnings = F)
io$outdir.qc <- paste0(io$basedir,"/results/rna/qc"); dir.create(io$outdir.qc, recursive = T, showWarnings = F)

# Options
# opts$samples <- opts$samples[1]
opts$trim.barcode <- FALSE
opts$sample_cell_separator <- "#"
opts$min_ngenes <- 1500
opts$max_ngenes <- 10000
opts$mitochondrial_percent_RNA <- 25
opts$ribosomal_percent_RNA <- 40

###############
## Load data ##
###############

count_mtx <- list()
cell.info <- list()

# i <- opts$samples[1]
for (i in opts$samples) {
  print(i)
    
  # Load cell metadata
  barcode.loc <- sprintf("%s/%s/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",io$inputdir,i)
  cell.info[[i]] <- fread(barcode.loc, header=F) %>% setnames("barcode") %>%
    .[,barcode:=ifelse(rep(opts$trim.barcode,.N),gsub("-1","",barcode),barcode)] %>%
    .[,c("sample","cell"):=list(i,sprintf("%s%s%s",i,opts$sample_cell_separator,barcode))]
  dim(cell.info[[i]])
  
  # Load matrix  
  count_mtx[[i]] <- Seurat::Read10X(sprintf("%s/%s/outs/filtered_feature_bc_matrix",io$inputdir,i))#[["Gene Expression"]]
  
  # Basic filtering
  # count_mtx[[i]] <- count_mtx[[i]][,colSums(count_mtx[[i]])>=opts$min_counts]

  barcodes <- intersect(colnames(count_mtx[[i]]), cell.info[[i]]$barcode)
  count_mtx[[i]] <- count_mtx[[i]][,barcodes]
  cell.info[[i]] <- cell.info[[i]] %>% .[barcode%in%barcodes] %>% setkey(barcode) %>% .[barcodes]
  stopifnot(cell.info[[i]]$barcode==colnames(count_mtx[[i]]))
  
}

print(lapply(count_mtx,dim))

#######################
## Keep common genes ##
#######################

genes <- Reduce("intersect",lapply(count_mtx,rownames))
for (i in 1:length(count_mtx)) {
  count_mtx[[i]] <- count_mtx[[i]][genes,]
}

stopifnot(length(unique(lapply(count_mtx,nrow)))==1)
stopifnot(length(unique(lapply(count_mtx,rownames)))==1)

#################
## Concatenate ##
#################

# Concatenate cell metadata
cell.info <- rbindlist(cell.info)
rownames(cell.info) <- cell.info$cell

# Concatenate matrices
count_mtx <- do.call("cbind",count_mtx)
colnames(count_mtx) <- cell.info$cell

# Remove duplicated cells
# count_mtx <- count_mtx[,!duplicated(colnames(count_mtx))]

##################
## Filter genes ##
##################

# Remove duplicated genes
count_mtx <- count_mtx[!duplicated(rownames(count_mtx)),]

# Keep genes that are found in the gene metadata file
gene_metadata.dt <- fread(io$gene_metadata) %>% .[!duplicated(symbol)]
genes.to.use <- intersect(rownames(count_mtx),gene_metadata.dt$symbol)
count_mtx <- count_mtx[genes.to.use,]
gene_metadata.dt <- gene_metadata.dt[symbol%in%genes.to.use] %>% setkey(symbol) %>% .[genes.to.use]

stopifnot(gene_metadata.dt$symbol==rownames(count_mtx))

########################################
## Create SingleCellExperiment object ##
########################################

# Sanity checks
stopifnot(sum(duplicated(rownames(count_mtx)))==0)
stopifnot(sum(duplicated(colnames(count_mtx)))==0)
stopifnot(colnames(count_mtx)==cell.info$cell)

# prepare colData
# coldata <- cell.info[cell%in%colnames(count_mtx)] %>% setkey(cell) %>% .[colnames(count_mtx)] %>% tibble::column_to_rownames("cell") %>% DataFrame
coldata <- cell.info %>% as.data.frame %>% tibble::column_to_rownames("cell") %>% DataFrame

# prepare rowData
rowdata <- gene_metadata.dt[,c("symbol","ens_id")] %>% as.data.frame %>% tibble::column_to_rownames("symbol") %>% DataFrame

# create SingleCellExperiemnt
sce <- SingleCellExperiment(
    assays = list(counts=count_mtx),
    rowData = rowdata,
    colData = coldata
    )
sce

#####################
## Quality control ##
#####################

# Number of UMIs
sce[["nUMI_RNA"]] <- colSums(assay(sce,"counts"))

# Add mitochondrial percenatge
mt.genes <- grep("^MT-",rownames(sce),value=T)
sce[["mitochondrial_percent_RNA"]] <- 100*(colSums(assay(sce[mt.genes],"counts")) / sce[["nUMI_RNA"]]) %>% round(4)

# Add ribosomal RNA content
ribo.genes <- grep("^RP[L|S]",rownames(sce),value=T)
sce[["ribosomal_percent_RNA"]] <- 100*(colSums(assay(sce[ribo.genes],"counts")) / sce[["nUMI_RNA"]]) %>% round(4)

# Number of expressed genes
sce[["nFeature_RNA"]] <- colSums(assay(sce,"counts")>0)

# Extract cell metadata
cell_metadata.dt <- colData(sce) %>% as.data.table(keep.rownames=T) %>% setnames("rn","cell")

# Update cell metadata
cell_metadata.dt[,alias:=opts$sample2alias[sample]]
table(cell_metadata.dt$alias)
stopifnot(!any(is.na(cell_metadata.dt$alias)))

cell_metadata.dt[,stage:=opts$sample2stage[sample]]
table(cell_metadata.dt$stage)
stopifnot(!any(is.na(cell_metadata.dt$stage)))

# QC calls
cell_metadata.dt %>% .[,pass_qc:=nFeature_RNA<=opts$max_ngenes & nFeature_RNA>=opts$min_ngenes & mitochondrial_percent_RNA<opts$mitochondrial_percent_RNA & ribosomal_percent_RNA<opts$ribosomal_percent_RNA]
table(cell_metadata.dt$pass_qc)

#####################
## Plot QC metrics ##
#####################

to.plot <- cell_metadata.dt[pass_qc==TRUE] %>%
    melt(id.vars=c("alias","cell"), measure.vars=c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))

facet.labels <- c("nFeature_RNA" = "Num. of genes", "mitochondrial_percent_RNA" = "Mitochondrial %", "ribosomal_percent_RNA" = "Ribosomal %")
    
# Box plot 
p <- ggplot(to.plot, aes(x=alias, y=value, fill=alias)) +
    geom_boxplot(outlier.shape=NA, coef=1) +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        axis.text.x = element_text(colour="black",size=rel(0.50), angle=20, hjust=1, vjust=1),
        axis.title.x = element_blank()
    )

pdf(sprintf("%s/qc_metrics_boxplot.pdf",io$outdir.qc), width=9, height=5)
print(p)
dev.off()

# histogram 
to.plot <- cell_metadata.dt %>% 
  melt(id.vars=c("alias","cell"), measure.vars=c("nFeature_RNA","mitochondrial_percent_RNA","ribosomal_percent_RNA"))
tmp <- data.table(
    variable = c("nFeature_RNA", "mitochondrial_percent_RNA", "ribosomal_percent_RNA"),
    value = c(opts$min_ngenes, opts$mitochondrial_percent_RNA, opts$ribosomal_percent_RNA)
)
p <- gghistogram(to.plot, x="value", fill="alias", bins=50) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free", nrow=1) +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.75))
    )
    
pdf(sprintf("%s/qc_metrics_histogram.pdf",io$outdir.qc), width=13, height=6)
print(p)
dev.off()


# Barplot fraction of cells that pass QC for each alias
to.plot <- cell_metadata.dt %>%
    .[,mean(pass_qc,na.rm=T),by=c("alias")]

p <- ggbarplot(to.plot, x="alias", y="V1", fill="alias") +
    labs(x="", y="Fraction of cells that pass QC (RNA)") +
    theme(
        legend.position = "none",
        axis.text.y = element_text(colour="black",size=rel(0.8)),
        axis.text.x = element_text(colour="black",size=rel(0.50), angle=20, hjust=1, vjust=1),
    )

pdf(sprintf("%s/qc_metrics_barplot.pdf",io$outdir.qc), width=6, height=5)
print(p)
dev.off()

#################################
## Update SingleCellExperiment ##
#################################

cell_metadata_qc.dt <- cell_metadata.dt[pass_qc==TRUE] %>% .[,pass_qc:=NULL]

sce <- sce[,cell_metadata_qc.dt$cell]
colData(sce) <- cell_metadata_qc.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

##########################
## Compute size factors ##
##########################

# clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100, BPPARAM = mcparam))
# min.clust <- min(table(clusts))/2
# new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
# sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

clusts <- as.numeric(quickCluster(sce))
sce <- computeSumFactors(sce, clusters = clusts)

##########
## Save ##
##########

fwrite(cell_metadata_qc.dt, file.path(io$outdir.data,"cell_metadata.txt.gz"), quote=F, na="NA", sep="\t")
saveRDS(sce, file.path(io$outdir.data,"SingleCellExperiment.rds"))

