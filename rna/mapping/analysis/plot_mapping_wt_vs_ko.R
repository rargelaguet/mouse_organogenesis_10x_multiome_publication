# here::i_am("mapping/analysis/plot_mapping_umap.R")

source(here::here("settings.R"))
source(here::here("rna/mapping/analysis/plot_utils.R"))

io$query_metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$atlas_metadata <- file.path(io$atlas.basedir,"sample_metadata.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/mapping/pdf/fig"); dir.create(io$outdir, showWarnings = F)

#####################
## Define settings ##
#####################

# Options
opts$remove_ExE_cells <- FALSE
opts$subset_atlas <- TRUE

# Dot size
opts$size.mapped <- 0.22
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.75
opts$alpha.nomapped <- 0.35

opts$samples <- c(
  "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

#########################
## Load query metadata ##
#########################

sample_metadata <- fread(io$query_metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%opts$samples & !is.na(closest.cell)] %>%
  .[,sample:=factor(sample,levels=opts$samples)]

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>% .[!celltype%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

################
## Load atlas ##
################

# Load atlas cell metadata
meta_atlas <- fread(io$atlas_metadata) %>%
  .[stripped==F & doublet==F]

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  meta_atlas <- meta_atlas %>% .[!celltype%in%c("ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

# Subset cells that are not neighbours
if (opts$subset_atlas) {
  meta_atlas <- rbind(
    meta_atlas[cell%in%sample_metadata$closest.cell],
    meta_atlas[!cell%in%sample_metadata$closest.cell][sample.int(50000)]
  )
  
}

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

###############################
## Plot one sample at a time ##
###############################

# to.plot <- opts$samples %>% map(function(i) {
#   # sample_metadata.subset <- sample_metadata[alias==i] %>% .[sample.int(min(ncells.to.subset,nrow(.)))]
#   sample_metadata.subset <- sample_metadata[alias==i]
#   umap.dt %>% copy %>%
#     .[,index:=match(cell, sample_metadata.subset[alias==i,closest.cell] )] %>%
#     .[,mapped:=as.factor(!is.na(index))] %>%
#     .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas","Query"))] %>%
#     .[,sample:=factor(i,levels=opts$samples)] %>%
#     setorder(mapped)
# }) %>% rbindlist
# 
# p <- plot.dimred(to.plot, query.label = "Query", atlas.label = "Atlas") +
#   facet_wrap(~sample, nrow=3) +
#   theme(
#     legend.position = "none",
#     strip.background = element_blank(),
#     strip.text = element_text(color="black", size=rel(0.8))
#   )
# 
# pdf(file.path(io$outdir,"umap_mapped_E8.5_samples.pdf"), width=8, height=12)
# print(p)
# dev.off()

#############################
## Plot WT and KO together ##
#############################

opts$wt.label <- "WT"
opts$ko.label <- "T_KO"

# Subsample query cells to have the same N per class
sample_metadata_subset <- sample_metadata %>% .[,.SD[sample.int(n=.N, size=6000)], by=c("stage","genotype")]

to.plot <- umap.dt %>% copy %>%
  .[,index.wt:=match(cell, sample_metadata[genotype==opts$wt.label,closest.cell] )] %>%
  .[,index.ko:=match(cell, sample_metadata[genotype==opts$ko.label,closest.cell] )] %>%
  .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
  .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
  .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
  .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas",opts$ko.label,opts$wt.label))] %>% setorder(mapped)

p <- plot.dimred.wtko(to.plot, wt.label = opts$wt.label, ko.label = opts$ko.label, nomapped.label = "Atlas") +
  theme(legend.position = "top", axis.line = element_blank())

pdf(sprintf("%s/umap_mapped_WT_and_KO.pdf",io$outdir), width=6, height=5)
print(p)
dev.off()
