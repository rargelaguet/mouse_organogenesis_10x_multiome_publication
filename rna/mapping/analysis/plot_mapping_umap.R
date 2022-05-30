# NOTE THAT THIS IS CURRENTLY IMPLEMENTED ONLY FOR THE MNN MAPPING

suppressPackageStartupMessages(library(argparse))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--atlas_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$query_metadata <- "/Users/argelagr/data/gastrulation_histones/results/rna/mapping/sample_metadata_after_mapping.txt.gz"
# args$atlas_metadata <- "/Users/argelagr/data//gastrulation10x/sample_metadata.txt.gz"
# args$outdir <- "/Users/argelagr/data/gastrulation_histones/results/rna/mapping/pdf"
## END TEST ##

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))

# Options

# Dot size
opts$size.mapped <- 0.18
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.65
opts$alpha.nomapped <- 0.35

###################
## Load metadata ##
###################

sample_metadata <- fread(args$query_metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%args$samples & !is.na(closest.cell)]

stopifnot("closest.cell"%in%colnames(sample_metadata))

####################
## Load 10x atlas ##
####################

# Load atlas cell metadata
meta_atlas <- fread(args$atlas_metadata) %>%
  # .[celltype%in%opts$celltypes] %>%
  .[stripped==F & doublet==F]

# Extract precomputed dimensionality reduction coordinates
umap.dt <- meta_atlas %>%
  .[,c("cell","umapX","umapY","celltype")] %>%
  setnames(c("umapX","umapY"),c("V1","V2"))

##############################
## Define plotting function ##
##############################


plot.dimred <- function(plot_df, query.label, atlas.label = "Atlas") {
  
  # Define dot size  
  size.values <- c(opts$size.mapped, opts$size.nomapped)
  names(size.values) <- c(query.label, atlas.label)
  
  # Define dot alpha  
  alpha.values <- c(opts$alpha.mapped, opts$alpha.nomapped)
  names(alpha.values) <- c(query.label, atlas.label)
  
  # Define dot colours  
  colour.values <- c("red", "lightgrey")
  names(colour.values) <- c(query.label, atlas.label)
  
  # Plot
  ggplot(plot_df, aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = size.values) +
    scale_alpha_manual(values = alpha.values) +
    scale_colour_manual(values = colour.values) +
    # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}


####################
## Plot all cells ##
####################

to.plot <- umap.dt %>% copy %>%
  .[,index:=match(cell, sample_metadata[,closest.cell] )] %>% 
  .[,mapped:=as.factor(!is.na(index))] %>% 
  .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("scRNA-seq atlas","10x Multiome"))] %>%
  setorder(mapped) 

p <- plot.dimred(to.plot, query.label = "10x Multiome", atlas.label = "scRNA-seq atlas")

pdf(sprintf("%s/umap_mapped_allcells.pdf",args$outdir), width=8, height=6.5)
print(p)
dev.off()

###############################
## Plot one sample at a time ##
###############################

for (i in args$samples) {
  
  to.plot <- umap.dt %>% copy %>%
    .[,index:=match(cell, sample_metadata[sample==i,closest.cell] )] %>% 
    .[,mapped:=as.factor(!is.na(index))] %>% 
    .[,mapped:=plyr::mapvalues(mapped, from = c("FALSE","TRUE"), to = c("Atlas",i))] %>%
    setorder(mapped) 
  
  p <- plot.dimred(to.plot, query.label = i, atlas.label = "Atlas")
  
  pdf(sprintf("%s/umap_mapped_%s.pdf",args$outdir,i), width=8, height=6.5)
  print(p)
  dev.off()
}
