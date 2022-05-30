here::i_am("rna/mapping/trajectories/plot_mapping_dimred.R")

source(here::here("settings.R"))
source(here::here("rna/mapping/analysis/plot_utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--query_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--atlas_metadata',        type="character",                               help='Cell metadata (after mapping)')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
args$query_metadata <- file.path(io$basedir,"results/rna/mapping/trajectories/nmp_somitic_spinal/sample_metadata_after_mapping.txt.gz")
args$atlas_metadata <- file.path(io$atlas.basedir,"results/trajectories/nmp_somitic_spinal/nmp_trajectory.txt.gz")
args$outdir <- file.path(io$basedir,"results/rna/mapping/trajectories/nmp_somitic_spinal/pdf")
## END TEST ##

dir.create(args$outdir, showWarnings = F, recursive = T)

#####################
## Define settings ##
#####################

# Options

# Dot size
opts$size.mapped <- 1.20
opts$size.nomapped <- 0.1

# Transparency
opts$alpha.mapped <- 0.80
opts$alpha.nomapped <- 0.35

#########################
## Load query metadata ##
#########################

sample_metadata <- fread(args$query_metadata) %>%
  .[!is.na(closest.cell)]

stopifnot("closest.cell"%in%colnames(sample_metadata))

###########################
## Load atlas trajectory ##
###########################

meta_atlas <- fread(args$atlas_metadata) %>%
  setnames(c("cell","V1","V2"))

##########
## Plot ##
##########

# i <- "E7.5"
to.plot <- meta_atlas %>% copy %>%
  .[,index.wt:=match(cell, sample_metadata[genotype=="WT",closest.cell] )] %>%
  .[,index.ko:=match(cell, sample_metadata[genotype=="T_KO",closest.cell] )] %>%
  .[,mapped.wt:=c(0,-10)[as.numeric(as.factor(!is.na(index.wt)))]] %>%
  .[,mapped.ko:=c(0,10)[as.numeric(as.factor(!is.na(index.ko)))]] %>%
  .[,mapped:=factor(mapped.wt + mapped.ko, levels=c("0","-10","10"))] %>%
  .[,mapped:=plyr::mapvalues(mapped, from = c("0","-10","10"), to = c("Atlas","WT","T_KO"))] %>% setorder(mapped)

p <- plot.dimred.wtko(to.plot, wt.label = "WT", ko.label = "T_KO", nomapped.label = "Atlas") +
  theme(legend.position = "none", axis.line = element_blank())

pdf(file.path(args$outdir,"umap_mapped_trajectory_WT_and_KO.pdf"), width=4.5, height=5)
print(p)
dev.off()
