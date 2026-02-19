data_folder <- "data"


###########################
## Load global variables ##
###########################

source("R/global.R")

genes <- fread(file.path(data_folder,"genes.txt"), header=F)[[1]]
gene_signatures_names <- fread(file.path(data_folder,"gene_signatures/gene_signatures.txt"), header=F)[[1]]

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(file.path(data_folder,"sample_metadata.txt")) %>% 
  .[,time:=factor(time, levels=timepoints)] %>%
  .[,cell_line:=factor(cell_line, levels=cell_lines)] %>%
  .[,surface_marker:=factor(surface_marker, levels=surface_markers)] %>%
  .[condition_surface_marker=="Ctrl_None",condition_surface_marker:="Ctrl"] %>%
  .[,condition_surface_marker:=factor(condition_surface_marker,levels=c("Ctrl", "Ctrl_CD13", "Dox_CD13", "Dox_CD70", "Dox_SSEA4"))]

sample_metadata.dt[,condition_surface_marker_time:=sprintf("%s_%s_%s",condition,surface_marker,time)]
sample_metadata.dt[,donor:=cell_line]

sample_metadata.dt <- sample_metadata.dt[!surface_marker_time=="None_12ON-14OFF"]
sample_metadata.dt <- sample_metadata.dt[!surface_marker_time=="CD70_22ON"]

####################
## Load gene sets ##
####################

gene_sets <- list()
gene_sets[["transcription_factors"]] <- fread(file.path(data_folder,"gene_sets/transcription_factors.txt"))[[1]]
gene_sets[["surface_markers"]] <- fread(file.path(data_folder,"gene_sets/surface_markers.txt"))[[1]]
gene_sets[["epigenetic_regulators"]] <- fread(file.path(data_folder,"gene_sets/epigenetic_factors.txt"))[[1]]
gene_sets[["dna_repair"]] <- fread(file.path(data_folder,"gene_sets/dna_repair.txt"))[[1]]
gene_sets[["ecm"]] <- fread(file.path(data_folder,"gene_sets/ECM_genes.txt"))[[1]]
gene_sets[["mitochondrial"]] <- fread(file.path(data_folder,"gene_sets/mitochondrial_genes.txt"))[[1]]
gene_sets[["ribosomal"]] <- fread(file.path(data_folder,"gene_sets/ribosome_genes.txt"))[[1]]
gene_sets[["cell_cycle"]] <- fread(file.path(data_folder,"gene_sets/cell_cycle_g2m.txt"))[[1]]
gene_sets[["transposable_elements"]] <- fread(file.path(data_folder,"gene_sets/transposable_elements.txt"))[[1]]

gene_sets_names <- c("Transcription Factors","Surface markers","Cell cycle genes","Epigenetic regulators","DNA repair","ECM genes","Mitochondrial","Ribosomal","Transposable Elements")

##########################
## Load gene expression ##
##########################

expr.sce <- readRDS(file.path(data_folder,"SingleCellExperiment.rds"))[,sample_metadata.dt$sample]
ucell.sce <- readRDS(file.path(data_folder,"gene_signatures/UCell_SingleCellExperiment.rds"))[,sample_metadata.dt$sample]

#################
## Load clocks ##
#################

transcriptional_clocks.dt <- fread(file.path(data_folder,"clocks/transcriptional_clocks.txt"))
transcriptional_clocks_names <- fread(file.path(data_folder,"clocks/transcriptional_clocks_names.txt"),header = F)[[1]]
# epigenetic_clocks.dt <- fread(file.path(data_folder,"clocks/epigenetic_clocks.txt"))

