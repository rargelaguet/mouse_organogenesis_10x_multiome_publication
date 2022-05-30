here::i_am("rna/celltype_proportions/compare_celltype_proportions_paga.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
io$outdir <- file.path(io$basedir,"results/rna/celltype_proportions/comparisons"); dir.create(io$outdir, showWarnings = F)

####################
## Define options ##
####################

opts$samples <- c(
  "E8.5_CRISPR_T_KO",
  "E8.5_CRISPR_T_WT"
)

opts$celltypes <- c(
  "Epiblast",
  "Primitive_Streak",
  "Caudal_epiblast",
  "PGC",
  "Anterior_Primitive_Streak",
  "Notochord",
  "Def._endoderm",
  "Gut",
  "Nascent_mesoderm",
  "Mixed_mesoderm",
  "Intermediate_mesoderm",
  "Caudal_Mesoderm",
  "Paraxial_mesoderm",
  "Somitic_mesoderm",
  "Pharyngeal_mesoderm",
  "Cardiomyocytes",
  "Allantois",
  "ExE_mesoderm",
  "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm",
  "Visceral_endoderm",
  "ExE_endoderm",
  "ExE_ectoderm",
  "Parietal_endoderm"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & sample%in%opts$samples & celltype%in%opts$celltypes]

# Print statistics
print(table(sample_metadata$sample))
print(table(sample_metadata$celltype))

# Filter cell types with small N
tmp <- sample_metadata[,.N,by=c("celltype","genotype")] %>% 
  dcast(celltype~genotype,value.var="N", fill=0)
celltypes.to.use <- tmp[T_KO+WT>=50,celltype]
sample_metadata <- sample_metadata[celltype%in%celltypes.to.use]

####################################
## Calculate celltype proportions ##
####################################

# Calculate celltype proportions for WT
wt_proportions.dt <- sample_metadata %>%
  .[genotype=="WT"] %>%
  .[,ncells:=.N] %>%
  .[,.(proportion=.N/unique(ncells), N=round(.N/length(unique(sample)))), by="celltype"]

# Calculate celltype proportions for KO
ko_proportions_per_genotype.dt <- sample_metadata %>%
  .[genotype=="T_KO"] %>%
  setkey(celltype,genotype) %>%
  .[CJ(celltype,genotype, unique = TRUE), .N, by = .EACHI] %>%
  .[,ncells:=sum(N), by="genotype"] %>% .[,proportion:=(N+1)/ncells]

# Merge
proportions_per_genotype.dt <- merge(
  ko_proportions_per_genotype.dt, 
  wt_proportions.dt, 
  by = c("celltype"), allow.cartesian=T, suffixes = c(".ko",".wt")
) %>% .[,c("diff_proportion"):=list(log2(proportion.ko/proportion.wt))] 

#####################
## Load PAGA graph ##
#####################

source(here::here("load_paga_graph.R"))

cellype.order <- rownames(connectivity.mtx)

# Add differential abundance values
diff.values <- rep(0,length(opts$celltypes)); names(diff.values) <- opts$celltypes
diff.values[proportions_per_genotype.dt$celltype] <- proportions_per_genotype.dt$diff_proportion
diff.values[diff.values>=2.5] <- 2.5
diff.values[diff.values<=(-5)] <- (-5)

igraph.paga.tbl <- igraph.paga.tbl %>% activate(nodes) %>% mutate(diff=diff.values[cellype.order])

##########
## Plot ##
##########
  
p <- ggraph(igraph.paga.tbl, x = x, y = y) +
  # geom_edge_link0(aes(width = weight), edge_colour = "grey66", edge_alpha=0.2) +
  geom_edge_link(edge_colour = "grey66", edge_alpha=1, edge_width=0.75) +
  geom_node_point(aes(fill=diff), size=8, shape=21) +
  scale_fill_brewer(palette = "PuOr") +
  # scale_fill_manual(values=opts$celltype.colors) +
  # geom_node_point(aes(fill = colour_by, size = size), shape=21, stroke=node_stroke) +
  scale_fill_gradient2() +
  # scale_edge_width(range = c(0.2,3), name="overlap size") +
  theme_classic(base_size=14) +
  theme(
    axis.line = element_blank(), 
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    axis.title = element_blank(),
    legend.position = "right"
  )

pdf(file.path(io$outdir,"paga_coloured_by_diff_abundance_legend.pdf"), width=4, height=4)
print(p)
dev.off()