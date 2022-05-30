here::i_am("rna/celltype_proportions/compare_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
io$basedir <- file.path(io$basedir,"test")
args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
args$celltype_label <- "celltype"
args$outdir <- file.path(io$basedir,"results/rna/celltype_proportions/comparisons")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)

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
  # "PGC",
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
  # "Mesenchyme",
  "Haematoendothelial_progenitors",
  "Endothelium",
  "Blood_progenitors_1",
  "Blood_progenitors_2",
  "Erythroid1",
  "Erythroid2",
  "Erythroid3",
  "NMP",
  "Rostral_neurectoderm",
  # "Caudal_neurectoderm",
  "Neural_crest",
  "Forebrain_Midbrain_Hindbrain",
  "Spinal_cord",
  "Surface_ectoderm"
  # "Visceral_endoderm"
  # "ExE_endoderm",
  # "ExE_ectoderm"
  # "Parietal_endoderm"
)

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & doublet_call==FALSE & genotype%in%c("WT","T_KO") & celltype%in%opts$celltypes]

# Print statistics
print(table(sample_metadata$sample))
print(table(sample_metadata$celltype))

####################################
## Calculate celltype proportions ##
####################################

# Calculate celltype proportions for WT
wt_proportions.dt <- sample_metadata %>%
  .[genotype=="WT"] %>%
  .[,ncells:=.N] %>%
  .[,.(proportion=.N/unique(ncells), N=round(.N/length(unique(sample)))), by="celltype"]


# Calculate celltype proportions for KO samples
ko_proportions_per_sample.dt <- sample_metadata %>%
  .[genotype=="T_KO"] %>%
  setkey(celltype,sample) %>%
  .[CJ(celltype,sample, unique = TRUE), .N, by = .EACHI] %>%
  merge(unique(sample_metadata[,c("sample","genotype")]), by="sample") %>%
  .[,ncells:=sum(N), by="sample"] %>% .[,proportion:=(N+1)/ncells]

# Calculate celltype proportions for KO genotypes
ko_proportions_per_genotype.dt <- sample_metadata %>%
  .[genotype=="T_KO"] %>%
  setkey(celltype,genotype) %>%
  .[CJ(celltype,genotype, unique = TRUE), .N, by = .EACHI] %>%
  .[,ncells:=sum(N), by="genotype"] %>% .[,proportion:=(N+1)/ncells]

# Merge
proportions_per_sample.dt <- merge(
  ko_proportions_per_sample.dt, 
  wt_proportions.dt, 
  by = c("celltype"), allow.cartesian=T, suffixes = c(".ko",".wt")
) %>% .[,c("diff_proportion"):=list(log2(proportion.ko/proportion.wt))] 

proportions_per_genotype.dt <- merge(
  ko_proportions_per_genotype.dt, 
  wt_proportions.dt, 
  by = c("celltype"), allow.cartesian=T, suffixes = c(".ko",".wt")
) %>% .[,c("diff_proportion"):=list(log2(proportion.ko/proportion.wt))] 

#########################
## Barplots per sample ##
#########################

ylimits <- max(abs(proportions_per_sample.dt[!is.infinite(diff_proportion),diff_proportion]))

# i <- "E8.5_CRISPR_T_KO"
for (i in unique(proportions_per_sample.dt$sample)) {
  
  to.plot <- proportions_per_sample.dt %>%
    .[sample==i] %>% 
    .[N.ko+N.wt>=25]
  
  celltype.order <- to.plot %>%
    .[,mean(diff_proportion),by="celltype"] %>% setorder(-V1) %>% .$celltype
  to.plot <- to.plot %>% .[,celltype:=factor(celltype,levels=celltype.order)]
  
  p <- ggplot(to.plot, aes(x=factor(celltype), y=diff_proportion)) +
    geom_point(aes(fill = celltype), shape=21, size=1) +
    geom_bar(aes(fill = celltype), stat="identity", alpha=0.5, color="black") +
    # geom_text(y=-ylim, aes(label=N.wt), size=2.5) +
    # geom_text(y=ylim, aes(label=N.ko), size=2.5) +
    coord_flip(ylim=c(-ylimits,ylimits)) +
    geom_hline(yintercept=0, linetype="dashed", size=0.5) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    theme_classic() +
    labs(y="Difference in proportions (log2)", x="") +
    theme(
      legend.position = "none",
      # axis.title = element_blank(),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black")
    )
  
  pdf(sprintf("%s/%s_barplots.pdf",args$outdir,i))
  print(p)
  dev.off()
}

###########################
## Boxplots per genotype ##  
###########################

# ylimits <- max(abs(proportions_per_sample.dt$diff_proportion)) + 0.25
# 
# for (i in unique(proportions_per_sample.dt$genotype)) {
#   
#   celltypes.to.plot <- proportions_per_sample.dt %>%
#     .[genotype==i,.(N=sum(N.ko)+sum(N.wt)),by=c("genotype","celltype")] %>% 
#     .[N>=50,celltype] %>% as.character
#     # .[N>=5,celltype] %>% as.character
#   
#   to.plot <- proportions_per_sample.dt %>%
#     .[genotype==i & celltype%in%celltypes.to.plot] %>%
#     merge(unique(sample_metadata[,c("sample")]),by="sample")
#   
#   celltype.order <- to.plot %>%
#     .[,mean(diff_proportion),by="celltype"] %>% setorder(-V1) %>% .$celltype
#   to.plot <- to.plot %>% .[,celltype:=factor(celltype,levels=celltype.order)]
#   
#   # text.dt <- proportions_per_genotype.dt %>% 
#   #   .[genotype==i & celltype%in%celltype.order] %>% 
#   #   .[,celltype:=factor(celltype,levels=celltype.order)]
#   
#   p <- ggplot(to.plot, aes(x=celltype, y=diff_proportion)) +
#     geom_point(aes(fill = celltype), shape=21, size=1) +
#     geom_boxplot(aes(fill = celltype), alpha=0.5) +
#     # geom_text(y=-ylimits, aes(label=N.wt), size=2.5, data=text.dt) +
#     # geom_text(y=ylimits, aes(label=N.ko), size=2.5, data=text.dt) +
#     coord_flip(ylim=c(-ylimits,ylimits)) +
#     geom_hline(yintercept=0, linetype="dashed", size=0.5) +
#     scale_fill_manual(values=opts$celltype.colors, drop=F) +
#     theme_classic() +
#     labs(y="Difference in proportions (log2)", x="", title=i) +
#     theme(
#       legend.position = "none",
#       # axis.title = element_blank(),
#       plot.title = element_text(size=rel(1.25), hjust=0.5, color="black"),
#       axis.text.y = element_text(color="black"),
#       axis.text.x = element_text(color="black")
#     )
#   
#   pdf(sprintf("%s/%s_boxplots.pdf",args$outdir,i), width=9, height=7)
#   print(p)
#   dev.off()
# }

############################
## Polar plots per sample ##
############################

to.plot.wt_line <- data.table(
  celltype = unique(proportions_per_sample.dt$celltype),
  diff_proportion = log2(1)
)

# ylimits <- max(abs(proportions_per_sample.dt[!is.infinite(diff_proportion),diff_proportion]))
ylimits <- 6

for (i in unique(proportions_per_sample.dt$sample)) {
  
  to.plot <- proportions_per_sample.dt %>%
    .[sample==i] %>% 
    .[N.ko+N.wt>=25] %>%
    .[diff_proportion>=ylimits,diff_proportion:=ylimits] %>%
    .[diff_proportion<=(-ylimits),diff_proportion:=(-ylimits)]
  
  p <- ggplot(to.plot, aes(x=celltype, y=-diff_proportion)) +
    geom_jitter(aes(fill = celltype), size=5, shape=21, alpha=0.9, width=0.05, height=0.1) +
    # geom_boxplot(aes(fill = celltype)) +
    geom_polygon(group=1, color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line[celltype%in%to.plot$celltype]) +
    scale_fill_manual(values=opts$celltype.colors, drop=F) +
    coord_polar() + ylim(ylimits,-ylimits) +
    guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
    scale_size(guide = 'none') +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text =  element_text(size=rel(0.5)),
      legend.text = element_text(size=rel(0.75)),
      legend.title = element_blank(),
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.line=element_blank(),
      axis.text.x = element_blank()
    )
  
  pdf(sprintf("%s/%s_polar_plot.pdf",args$outdir,i), width = 5, height=5)
  print(p)
  dev.off()
}


##############################
## Polar plots per genotype ##
##############################

# # ylimits <- max(abs(proportions_per_genotype.dt[!is.infinite(diff_proportion),diff_proportion]))
# ylimits <- 6
# 
# for (i in unique(proportions_per_genotype.dt$genotype)) {
#   
#   celltypes.to.plot <- proportions_per_sample.dt %>%
#     .[genotype==i,.(N=sum(N.ko)+sum(N.wt)),by=c("genotype","celltype")] %>% 
#     .[N>=50,celltype] %>% as.character
#   
#   to.plot <- proportions_per_sample.dt %>%
#     .[genotype==i & celltype%in%celltypes.to.plot]  %>%
#     .[diff_proportion>=ylimits,diff_proportion:=ylimits] %>%
#     .[diff_proportion<=(-ylimits),diff_proportion:=(-ylimits)]
#   
#   p <- ggplot(to.plot, aes(x=factor(celltype), y=-diff_proportion, group=1)) +
#     geom_point(aes(fill = celltype), size=3, stat = 'identity', shape=21) +
#     geom_polygon(color="black", fill=NA, alpha=0.5, linetype="dashed", data=to.plot.wt_line[celltype%in%to.plot$celltype]) +
#     scale_fill_manual(values=opts$celltype.colors, drop=F) +
#     coord_polar() + ylim(ylimits,-ylimits) +
#     guides(colour = guide_legend(override.aes = list(size=2), ncol=1)) +
#     scale_size(guide = 'none') +
#     theme_bw() +
#     theme(
#       legend.position = "none",
#       strip.text =  element_text(size=rel(0.5)),
#       legend.text = element_text(size=rel(0.75)),
#       legend.title = element_blank(),
#       axis.title=element_blank(),
#       axis.text.y=element_blank(),
#       axis.ticks.y=element_blank(),
#       axis.line=element_blank(),
#       axis.text.x = element_blank()
#     )
#   
#   pdf(sprintf("%s/polar_plots/per_genotype/%s_polar_plot.pdf",args$outdir,i), width = 5, height=5)
#   print(p)
#   dev.off()
# }

# Completion token
#file.create(file.path(args$outdir,"completed.txt"))