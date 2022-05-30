here::i_am("rna/celltype_proportions/plot_celltype_proportions.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata',        type="character",                               help='Cell metadata file')
p$add_argument('--samples',         type="character",       nargs="+",   help='Samples')
p$add_argument('--celltype_label', type="character", help='Cell type label')
p$add_argument('--outdir',          type="character",                               help='Output file')

args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args$metadata <- file.path(io$basedir,"results/rna/mapping/sample_metadata_after_mapping.txt.gz")
# args$samples <- opts$samples
# args$celltype_label <- "celltype.mapped"
# args$outdir <- file.path(io$basedir,"results/rna/celltype_proportions")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings = F)
dir.create(file.path(args$outdir,"per_sample"), showWarnings = F)

# Options
opts$remove_ExE_cells <- FALSE

##########################
## Load sample metadata ##
##########################

sample_metadata <- fread(args$metadata) %>%
  .[pass_rnaQC==TRUE & sample%in%args$samples & !is.na(eval(as.name(args$celltype_label)))] %>%
  setnames(args$celltype_label,"celltype")

table(sample_metadata$sample)

if (opts$remove_ExE_cells) {
  print("Removing ExE cells...")
  sample_metadata <- sample_metadata %>%
    .[!celltype%in%c("Visceral_endoderm","ExE_endoderm","ExE_ectoderm","Parietal_endoderm")]
}

#####################################
## Calculate cell type proportions ##
#####################################

celltype_proportions <- sample_metadata %>%
  .[,N:=.N,by="sample"] %>%
  .[,.(N=.N, celltype_proportion=.N/unique(N)),by=c("sample","stage","celltype")] %>%
  setorder(sample) %>% .[,sample:=factor(sample,levels=args$samples)]

#####################
## Plot per sample ##
#####################

to.plot <- celltype_proportions

# Define colours and cell type order
opts$celltype.colors <- opts$celltype.colors[names(opts$celltype.colors) %in% unique(to.plot$celltype)]
to.plot[,celltype:=factor(celltype, levels=rev(names(opts$celltype.colors)))]

samples.to.plot <- args$samples

for (i in samples.to.plot) {
  p <- ggplot(to.plot[sample==i], aes(x=celltype, y=N)) +
    geom_bar(aes(fill=celltype), stat="identity", color="black") +
    scale_fill_manual(values=opts$celltype.colors) +
    coord_flip() +
    labs(y="Number of cells") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(color="black", size=rel(0.9)),
      axis.title.x = element_text(color="black", size=rel(0.9)),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1), color="black"),
      axis.text.x = element_text(size=rel(1), color="black")
    )
  
  pdf(file.path(args$outdir,sprintf("per_sample/%s_celltype_proportions_horizontal_barplots.pdf",i)), width=7, height=5)
  print(p)
  dev.off()
}

######################
## Stacked barplots ##
######################

to.plot <- celltype_proportions

p <- ggplot(to.plot, aes(x=sample, y=celltype_proportion)) +
  geom_bar(aes(fill=celltype), stat="identity", color="black") +
  # facet_wrap(~stage, scales = "free_x", nrow=1) +
  scale_fill_manual(values=opts$celltype.colors) +
  guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color="black", size=rel(0.85)),
    # axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

pdf(sprintf("%s/celltype_proportions_stacked_barplots.pdf",args$outdir), width=7, height=5)
print(p)
dev.off()

