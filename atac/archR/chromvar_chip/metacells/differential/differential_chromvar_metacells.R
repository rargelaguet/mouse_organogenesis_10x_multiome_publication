here::i_am("atac/archR/chromvar_chip/metacells/differential/differential_chromvar_metacells.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',        type="character",                               help='')
p$add_argument('--chromvar_matrix_file',        type="character",                               help='')
p$add_argument('--metadata',        type="character",                               help='')
p$add_argument('--celltypes',    type="character",    default="all", nargs="+", help='Celltypes to use')
p$add_argument('--samples',    type="character",    default="all", nargs="+", help='Samples to use')
p$add_argument('--group_variable',          type="character",   help='Group variable')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--outfile',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$motif_annotation <- "CISBP"
# args$chromvar_matrix_file <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/metacells/chromVAR_chip_%s_archr.rds",args$motif_annotation))
# args$metadata <- file.path(io$basedir,sprintf("results/atac/archR/metacells/all_cells/PeakMatrix/metacells_metadata.txt.gz"))
# args$group_variable <- "celltype"
# args$samples <- "all"
# args$celltypes <- c("NMP","Neural_crest")
# args$groupA <- "NMP"
# args$groupB <- "Neural_crest"
# args$outfile <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/metacells/differential/%s/chromVAR_%s_vs_%s.txt.gz",args$group_variable,args$groupA,args$groupB))
## END TEST ##

#####################
## Parse arguments ##
#####################

# Define cell types
if (args$celltypes[1]=="all") {
  args$celltypes <- opts$celltypes
} else {
  stopifnot(args$celltypes%in%opts$celltypes)
}

# Define samples
if (args$samples[1]=="all") {
  args$samples <- opts$samples
} else {
  stopifnot(args$samples%in%opts$samples)
}

#####################
## Define settings ##
#####################

# I/O
dir.create(dirname(args$outfile), showWarnings=F, recursive = T)

# Options
opts$groups <- c(args$groupA,args$groupB)

###################
## Load metadata ##
###################

print("Loading metacells...")

metacell_metadata.dt <- fread(args$metadata) %>%
  .[celltype%in%args$celltypes & sample%in%args$samples] %>%
  .[,celltype_genotype:=sprintf("%s_%s",celltype,genotype)]

stopifnot(args$group_variable%in%colnames(metacell_metadata.dt))

metacell_metadata.dt <- metacell_metadata.dt %>%
  .[,group:=eval(as.name(args$group_variable))] %>%
  .[group%in%opts$groups] %>%
  .[,group:=factor(group,levels=opts$groups)] %>% setorder(group) # Sort cells so that groupA comes before groupB

print(table(metacell_metadata.dt$sample))
print(table(metacell_metadata.dt$celltype))
print(table(metacell_metadata.dt$group))

######################
## Load ATAC Matrix ##
######################

print(sprintf("Fetching chromVAR matrix: '%s'...",args$chromvar_matrix_file))

# Load 
atac_chromvar_metacells.se <- readRDS(args$chromvar_matrix_file)[,metacell_metadata.dt$metacell]

# Add metadata
colData(atac_chromvar_metacells.se) <- metacell_metadata.dt %>% tibble::column_to_rownames("metacell") %>% DataFrame

# Create data.table
atac_chromvar.dt <- assay(atac_chromvar_metacells.se,"z") %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","metacell") %>%
  melt(id.vars=c("metacell"), variable.name="gene", value.name="chromvar_zscore") %>%
  merge(metacell_metadata.dt[,c("metacell","group")])

##########################
## Differential testing ##
##########################

out <- atac_chromvar.dt %>% .[,.(
    diff = mean(.SD[group==opts$groups[2],chromvar_zscore]) - mean(.SD[group==opts$groups[1],chromvar_zscore]),
    p.value = t.test(x=.SD[group==opts$groups[1],chromvar_zscore], y=.SD[group==opts$groups[2],chromvar_zscore])[["p.value"]]
  ), by="gene"] %>%
  .[,padj:=p.adjust(p.value,method="fdr")] %>% .[,p.value:=NULL] %>%
  setorder(padj, na.last=T) %>%
  .[,c("diff","padj"):=list(round(diff,2),format(padj,digits=3))]

##################
## Save results ##
##################

fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
