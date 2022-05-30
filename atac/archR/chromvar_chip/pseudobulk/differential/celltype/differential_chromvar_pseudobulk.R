here::i_am("atac/archR/chromvar_chip/pseudobulk/differential/celltype/differential_chromvar_pseudobulk.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',        type="character",                               help='')
p$add_argument('--chromvar_chip',        type="character",                               help='')
p$add_argument('--groupA',    type="character",    help='group A')
p$add_argument('--groupB',    type="character",    help='group B')
p$add_argument('--outfile',          type="character",                               help='Output directory')

args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# io$basedir <- file.path(io$basedir,"test")
# args <- list()
# args$motif_annotation <- "JASPAR"
# args$chromvar_chip <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk_with_replicates/chromVAR_chip_%s_archr.rds",args$motif_annotation))
# args$groupA <- "ExE_ectoderm"
# args$groupB <- "Caudal_neurectoderm"
# args$outfile <- file.path(io$basedir,sprintf("results/atac/archR/chromvar_chip/pseudobulk_with_replicates/differential/celltype/chromVAR_%s_vs_%s.txt.gz",args$groupA,args$groupB))
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
dir.create(dirname(args$outfile), showWarnings=F, recursive = T)

# Options
opts$groups <- c(args$groupA,args$groupB)

# stupid stuff but otherwise the snakemake pipeline doesn't work
if (args$groupA==args$groupB) {
  out <- data.table(feature=NA, diff=NA, padj=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  warning("groupA and groupB are the same, saving an empty file...")
  quit(status=0)
}

###########################
## Load chromVAR results ##
###########################

print(sprintf("Fetching chromVAR results: '%s'...",args$chromvar_chip))

# Load 
atac_chromvar.se <- readRDS(args$chromvar_chip)

# parse
if (!"celltype"%in%colnames(colData(atac_chromvar.se))) {
  atac_chromvar.se$celltype <- colnames(atac_chromvar.se) %>% strsplit("_rep") %>% map_chr(1)
}

atac_chromvar.se <- atac_chromvar.se[,atac_chromvar.se$celltype%in%opts$groups]
atac_chromvar.se$celltype <- factor(atac_chromvar.se$celltype, levels=opts$groups)

# check that we have pseudobulk replicates for both cell types
if (any(!opts$groups%in%unique(atac_chromvar.se$celltype))) {
  warning("groups not found, saving an empty file...")
  out <- data.table(feature=NA, diff=NA, padj=NA)
  fwrite(out, args$outfile, sep="\t", na="NA", quote=F)
  quit(status=0)
}

# Create data.table
atac_chromvar.dt <- assay(atac_chromvar.se,"z") %>% t %>%
  as.data.table(keep.rownames = T) %>%
  setnames("rn","sample") %>%
  melt(id.vars=c("sample"), variable.name="gene", value.name="chromvar_zscore")
  
tmp <- data.table(sample=colnames(atac_chromvar.se), group=atac_chromvar.se$celltype)
atac_chromvar.dt <- atac_chromvar.dt %>% merge(tmp[,c("sample","group")])

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
