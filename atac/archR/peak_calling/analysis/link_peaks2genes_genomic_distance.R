here::i_am("atac/archR/peak_calling/analysis/link_peaks2genes_genomic_distance.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--gene_metadata',  type="character",                help='Gene metadata') 
p$add_argument('--peak_metadata',  type="character",                help='Peak metadata') 
p$add_argument('--gene_window',  type="integer",  default=1e5,               help='Genomic window size') 
p$add_argument('--outdir',  type="character",                help='Output directory') 
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$gene_metadata <- io$gene_metadata
# args$peak_metadata <- file.path(io$basedir,"processed_new/atac/archR/PeakCalls/peak_metadata.tsv.gz")
# args$gene_window <- 5e5  # maximum window length for the overlap
# args$outdir <- file.path(io$basedir,"results_new/atac/archR/peak_calling/peaks2genes")
## END TEST ##

# I/O
dir.create(args$outdir, showWarnings=F)

###############
## Load data ##
###############

# Load gene metadata
gene_metadata <- fread(args$gene_metadata) %>% 
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  setnames("symbol","gene") %>%
  .[, c("chr","start","end","gene","ens_id","strand")]

# Load peak metadata
peakSet.dt <- fread(args$peak_metadata) %>%
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  .[,c("chr","start","end")] %>%
  .[,peak:=sprintf("chr%s:%s-%s",chr,start,end)] %>%
  setkey(chr,start,end)

#############
## Overlap ##
#############

gene_metadata.ov <- copy(gene_metadata) %>%
  .[strand=="+",c("gene.start","gene.end"):=list(start,end)] %>%
  .[strand=="-",c("gene.start","gene.end"):=list(end,start)] %>%
  .[strand=="+",c("start","end"):=list (gene.start-args$gene_window, gene.end+args$gene_window)] %>%
  .[strand=="-",c("end","start"):=list (gene.start+args$gene_window, gene.end-args$gene_window)] %>% 
  # .[,strand:=NULL] %>% 
  setkey(chr,start,end)

stopifnot((gene_metadata.ov$end-gene_metadata.ov$start)>0)

ov <- foverlaps(
  peakSet.dt,
  gene_metadata.ov,
  nomatch = NA
) %>%  .[,c("start","end"):=NULL] %>%
  setnames(c("i.start","i.end"),c("peak.start","peak.end")) %>%
  .[,peak.mean:=(peak.start+peak.end)/2] %>%
  # calculate distance from the peak to the genebody
  .[,dist:=min(abs(gene.end-peak.mean), abs(gene.start-peak.mean)), by=c("gene","ens_id","peak","strand")] %>%
  .[strand=="+" & peak.mean>gene.start & peak.mean<gene.end,dist:=0] %>%
  .[strand=="-" & peak.mean<gene.start & peak.mean>gene.end,dist:=0]


# ov[peak=="18_64485555_64486155"]
# gene_metadata[gene=="Fech"]
# gene_metadata.ov[gene=="Fech"]
# ov[gene=="Fech" & peak=="chr18:64485555-64486155"]
# ov[peak=="chr7:103850833-103851433"]

# Select nearest gene
ov_nearest <- ov %>%
  .[.[,.I[dist==min(dist)], by=c("peak")]$V1] %>%
  .[complete.cases(.)] %>%
  .[!duplicated(peak)]

# Sanity check  
# ov_nearest$gene[(duplicated(ov_nearest$peak))]

##########
## Save ##
##########

fwrite(ov, file.path(args$outdir,"peaks2genes_all.txt.gz"), sep="\t", na="NA")
fwrite(ov_nearest, file.path(args$outdir,"peaks2genes_nearest.txt.gz"), sep="\t", na="NA")
