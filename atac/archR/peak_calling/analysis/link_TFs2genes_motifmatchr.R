
stop("TO FINISH")

#####################
## Define settings ##
#####################

if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/gastrulation_multiome_10x/settings.R")
  source("/Users/ricard/gastrulation_multiome_10x/utils.R")
} else if (grepl("ebi",Sys.info()['nodename'])) {
  source("/homes/ricard/gastrulation_multiome_10x/settings.R")
  source("/homes/ricard/gastrulation_multiome_10x/utils.R")
} else {
  stop("Computer not recognised")
}

# I/O
io$outdir <- paste0(io$basedir,"/results/atac/archR/peak_calling/TF2genes")

# Options
opts <- list()
opts$gene_window <- 1e5 # window length for the overlap

###############
## Load data ##
###############

# Load gene metadata
gene_metadata <- fread(io$gene_metadata) %>% 
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  setnames("symbol","gene") %>%
  .[, c("chr","start","end","gene")] %>%
  setkey(chr,start,end)

# Load peak metadata
peakSet.dt <- fread(io$archR.peak.metadata) %>%
  .[,chr:=as.factor(sub("chr","",chr))] %>%
  .[,c("chr","start","end")] %>%
  .[,peak:=sprintf("%s_%s_%s",chr,start,end)] %>%
  setkey(chr,start,end)

##########################
## Load peak2gene links ##
##########################

#################################
## Load motifmatcher Matches ##
#################################

motifmatcher.se <- readRDS(sprintf("%s/Annotations/Motif_cisbp-Matches-In-Peaks.rds",io$archR.directory))

# Rename TFs
colnames(motifmatcher.se) <- colnames(motifmatcher.se) %>% toupper %>% stringr::str_split(.,"_") %>% map_chr(1)
motifmatcher.se <- motifmatcher.se[,!duplicated(colnames(motifmatcher.se))]

# Rename peaks
tmp <- rowRanges(motifmatcher.se)
rownames(motifmatcher.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))

# Subset pekas
# motifmatcher.se <- motifmatcher.se[unique(cor_dt$peak),]

#############
## Overlap ##
#############

ov <- foverlaps(
  peakSet.dt,
  gene_metadata[, c("chr","start","end","gene")],
  nomatch = NA
) %>% 
  setnames(c("i.start","i.end"),c("peak.start","peak.end")) %>%
  setnames(c("start","end"),c("gene.start","gene.end")) %>%
  .[,c("gene.start","gene.end") := list (gene.start+opts$gene_window, gene.end-opts$gene_window)] %>%
  # .[,c("start_dist","end_dist"):=list( abs(gene.end-peak.start), abs(gene.start-peak.end))] %>%
  .[,c("start_dist","end_dist"):=list( gene.end-peak.start, gene.start-peak.end)] %>%
  .[,c("start_dist","end_dist"):=list( ifelse(end_dist<0 & start_dist>0,0,start_dist), ifelse(end_dist<0 & start_dist>0,0,end_dist) )] %>%
  .[,dist:=ifelse(abs(start_dist)<abs(end_dist),abs(start_dist),abs(end_dist))] %>% .[,c("start_dist","end_dist"):=NULL]

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

fwrite(ov, paste0(io$outdir,"/peaks2genes_all.txt.gz"), sep="\t", na="NA")
fwrite(ov_nearest, paste0(io$outdir,"/peaks2genes_nearest.txt.gz"), sep="\t", na="NA")
