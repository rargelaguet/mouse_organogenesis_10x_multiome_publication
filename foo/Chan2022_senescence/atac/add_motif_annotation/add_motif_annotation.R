source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("atac/add_motif_annotation/utils.R"))


suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(motifmatchr))
suppressPackageStartupMessages(library(TFBSTools))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation = "CISBP"
# opts$motif_annotation = "JASPAR"
opts$cutoff = 0.0001
opts$width = 7

# I/O
# io$peak_calls <- file.path(io$basedir,"processed/atac/archR/PeakSet.rds")
io$folder_manual_motifs <- "/share/rargelaguet/data/hg38_regulation/tf_motifs/cisbp"
# io$folder_manual_motifs <- "/share/rargelaguet/data/hg38_regulation/tf_motifs/jaspar"
io$outdir <- file.path(io$basedir,sprintf("processed/atac/motif_annotation/%s",opts$motif_annotation)); dir.create(io$outdir, showWarnings=F, recursive=T)
io$positions_outfile <- file.path(io$outdir, sprintf("%s-Positions.rds",opts$motif_annotation))   
io$scores_outfile <- file.path(io$outdir, sprintf("%s-Scores.rds",opts$motif_annotation))
io$motif2gene_outfile <- file.path(io$outdir, sprintf("%s_motif2gene.txt.gz",opts$motif_annotation))
io$peakAnnotations.output <- file.path(io$outdir, "peakAnnotation.rds")

###################
## Load peak set ##
###################

feature_metadata.dt <- fread("/Users/rargelaguet/data/Chan2022_senescence/processed/atac/jack/peak_calls.bed") %>%
  setnames(c("chr","start","end")) %>%
  .[chr%in%paste0("chr",c(1:22,"X","Y"))]
fwrite(feature_metadata.dt, "/Users/rargelaguet/data/Chan2022_senescence/processed/atac/ATAC_peaks.bed", sep="\t", quote=F, na="NA", col.names = F)

peakSet.gr <- makeGRangesFromDataFrame(feature_metadata.dt, keep.extra.columns = F)

# peakSet.gr <- readRDS(io$peak_calls)

#########################
## Load motif database ##
#########################

print("Loading motif annotations...")

if (grepl("JASPAR",toupper(opts$motif_annotation))) {
  library(JASPAR2020)

  # this returns a PFMatrixList
  # Homo sapiens: 633 TF motifs (MA0030.1, MA0031.1, etc.)
  # Mus musculus: 107 TF motifs
  # motifs <- getMatrixSet(JASPAR2020, list(species = "Mus musculus", collection = "CORE"))
  motifs <- getMatrixSet(JASPAR2020, list(species = "Homo sapiens", collection = "CORE"))

  # this returns a list with two entries: (1) motifs (PFMatrixList), (2) motifSummary (a data.frame). 
  # Note TFs are renamed to capitalised gene symbol (FOXF2_1, FOXD1_2, IRF2_3)
  obj <- .summarizeJASPARMotifs(motifs)
  motifs <- obj$motifs
  motifSummary <- obj$motifSummary
  motifSummary$name <- NULL

  # Filter out TF motifs
  tmp <- rownames(obj$motifSummary)
  tmp <- tmp[grep("\\.\\.",tmp,invert=TRUE)] # fusion proteins
  motifs <- motifs[tmp]
  motifSummary <- motifSummary[tmp,]

  # Rename TF motifs
  motifSummary$symbol <- strsplit(rownames(motifSummary),"\\_") %>% map_chr(1)
  rename <- c("(NKX[0-9])\\.([0-9])" = "\\1-\\2"); motifSummary$symbol <- stringr::str_replace_all(motifSummary$symbol,rename)
  motifSummary$symbol <- strsplit(motifSummary$symbol,"\\.") %>% map_chr(1)
  
  stopifnot(sum(duplicated(motifSummary$name))==0)

} else if (grepl("CISBP",toupper(opts$motif_annotation))) {

  library(chromVARmotifs)

  data("human_pwms_v2")
  motifs <- human_pwms_v2 # mouse_pwms_v2
  obj <- .summarizeChromVARMotifs(motifs)
  motifs <- obj$motifs
  motifSummary <- obj$motifSummary

  # Filter out TF motifs

  # Rename TF motifs
  motifSummary$symbol <- motifSummary$name; motifSummary$name <- NULL
  motifSummary$symbol <- toupper(motifSummary$symbol)
  rename <- c(
    TCFE = "TFE", NKX1 = "NKX1-", NKX2 = "NKX2-", NKX3 = "NKX3-", 
    NKX4 = "NKX4-", NKX5 = "NKX5-", NKX6 = "NKX6-", FOXF1A = "FOXF1", 
    HMGA1RS1 = "RS1", `MYCL1$` = "MYCL", `DUX$` = "DUXF3", `DUXBL$` = "DUXBL1", 
    `PIT1$` = "PROP1", ENSMUSG00000079994 = "SOX1", TCFAP = "TFAP"
  )
  motifSummary$symbol <- stringr::str_replace_all(motifSummary$symbol,rename)
  
} else {
  stop("Motif annotation not recognised")
}

stopifnot(names(motifs)==rownames(motifSummary))

########################
## Add manual matches ##
########################

if (!is.null(io$folder_manual_motifs)) {
  manual_motifs <- list.files(io$folder_manual_motifs, pattern = "*.jaspar") %>% gsub(".jaspar","",.)
  print(sprintf("Loading manual motifs found in folder %s...",io$folder_manual_motifs))
  
  # j <- "MA0009.1"
  for (j in manual_motifs) {
    print(sprintf("Processing %s...",j))
    pfm <- readJASPARMatrix(file.path(io$folder_manual_motifs, sprintf("%s.jaspar",j)))
    names(pfm) <- j
    
    if (opts$motif_annotation=="CISBP") {
      motifs[names(pfm)] <- toPWM(pfm, pseudocounts = 1)
      motif.metadata <- read.table(file.path(io$folder_manual_motifs, sprintf("%s.txt",j)), header = 1, colClasses = "character")[,c("ID","strand","symbol")]
    } else if (opts$motif_annotation=="JASPAR") {
      motifs[names(pfm)] <- pfm
      motif.metadata <- read.table(file.path(io$folder_manual_motifs, sprintf("%s.txt",j)), header = 1, colClasses = "character")
    }
    
    stopifnot(colnames(motif.metadata)==colnames(motifSummary))
    motifSummary <- rbind(motifSummary, motif.metadata)
  }
}

# Save motif2gene
motif2gene.dt <- motifSummary %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>% .[,strand:=NULL] %>% setnames("symbol","gene")

fwrite(motif2gene.dt, io$motif2gene_outfile, sep="\t", quote=F)

####################################
## Run motifmatchr: get positions ##
####################################

print("Running motifmatchr to get positions...")

motifmatcher_positions.out <- matchMotifs(
  pwms = motifs,
  subject = peakSet.gr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "positions",
  p.cutoff = opts$cutoff,
  w = opts$width
)

# Save
saveRDS(motifmatcher_positions.out, io$positions_outfile)

#################################
## Run motifmatchr: get scores ##
#################################

print("Running motifmatchr to get scores...")

motifmatcher_scores.out <- matchMotifs(
  pwms = motifs,
  subject = peakSet.gr,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  out = "scores",
  p.cutoff = opts$cutoff,
  w = opts$width
)
rownames(motifmatcher_scores.out) <- sprintf("%s:%s-%s",as.character(seqnames(peakSet.gr)), start(peakSet.gr), end(peakSet.gr))

# print
print("MotifScores:"); assay(motifmatcher_scores.out,"motifScores")[1:5,1:5]
print("MotifMatches:"); assay(motifmatcher_scores.out,"motifMatches")[1:5,1:5]
print("MotifCounts:"); assay(motifmatcher_scores.out,"motifCounts")[1:5,1:5]
print(sprintf("Matches with cutoff=%s and width=%s: %.2f%%", opts$cutoff, opts$width, 100*mean(assay(motifmatcher_scores.out,"motifMatches")>0)))

# Save
saveRDS(motifmatcher_scores.out, io$scores_outfile)

###########################
## Save motif annotation ##
###########################

# print("Saving peakAnnotations file...")
# 
# if (file.exists(io$peakAnnotations.output)) {
#   peakAnnotation <- readRDS(io$peakAnnotations.output)
# } else {
#   # peakAnnotation <- list()
#   peakAnnotation <- S4Vectors::SimpleList()
# 
# }
# 
# if (opts$motif_annotation%in%names(peakAnnotation)) {
#   print(sprintf("%s already exists in peakAnnotation.rds, replacing...",opts$motif_annotation))
# }
# 
# peakAnnotation[[opts$motif_annotation]]$Name <- opts$motif_annotation
# peakAnnotation[[opts$motif_annotation]]$motifs <- motifs
# peakAnnotation[[opts$motif_annotation]]$motifSummary <- motifSummary
# peakAnnotation[[opts$motif_annotation]]$Positions <- io$positions_outfile
# peakAnnotation[[opts$motif_annotation]]$Matches <- io$scores_outfile
# 
# saveRDS(peakAnnotation, io$peakAnnotations.output)
