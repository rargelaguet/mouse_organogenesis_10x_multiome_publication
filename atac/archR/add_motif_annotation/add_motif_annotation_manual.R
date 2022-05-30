here::i_am("atac/archR/add_motif_annotation/add_motif_annotation_manual.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressPackageStartupMessages(library(motifmatchr))
suppressPackageStartupMessages(library(TFBSTools))

######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--archr_directory',     type="character", help='ArchR directory')
p$add_argument('--peak_calls',     type="character", help='Peak calls file (.rds)')
p$add_argument('--cutoff',     type="double",    default=1e-4,    help='Motif match p-value cutoff')
p$add_argument('--width',     type="integer",    default=7,    help='Minimum motif length')
p$add_argument('--motif_annotation',     type="character", default="JASPAR", help='Motif annotation (CISBP or JASPAR)')
p$add_argument('--folder_manual_motifs',     type="character", help='Folder with manually defined motifs encoded in a .jaspar format')
p$add_argument('--test',  action="store_true",  help='Test mode?')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$archr_directory <- file.path(io$basedir,"processed/atac/archR")
# args$peak_calls <- file.path(io$basedir,"processed/atac/archR/PeakCalls/PeakSet.rds")
# args$motif_annotation = "CISBP"
# args$cutoff = 0.0001
# args$width = 7
# args$test <- FALSE
# args$folder_manual_motifs <- "/Users/argelagr/data/hg38_regulation/tf_motifs/jaspar"
## END TEST ##

#####################
## Define settings ##
#####################

# I/O
io$positions_outfile <- file.path(args$archr_directory, sprintf("Annotations/%s-Positions.rds",args$motif_annotation))   
io$scores_outfile <- file.path(args$archr_directory, sprintf("Annotations/%s-Scores.rds",args$motif_annotation))
io$motif2gene_outfile <- file.path(args$archr_directory, sprintf("Annotations/%s_motif2gene.txt.gz",args$motif_annotation))
io$archR.projectMetadata <- file.path(args$archr_directory,"projectMetadata.rds")
io$peakAnnotations.output <- file.path(args$archr_directory, "Annotations/peakAnnotation.rds")

########################
## Load ArchR Project ##
########################

setwd(args$archr_directory)

addArchRGenome("mm10")
addArchRThreads(threads = 1) 

ArchRProject <- loadArchRProject(args$archr_directory)

# Load ArchR projectMetadata
ArchRProject@projectMetadata <- readRDS(io$archR.projectMetadata)

# Load peaks
ArchRProject <- addPeakSet(ArchRProject, peakSet = readRDS(args$peak_calls), force = TRUE)
names(ArchRProject@peakSet) <- sprintf("%s:%s-%s",seqnames(ArchRProject@peakSet), start(ArchRProject@peakSet), end(ArchRProject@peakSet))

# Sanity checks
ArchRProject

#########################
## Load motif database ##
#########################

print("Loading motif annotations...")

if (grepl("JASPAR",toupper(args$motif_annotation))) {
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

} else if (grepl("CISBP",toupper(args$motif_annotation))) {

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

if (args$test) {
  motifs <- motifs[1:10]
  motifSummary <- motifSummary[1:10,]
}

stopifnot(names(motifs)==rownames(motifSummary))

########################
## Add manual matches ##
########################

if (!is.null(args$folder_manual_motifs)) {
  manual_motifs <- list.files(args$folder_manual_motifs, pattern = "*.jaspar") %>% gsub(".jaspar","",.)
  print(sprintf("Loading manual motifs found in folder %s...",args$folder_manual_motifs))
  
  # j <- "MA0009.1"
  for (j in manual_motifs) {
    print(sprintf("Processing %s...",j))
    pfm <- readJASPARMatrix(file.path(args$folder_manual_motifs, sprintf("%s.jaspar",j)))
    names(pfm) <- j
    
    if (args$motif_annotation=="CISBP") {
      motifs[names(pfm)] <- toPWM(pfm, pseudocounts = 1)
      motif.metadata <- read.table(file.path(args$folder_manual_motifs, sprintf("%s.txt",j)), header = 1, colClasses = "character")[,c("ID","strand","symbol")]
    } else if (args$motif_annotation=="JASPAR") {
      motifs[names(pfm)] <- pfm
      motif.metadata <- read.table(file.path(args$folder_manual_motifs, sprintf("%s.txt",j)), header = 1, colClasses = "character")
    }
    
    stopifnot(colnames(motif.metadata)==colnames(motifSummary))
    motifSummary <- rbind(motifSummary, motif.metadata)
  }
}

# Save motif2gene
motif2gene.dt <- motifSummary %>%
  as.data.table(keep.rownames = T) %>% setnames("rn","motif") %>% .[,strand:=NULL] %>% setnames("symbol","gene")
fwrite(motif2gene.dt, io$motif2gene_outfile, sep="\t", quote=F)

###################################
## Run motifmatchr: get matches ##
###################################

# NOT NECESSARY: SCORES ALSO OUTPUTS MATCHES

# motifmatcher_matches.out <- matchMotifs(
#   pwms = motifs,
#   subject = ArchRProject@peakSet,
#   genome = BSgenome.Mmusculus.UCSC.mm10,
#   out = "matches",
#   p.cutoff = args$cutoff,
#   w = args$width
# )
# 
# io$outfile <- file.path(args$archr_directory, sprintf("Annotations/%s_cutoff%s_width%s-Matches.rds",args$motif_annotation,args$cutoff,args$width))   
# saveRDS(motifmatcher_matches.out, io$outfile)

###################################
## Run motifmatchr: get positions ##
###################################

print("Running motifmatchr to get positions...")

motifmatcher_positions.out <- matchMotifs(
  pwms = motifs,
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  out = "positions",
  p.cutoff = args$cutoff,
  w = args$width
)

# Save
saveRDS(motifmatcher_positions.out, io$positions_outfile)

#################################
## Run motifmatchr: get scores ##
#################################

print("Running motifmatchr to get scores...")

motifmatcher_scores.out <- matchMotifs(
  pwms = motifs,
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  out = "scores",
  p.cutoff = args$cutoff,
  w = args$width
)

# unique(assay(motifmatcher_scores.out[,rownames(motifSummary[motifSummary$symbol=="GATA1",])], "motifScores")[,1])
# unique(assay(motifmatcher_scores.out[,rownames(motifSummary[motifSummary$symbol=="NKX2-5",])], "motifScores")[,1])

# print
print("MotifScores:"); assay(motifmatcher_scores.out,"motifScores")[1:5,1:5]
print("MotifMatches:"); assay(motifmatcher_scores.out,"motifMatches")[1:5,1:5]
print("MotifCounts:"); assay(motifmatcher_scores.out,"motifCounts")[1:5,1:5]
print(sprintf("Matches with cutoff=%s and width=%s: %.2f%%", args$cutoff, args$width, 100*mean(assay(motifmatcher_scores.out,"motifMatches")>0)))

# Save
saveRDS(motifmatcher_scores.out, io$scores_outfile)


###########################
## Save motif annotation ##
###########################

print("Saving peakAnnotations file...")

if (file.exists(io$peakAnnotations.output)) {
  peakAnnotation <- readRDS(io$peakAnnotations.output)
} else {
  # peakAnnotation <- list()
  peakAnnotation <- S4Vectors::SimpleList()

}

if (args$motif_annotation%in%names(peakAnnotation)) {
  print(sprintf("%s already exists in peakAnnotation.rds, replacing...",args$motif_annotation))
}

peakAnnotation[[args$motif_annotation]]$Name <- args$motif_annotation
peakAnnotation[[args$motif_annotation]]$motifs <- motifs
peakAnnotation[[args$motif_annotation]]$motifSummary <- motifSummary
peakAnnotation[[args$motif_annotation]]$Positions <- io$positions_outfile
peakAnnotation[[args$motif_annotation]]$Matches <- io$scores_outfile

saveRDS(peakAnnotation, io$peakAnnotations.output)

##########
## TEST ##
##########

# assay(motifmatcher_scores.out["chr6:97864440-97865040","MA0009.1"])
# mean(assay(motifmatcher_scores.out[,"MA0009.1"],"motifCounts")[,1])

# motifmatcher_scores.out <- readRDS(file.path(args$archr_directory,"Annotations/CISBP-Scores.rds"))
# motifmatcher_positions.out <- readRDS(file.path(args$archr_directory,"Annotations/CISBP-Positions.rds"))
# stopifnot(colnames(motifmatcher_scores.out)==names(motifmatcher_positions.out))

# i <- "chr1:3035600-3036200"
# peak.gr <- rowRanges(motifmatcher_scores.out[i,])
# TFs <- which(assay(motifmatcher_scores.out[i,],"motifMatches")[1,]) %>% names
# 
# j <- "TFAP2D_2"
# assay(motifmatcher_scores.out[i,j],"motifCounts")
# assay(motifmatcher_scores.out[i,j],"motifScores")
# 
# tmp <- motifmatcher_positions.out[[j]][seqnames(motifmatcher_positions.out[[j]])=="chr1"]
# bar <- tmp[start(tmp)>=start(peak.gr) & end(tmp)<=end(peak.gr)]
# 
# findOverlaps(bar,tmp)
# tmp <- TFs %>% map(function(j) {
#   print(j)
#   hits <- findOverlaps(
#     query = peak.gr,
#     subject = motifmatcher_positions.out[[j]],
#     ignore.strand=TRUE
#   ) %>% subjectHits()
#   if (length(hits)>0) {
#     tmp <- motifmatcher_positions.out[[j]][hits]
#     tmp$TF <- j
#     return(tmp)
#   } else {
#     print(sprintf("No hits found for %s",j))
#     return(NULL)
#   }
# })
# motif_locations.dt <- tmp[!sapply(tmp,is.null)] %>% as(., "GRangesList") %>% unlist %>%
#   as.data.table %>% setnames("score","motif_score") %>% setorder(start)