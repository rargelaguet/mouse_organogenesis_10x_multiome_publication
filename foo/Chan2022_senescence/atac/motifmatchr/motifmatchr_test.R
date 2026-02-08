library(chromVARmotifs)
library(BSgenome.Mmusculus.UCSC.hg38)
library(motifmatchr)

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

###############
## Load data ##
###############

peakSet.gr <- readRDS("/Users/argelagr/data/gastrulation_multiome_10x/processed/atac/archR/PeakCalls/PeakSet.rds")
length(peakSet.gr)

#####################
## Run motifmatchr ##
#####################

# Get PWM List adapted from chromVAR!
data("human_pwms_v2")
# motifs <- mouse_pwms_v2
motifs <- human_pwms_v2
obj <- .summarizeChromVARMotifs(motifs)
motifs <- obj$motifs
motifSummary <- obj$motifSummary
  
stopifnot(names(motifs)==rownames(motifSummary))

#########################
## Get motif positions ##
#########################

motifPositions <- matchMotifs(
  pwms = motifs["T_789"],
  subject = peakSet.gr,
  genome = BSgenome.Mmusculus.UCSC.hg38,
  out = "positions",
  p.cutoff = 1e-3,
  w = 7
)

tmp <- motifPositions[[1]]
tmp[seqnames(tmp)=="chr15" & start(tmp)>36706260 & end(tmp)<=36706860]

tmp[tmp$score>=10]

######################
## Get motif scores ##
######################

motifScores <- matchMotifs(
  pwms = motifs[1:10],
  subject = ArchRProject@peakSet,
  genome = BSgenome.Mmusculus.UCSC.hg38, 
  out = "scores", 
  p.cutoff = 1e-5, 
  w = 7
)

class(assay(motifScores,"motifScores"))
class(assay(motifScores,"motifMatches"))
class(assay(motifScores,"motifCounts"))

