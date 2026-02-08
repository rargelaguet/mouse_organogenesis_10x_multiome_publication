source(here::here("settings.R"))
source(here::here("utils.R"))

library(furrr)

#####################
## Define settings ##
#####################

# I/O
io$peak_annotation <- file.path(io$basedir,"processed/atac/archR/Annotations/CISBP-Scores.rds")
io$motif2gene <- file.path(io$basedir,"processed/atac/archR/Annotations/CISBP_motif2gene.txt.gz")
io$feature_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")
io$diff_results <- file.path(io$basedir,"results/atac/archR/differential/metacells/alias/PeakMatrix/parsed/diff_results.txt.gz")
io$outdir <- file.path(io$basedir,"results/atac/archR/differential/metacells/alias/PeakMatrix/parsed/pdf/motif_enrichment")
dir.create(io$outdir, showWarnings = F)

# Options
opts$groups <- opts$aliases

opts$group_variable <- "alias"

opts$conditions <- c(
  "8OFF",
  "6OFF-2ON",
  "4OFF-4ON",
  "2OFF-2ON-4OFF",
  "4ON-4OFF",
  "8ON"
)

opts$cell_lines <- c(
  "4F",
  "2F"
)

plan(multisession, workers = 4)

opts$mem_gb <- 3
options(future.globals.maxSize = opts$mem_gb * 1024 ^ 3)

###########################
## Load feature metadata ##
###########################

feature_metadata.dt <- fread(io$feature_metadata)

############################
## Load motif annotations ##
############################

peak_annotation.se <- readRDS(io$peak_annotation)

motif2gene.dt <- fread(io$motif2gene)

###############################
## Load differential results ##
###############################

diff.dt <- fread(io$diff_results) %>%
  .[,sig:=padj_fdr<0.01 & abs(logFC)>=1] %>%
  .[,sign:=ifelse(logFC>0,"Open","Closed")]

###########################################################
## Run motif enrichment for each differential comparison ##
###########################################################

motifs.to.use <- colnames(peak_annotation.se)

peak_annnotation_match.mtx <- assay(peak_annotation.se,"motifMatches")

motif_enrichment.list <- list()

# i <- opts$groups[7]; j <- opts$groups[9]; m <- motifs.to.use[1]
for (i in opts$groups) {
  for (j in opts$groups) {
    if (i!=j) {
      
      print(sprintf("%s vs %s",i,j))
      
      diff_filt.dt <- diff.dt[groupA==i & groupB==j]
      motif_enrichment.list[[sprintf("%s_vs_%s",i,j)]] <- motifs.to.use %>% future_map(function(m) { c("Open","Closed") %>% map(function(k) { 
          
        open.peaks <- diff_filt.dt[sign==k & sig==TRUE,feature]
        closed.peaks <- diff_filt.dt[sign=="Closed" & sig==TRUE,feature]
        background.peaks <- diff_filt.dt[sig==FALSE,feature]
        
        # foreground.nmatches <- assay(peak_annotation.se,"motifMatches")[,i][foreground.peaks] %>% sum
        open.nmatches <- peak_annnotation_match.mtx[,m][open.peaks] %>% sum
        background.nmatches <- peak_annnotation_match.mtx[,m][background.peaks] %>% sum
        
        p.value <- phyper(
          open.nmatches-1, 
          background.nmatches, 
          length(background.peaks)-background.nmatches,
          length(open.peaks), 
          lower.tail = F
        )
        
        data.table(motif=m, groupA=i, groupB=j, sign=k, pval=format(p.value,digits=3))
      }) %>% rbindlist }) %>% rbindlist
    }
  }
}
motif_enrichment.dt <- rbindlist(motif_enrichment.list)

fwrite(motif_enrichment.dt, file.path(io$outdir,"motif_enrichment_diff_acc.txt.gz"), sep="\t", na="NA", quote=F)
