source(here::here("settings.R"))
source(here::here("utils.R"))

library(furrr)

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# io$peak_annotation <- file.path(io$basedir,sprintf("processed/atac/bigwig/deeptools/ATAC_peaks/motif_annotation/%s/%s-Scores.rds",opts$motif_annotation,opts$motif_annotation))
io$peak_annotation <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s-Scores.rds",opts$motif_annotation))
io$motif2gene <- file.path(io$basedir,sprintf("processed/atac/bigwig/deeptools/ATAC_peaks/motif_annotation/%s/%s_motif2gene.txt.gz",opts$motif_annotation,opts$motif_annotation))
io$feature_metadata <- file.path(io$basedir,"processed/atac/archR/peak_metadata.tsv.gz")
io$limma_results <- file.path(io$basedir,"results/atac/archR/limma/limma_results.rds"); 
io$outdir <- file.path(io$basedir,sprintf("results/atac/archR/motif_enrichment/%s",opts$motif_annotation)); dir.create(io$outdir, showWarnings = F)

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

limma_output.list <- readRDS(io$limma_results)

############################################
## Run motif enrichment for each contrast ##
############################################

contrasts <- names(limma_output.list)

motifs.to.use <- colnames(peak_annotation.se)

peak_annnotation_match.mtx <- assay(peak_annotation.se,"motifMatches")

# i <- contrasts[1]; m <- motifs.to.use[1]; k <- "Open"
motif_enrichment.list <- contrasts %>% future_map(function(i) {
  
  limma_output.dt <- limma_output.list[[i]] %>%
    .[,sig:=pval<0.01] %>%
    .[,sign:=ifelse(logFC>0,"Open","Closed")]
  
  background.peaks <- limma_output.dt[pval>0.30,feature]
  
  motifs.to.use %>% future_map(function(m) { c("Open","Closed") %>% map(function(k) { 
    
    foreground.peaks <- limma_output.dt[sign==k & sig==TRUE,feature]
    foreground.nmatches <- peak_annnotation_match.mtx[,m][foreground.peaks] %>% sum
    background.nmatches <- peak_annnotation_match.mtx[,m][background.peaks] %>% sum
    
    p.value <- phyper(
      foreground.nmatches-1, 
      background.nmatches, 
      length(background.peaks)-background.nmatches,
      length(foreground.peaks), 
      lower.tail = F
    )
    data.table(
      contrast = i, motif = m, sign = k, pval = p.value, 
      foreground_peaks = length(foreground.peaks), background_peaks = length(background.peaks), 
      foreground_nmatches = foreground.nmatches, background_nmatches = background.nmatches
    )
  }) %>% rbindlist }) %>% rbindlist
})

motif_enrichment.dt <- rbindlist(motif_enrichment.list)

fwrite(motif_enrichment.dt, file.path(io$outdir,"motif_enrichment_limma_acc.txt.gz"), sep="\t", na="NA", quote=F)

