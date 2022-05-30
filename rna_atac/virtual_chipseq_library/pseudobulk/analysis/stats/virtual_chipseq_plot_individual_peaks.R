library(GenomicRanges)

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

# I/O
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/individual_peaks",opts$motif_annotation)); dir.create(io$outdir, showWarnings = F)

io$motifmatcher.se <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s-Scores.rds",opts$motif_annotation))
io$motifmatcher_positions.se <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s-Positions.rds",opts$motif_annotation))
io$motif2gene <- file.path(io$basedir,sprintf("processed/atac/archR/Annotations/%s_motif2gene.txt.gz",opts$motif_annotation))
io$motifmatcher_virtual_chip.se <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/motifmatchr_virtual_chip.rds",opts$motif_annotation))
io$peak_metadata <- file.path(io$basedir,"processed/atac/archR/PeakCalls/peak_metadata.tsv.gz")

######################
## Load motifmatchR ##
######################

motifmatcher.se <- readRDS(io$motifmatcher_virtual_chip.se)
motifmatcher_positions.se <- readRDS(io$motifmatcher_positions.se)

# stopifnot(colnames(motifmatcher.se)==names(motifmatcher_positions.se))

################
## Subset TFs ##
################

motif2gene.dt <- fread(io$motif2gene)

motifs <- intersect(names(motifmatcher_positions.se),motif2gene.dt$motif)
genes <- intersect(colnames(virtual_chip.mtx),motif2gene.dt$gene)

motif2gene_filt.dt <- motif2gene.dt[motif%in%motifs & gene%in%genes]
motifs <- motif2gene_filt.dt$motif
genes <- motif2gene_filt.dt$gene

tmp <- genes; names(tmp) <- motifs

stopifnot(motif2gene_filt.dt$motif%in%names(motifmatcher_positions.se))
stopifnot(motif2gene_filt.dt$gene%in%colnames(virtual_chip.mtx))

motifmatcher.se <- motifmatcher.se[,genes]
motifmatcher_positions.se <- motifmatcher_positions.se[motifs]
names(motifmatcher_positions.se) <- tmp[names(motifmatcher_positions.se)]

########################
## Load peak metadata ##
########################

peak_metadata.dt <- fread(io$peak_metadata) %>%
  .[,c("chr","start","end","score")] %>%
  .[,idx:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[idx%in%rownames(motifmatcher.se)]

peak_metadata.gr <- makeGRangesFromDataFrame(peak_metadata.dt, keep.extra.columns = F)

###########################
## Plot individual peaks ##
###########################

# opts$min_motif_score <- 4

peaks.to.plot <- c("chr6:87617887-87618487")
peaks.to.plot <- fread("/Users/argelagr/data/gastrulation_multiome_10x/test/results/rna_atac/virtual_chipseq/pseudobulk/CISBP/TAL1.bed.gz") %>% .[,idx:=sprintf("%s:%s-%s",V1,V2,V3)] %>% .$idx %>% head(n=6)

stopifnot(peaks.to.plot%in%rownames(motifmatcher.se))

# i <- "chr1:133329196-133329796"
for (i in peaks.to.plot) {
  
  # Identify motif locations + scores
  peak.gr <- rowRanges(motifmatcher.se[i,])
  TFs <- which(assay(motifmatcher.se[i,],"motifMatches")[1,]) %>% names
  tmp <- TFs %>% map(function(j) {
    hits <- findOverlaps(
      query = peak.gr,
      subject = motifmatcher_positions.se[[j]],
      ignore.strand=TRUE
    ) %>% subjectHits()
    if (length(hits)>0) {
      tmp <- motifmatcher_positions.se[[j]][hits]
      tmp$TF <- j
      return(tmp)
    } else {
      print(sprintf("No hits found for %s",j))
      return(NULL)
    }
  }) 
  motif_locations.dt <- tmp[!sapply(tmp,is.null)] %>% as(., "GRangesList") %>% unlist %>%
    as.data.table %>% setnames("score","motif_score") %>% setorder(start) 
  
  # Add in silico ChIP-seq scores
  stopifnot(motif_locations.dt$TF%in%colnames(motifmatcher.se))
  insilico_chip_score <- assay(motifmatcher.se,"VirtualChipScores")[i,motif_locations.dt$TF]
  
  # Prepare data.table for plotting
  to.plot <- motif_locations.dt %>% .[,insilico_chip_score:=insilico_chip_score]
  
  # Remove palindromic motif duplicates
  to.plot[,N:=.N,by=c("seqnames","start","end","TF")] 
  to.plot <- rbind(to.plot[N==1],to.plot[N==2 & strand=="+"]) %>% .[,N:=NULL]
  
  # ignore negative binding scores
  to.plot <- to.plot[insilico_chip_score<=0,insilico_chip_score:=0]
  # Filter by mininum motif score
  # to.plot <- to.plot %>% .[score>=opts$min_motif_score]
  # to.plot <- motif_locations.dt %>% .[score<=opts$min_motif_score,score:=opts$min_motif_score]
  
  # Plot
  to.plot.text <- to.plot %>% .[insilico_chip_score>=0.25] %>% setorder(-insilico_chip_score)
  
  p <- ggplot(to.plot, aes_string(x="start", y="insilico_chip_score")) +
    # geom_jitter(aes(shape=strand, alpha=insilico_chip_score, size=insilico_chip_score), fill="gray70", width = 10, height=0.01) +
    ggrastr::geom_jitter_rast(aes(shape=strand, alpha=insilico_chip_score, size=insilico_chip_score), fill="gray70", width = 10, height=0.01) +
    scale_x_continuous(limits=c(start(peak.gr),end(peak.gr))) +
    # scale_fill_gradientn(colours = rev(terrain.colors(10))) +
    # scale_fill_gradient2(low = "gray50", high = "gray10") +
    # scale_fill_gradient2(low = "blue", mid="gray90", high = "red", limits=c(0,1)) +
    scale_shape_manual(values=c(21,24)) +
    scale_size_continuous(range=c(1,3.5)) +
    scale_alpha_continuous(range=c(0.25,1)) +
    ggrepel::geom_text_repel(aes(label=TF), size=3, max.overlaps=Inf, data=to.plot.text) +
    geom_hline(yintercept=0.25, linetype="dashed") +
    # geom_segment(aes_string(xend="score"), size=0.5, yend=0) +
    labs(x="Genomic position", y="in silico TF binding score", title=i) +
    theme_classic() +
    guides(size="none", alpha="none") +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position = "top",
      axis.text = element_text(size=rel(0.8), color="black"),
      axis.title = element_text(size=rel(0.9), color="black")
    )
  
  pdf(file.path(io$outdir,sprintf("peak_%s_motif_location_virtual_chip.pdf",sub(":","-",i))), width=6, height=4)
  print(p)
  dev.off()
}
