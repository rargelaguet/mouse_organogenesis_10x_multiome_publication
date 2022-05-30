suppressPackageStartupMessages(library(argparse))

here::i_am("atac/archR/chromvar_chip/pseudobulk/differential/plot_differential_chromvar_chip_pseudobulk.R")

################################
## Initialize argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--motif_annotation',  type="character",              help='Motif annotation') 
p$add_argument('--chromvar_diff_pseudobulk_dir',  type="character",              help='Motif annotation') 
p$add_argument('--outdir',  type="character",              help='Motif annotation') 
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

# load default setings
source(here::here("settings.R"))
source(here::here("utils.R"))

## START TEST ##
args$motif_annotation <- "CISBP"
args$chromvar_diff_pseudobulk_dir <- file.path(io$basedir,sprintf("results_new/atac/archR/chromvar_chip/pseudobulk/differential/%s",args$motif_annotation))
args$outdir <- file.path(io$basedir,sprintf("results_new/atac/archR/chromvar_chip/pseudobulk/differential/%s/pdf",args$motif_annotation))
## END TEST ##

dir.create(args$outdir, showWarnings = F)

########################################################
## Load precomputed differential chromVAR-ChIP scores ##
########################################################

diff.dt <- (1:length(opts$celltypes)) %>% map(function(i) {
  (i:length(opts$celltypes)) %>% map(function(j) {
    if (i!=j) {
      file <- file.path(args$chromvar_diff_pseudobulk_dir,sprintf("%s_vs_%s_chromVAR_chip_pseudobulk.txt.gz",opts$celltypes[[i]],opts$celltypes[[j]]))
      if (file.exists(file)) {
        fread(file) %>% 
          .[,groupA:=factor(opts$celltypes[[i]],levels=opts$celltypes)] %>% 
          .[,groupB:=factor(opts$celltypes[[j]],levels=opts$celltypes)] %>%
          return
      }
    }
  }) %>% rbindlist
}) %>% rbindlist

##########
## Plot ##
##########

# celltypes.to.plot <- c("Gut","Erythroid3")
# genes.to.plot <- c("TAL1")

opts$xlim.max <- 3
opts$xlim.min <- -3

# i <- "Gut"; j <- "Erythroid3"
for (i in opts$celltypes) {
  for (j in opts$celltypes) {
    
    to.plot <- diff.dt %>%
      .[groupA==i & groupB==j] %>% 
      .[,gene:=factor(gene,levels=rev(gene))] %>%
      .[diff>=opts$xlim.max,diff:=opts$xlim.max] %>%
      .[diff<=opts$xlim.min,diff:=opts$xlim.min]
    
    p <- ggplot(to.plot, aes(x=diff, y=gene)) +
      geom_jitter(aes(color=abs(diff), alpha=abs(diff)), width = 0.15) +
      ggrepel::geom_text_repel(data=head(to.plot[diff>0],n=10), aes(x=diff, y=gene, label=gene), size=4, max.overlaps=Inf) +
      ggrepel::geom_text_repel(data=head(to.plot[diff<0],n=10), aes(x=diff, y=gene, label=gene), size=4, max.overlaps=Inf) +
      scale_color_gradient(low = "gray80", high = "red") +
      scale_alpha_continuous(range=c(0.25,1)) +
      coord_cartesian(xlim=c(opts$xlim.min,opts$xlim.max)) +
      theme_classic() +
      labs(x="Differential motif accessibility (chromVAR)", y="") +
      # coord_flip() +
      annotate("text", x=opts$xlim.min/1.5, y=75, size=4, label=sprintf("(+) %s",i)) +
      annotate("text", x=opts$xlim.max/1.5, y=75, size=4, label=sprintf("(+) %s",j)) +
      geom_segment(x=0, xend=0, y=0, yend=nrow(to.plot), color="black", size=0.25, linetype="dashed") +
      theme(
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=rel(1.0), color="black")
      )
    
    
    pdf(file.path(args$outdir,sprintf("%s_vs_%s_%s_chromVAR_chip_pseudobulk_volcano.pdf",i,j,args$motif_annotation)), width=7, height=5)
    print(p)
    dev.off()
  }
}

