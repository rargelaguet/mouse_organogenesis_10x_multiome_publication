# TO-DO: USE OUTPUT OF SAVE ATAC MATRICES
here::i_am("atac/archR/feature_stats/archR_calculate_feature_stats.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(sparseMatrixStats))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--feature_stats',    type="character",    help='feature stats file')
p$add_argument('--outdir',     type="character",    help='Output directory')
args <- p$parse_args(commandArgs(TRUE))

## START TEST ##
# args$feature_stats <- file.path(io$basedir, "results/atac/archR/feature_stats/PeakMatrix_feature_stats.txt.gz")
# args$outdir <- file.path(io$basedir, "results/atac/archR/feature_stats/PeakMatrix/pdf")
## END TEST ##

dir.create(args$outdir, showWarnings=F, recursive = T)

########################
## Load feature stats ##
########################

##########
## Plot ##
##########

to.plot <- peak_metadata.dt %>% 
  .[score>=args$min_peak_score] %>%
  merge(peak_stats.dt,by="peak")

p <- ggboxplot(to.plot, x="peakType", y="mean_singlecell", fill="peakType", outlier.shape = NA) +
  coord_cartesian(ylim=c(0,1.0)) +
  labs(x="", y="Average chromatin accessibility") +
  theme(
    axis.text = element_text(size=rel(0.75)),
    legend.position = "none"
  )

pdf(file.path(args$outdir,"boxplots_atac_peak_type.pdf"), width = 7, height = 5)
print(p)
dev.off()


##########
## TEST ##
##########

# to.plot <- feature_stats.dt %>% head(n=1e4)
# to.plot <- feature_stats.dt[var_pseudobulk>=0.10 & var_metacells<1]
# 
# ggscatter(to.plot, x="var_metacells", y="var_pseudobulk", size=1) + 
#   geom_abline(slope=1, intercept=0)
# 
# ggscatter(to.plot, x="var_cells", y="var_pseudobulk", size=1) + 
#   geom_abline(slope=1, intercept=0)
# 
# ggscatter(to.plot, x="mean_metacells", y="var_metacells", size=1) + 
#   stat_smooth(method="loess")
# 
# ggscatter(to.plot, x="mean_pseudobulk", y="var_pseudobulk", size=0.5) + 
#   stat_smooth(method="loess")

