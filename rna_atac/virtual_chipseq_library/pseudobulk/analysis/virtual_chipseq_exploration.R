# here::i_am("rna_atac/virtual_chipseq_library/virtual_chipseq_plot_stats.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# Options
opts$motif_annotation <- "CISBP"

opts$TFs <- c("TAL1", "GATA1", "RUNX1", "FOXA2", "GATA4", "CDX2","NKX2-5","TBX5", "SOX10")
opts$TFs <- c("T")

# I/O
io$virtual_chip.dir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s",opts$motif_annotation))
io$virtual_chip.mtx <- file.path(io$virtual_chip.dir,"virtual_chip_matrix.rds")
io$outdir <- file.path(io$basedir,sprintf("results/rna_atac/virtual_chipseq/pseudobulk/%s/test",opts$motif_annotation)); dir.create(io$outdir, showWarnings = F)

###################################
## Load virtual ChIP-seq library ##
###################################

# Load detailed data.tables
virtual_chip.dt <- opts$TFs %>% map(function(i) {
    fread(sprintf("%s/%s.txt.gz",io$virtual_chip.dir,i)) %>%
    # .[,c("chr","start","end"):=NULL] %>%
    .[,tf:=i] %>%
    return
}) %>% rbindlist

# Load matrix
virtual_chip.mtx <- readRDS(io$virtual_chip.mtx)


#######################
## Explore Brachyury ##
#######################

to.plot <- virtual_chip.dt[abs(score)>=0.20 & motif_score>=0.30] %>%
    .[,sign:=as.factor(c("Repressor","Activator")[(correlation_score>0)+1])]

ggbarplot(to.plot[,.N,by=c("sign")], x="sign", y="N", fill="gray70") +
  labs(x="", y="Number of in silico T binding events")
