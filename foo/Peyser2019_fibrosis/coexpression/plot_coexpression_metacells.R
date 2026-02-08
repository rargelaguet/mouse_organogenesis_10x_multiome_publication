# here::i_am("rna/dimensionality_reduction/dimensionality_reduction_sce.R")

library(scran)
library(scater)

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$sce <- file.path(io$basedir,"results/rna/metacells/scTE/SingleCellExperiment_metacells.rds")
io$metadata <- file.path(io$basedir,"results/rna/metacells/scTE/metacells_metadata.txt.gz")
io$correlation.mtx <- file.path(io$basedir,"results/rna/coexpression/correlation_matrix.rds")
io$outdir <- paste0(io$basedir,"/results/rna/coexpression/pdf"); dir.create(io$outdir, showWarnings = F)

# Options
opts$groups <- c(
  # Baseline
  # "25F_day0",
  # "28F_day0",
  # "52F_day0",
  # "54F_day0",
  "63F_day0",
  # "70M_day0",
  # "FSFA_day0",
  "FSFB_day0",
  # 63F OSKM
  "63F_OSKM-T3G_day2",
  "63F_OSKM-T3G_day7",
  "63F_OSKM-T3G_day13",
  "63F_OSKM-T3G_day27_unsorted",
  # BJ_CHAMP OSKM
  "BJ_OSKM-CHAMP_day0",
  "BJ_OSKM-CHAMP_day2",
  "BJ_OSKM-CHAMP_day7",
  "BJ_OSKM-CHAMP_day13",
  "BJ_OSKM-CHAMP_day27_FTR",
  # FSF OSKM
  "FSF_OSKM-T3G_day2",
  "FSF_OSKM-T3G_day7",
  "FSF_OSKM-T3G_day13",
  "FSF_OSKM-T3G_day27_FTR",
  "FSF_OSKM-T3G_day27_TR",
  "FSF_OSKM-T3G_day27_unsorted"
)

opts$cell_lines <- c("FSF", "63F")
opts$celltypes <- c("Fibroblast", "Refractory", "Initiation", "Maturation", "FTR_Reversion", "TR_Reversion", "unsorted_Reversion")

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[group%in%opts$groups] %>% .[,group:=factor(group,levels=opts$groups)] %>%
  .[cell_line=="63F_OSKM-T3G",cell_line:="63F"] %>% .[cell_line%in%c("FSF_OSKM-T3G","FSFA","FSFB"),cell_line:="FSF"] %>%
  .[celltype_liu2020=="re6_unknown",c("celltype_manual","group_manual"):=list("re6_unknown",sprintf("%s#re6_unknown",group))] %>% 
  .[celltype_manual=="Reversion",celltype_manual:="unsorted_Reversion"] %>% 
  .[celltype_manual%in%opts$celltypes] %>% .[,celltype_manual:=factor(celltype_manual,levels=opts$celltypes)] %>%
  .[,cell_line:=factor(cell_line,levels=opts$cell_lines)] %>%
  .[,day:=factor(day,levels=opts$days)] %>%
  setnames("metacell","cell") %>%
  droplevels

table(cell_metadata.dt$cell_line)
table(cell_metadata.dt$celltype_manual)
table(cell_metadata.dt$group_manual)

########################
## Load gene metadata ##
########################

gene_metadata.dt <- fread(io$gene_metadata)

TEs <- fread(io$TEs, header = F)[[1]]

TFs <- fread(io$TFs)[[1]]

surface_markers <- unique(fread(io$cell_surface_markers)[[1]])

#############################################################
## Load RNA expression data as SingleCellExperiment object ##
#############################################################

sce <- load_SingleCellExperiment(io$sce, cells=cell_metadata.dt$cell, normalise = T)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

###############################
## Load correlation matrices ##
###############################

correlation.mtx <- readRDS(io$correlation.mtx)

correlation_liu2020.mtx <- readRDS(file.path(io$basedir_liu2020,"results/rna/coexpression/noNaive/correlation_matrix.rds"))

correlation_bj_champ.mtx <- readRDS(file.path(io$basedir_bj_scrna,"results/coexpression/correlation_matrix.rds"))

####################################################
## Plot correlation of HERVK with surface markers ##
####################################################

# i <- "HERVK-int"
i <- "LTR7Y"

cor.dt <- data.table(
  gene = colnames(correlation.mtx),
  cor = correlation.mtx[i,],
  abs_cor = abs(correlation.mtx[i,])
) %>% 
  .[gene!=i] %>%
  .[,transcript_type:="Genes"] %>%
  .[gene%in%TEs,transcript_type:="Transposable Elements"] %>%
  .[gene%in%TFs,transcript_type:="Transcription Factors"] %>%
  .[gene%in%surface_markers,transcript_type:="Surface markers"] %>%
  setorder(-abs_cor)

cor.dt <- cor.dt[transcript_type%in%c("Transcription Factors","Surface markers")]

to.plot.pos <- cor.dt[cor>0,.SD[1:20], by="transcript_type"] %>%
  .[,gene:=factor(gene,levels=rev(gene))]

# p.pos <- ggplot(to.plot.pos, aes(x=gene, y=abs_cor)) +
#   geom_point(size=1) +
#   facet_wrap(~transcript_type, nrow=1, scales="free_x") +
#   labs(x="", y="Correlation coefficient", title="Positive correlations") +
#   guides(x = guide_axis(angle = 90)) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust=0.5),
#     axis.text.x = element_text(size=rel(0.75), color="black"),
#     axis.text.y = element_text(size=rel(1.0), color="black"),
#     legend.position = "right",
#     legend.title = element_blank()
#   )

p.pos <- ggplot(to.plot.pos, aes(x=abs_cor, y=gene)) +
  geom_point(size=1) +
  geom_segment(x=0, aes(xend=abs_cor, y=gene, yend=gene),  linewidth=0.3) +
  facet_wrap(~transcript_type, nrow=1, scales="free_y") +
  labs(x="", y="Correlation coefficient", title=sprintf("(+) correlations with %s",i)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(0.90), color="black"),
    legend.position = "right",
    legend.title = element_blank()
  )


to.plot.neg <- cor.dt[cor<0,.SD[1:20], by="transcript_type"] %>%
  .[,gene:=factor(gene,levels=rev(gene))]

# p.neg <- ggplot(to.plot.neg, aes(x=gene, y=abs_cor)) +
#   geom_point(size=1) +
#   facet_wrap(~transcript_type, nrow=1, scales="free_x") +
#   labs(x="", y="Correlation coefficient", title="Negative correlations") +
#   guides(x = guide_axis(angle = 90)) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust=0.5),
#     axis.text.x = element_text(size=rel(0.75), color="black"),
#     axis.text.y = element_text(size=rel(1.0), color="black"),
#     legend.position = "right",
#     legend.title = element_blank()
#   )

p.neg <- ggplot(to.plot.neg, aes(x=abs_cor, y=gene)) +
  geom_point(size=1) +
  geom_segment(x=0, aes(xend=abs_cor, y=gene, yend=gene),  linewidth=0.3) +
  facet_wrap(~transcript_type, nrow=1, scales="free_y") +
  labs(x="", y="Correlation coefficient", title=sprintf("(-) correlations with %s",i)) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust=0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(0.80), color="black"),
    legend.position = "right",
    legend.title = element_blank()
  )

cowplot::plot_grid(plotlist=list(p.pos,p.neg), nrow = 2)

####################################
## Scatterplot of gene expression ##
####################################

# i <- "HERVK-int"; j <- "LTR5_Hs"
i <- "ANPEP"; j <- "ST6GAL1"
# i <- "HERVK-int"; j <- "CD55"
# i <- "LTR7Y"; j <- "HERVK-int"

to.plot <- data.table(
  geneA = i,
  geneB = j,
  expr_geneA = logcounts(sce)[i,],
  expr_geneB = logcounts(sce)[j,],
  cell_line = sce$cell_line,
  celltype_manual = sce$celltype_manual,
  group_manual = sce$group_manual
)

p <- ggplot(to.plot, aes(x=expr_geneA, y=expr_geneB)) +
  geom_point(aes(color=celltype_manual), size=1.5) +
  stat_smooth(method="lm", color="gray50") +
  coord_cartesian(ylim=c(0,max(to.plot$expr_geneB))) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  stat_cor(method = "pearson") +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(~cell_line, nrow=1) +
  labs(x=sprintf("%s expression (scMPTR)",i), y=sprintf("%s expression (scMPTR)",j)) +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(1.15), color="black"),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    legend.position = "right",
    legend.title = element_blank()
  )

pdf(file.path(io$outdir,sprintf("%s_vs_%s_coexpression_scatter_metacells_scMPTR.pdf",i,j)), width=10.5, height=5.75)
print(p)
dev.off()

########################################
## UpSet plots (Liu2020 and BJ CHAMP) ##
########################################

library(UpSetR)

# i <- "HERVK-int"
i <- "LTR7Y"

# scMPTR
cor_scmptr.dt <- data.table(
  gene = colnames(correlation.mtx),
  cor = correlation.mtx[i,],
  abs_cor = abs(correlation.mtx[i,])
) %>% setorder(-abs_cor) %>%
  .[gene!=i] %>%
  .[gene%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[transcript_type%in%c("Surface markers")]

# Liu2020
cor_liu2020.dt <- data.table(
  gene = colnames(correlation_liu2020.mtx),
  cor = correlation_liu2020.mtx[i,],
  abs_cor = abs(correlation_liu2020.mtx[i,])
) %>% setorder(-abs_cor) %>%
  .[gene!=i] %>%
  .[gene%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[transcript_type%in%c("Surface markers")]
  

# BJ CHAMP
cor_bj_champ.dt <- data.table(
  gene = colnames(correlation_bj_champ.mtx),
  cor = correlation_bj_champ.mtx[i,],
  abs_cor = abs(correlation_bj_champ.mtx[i,])
) %>% setorder(-abs_cor) %>%
  .[gene!=i] %>%
  .[gene%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[transcript_type%in%c("Surface markers")]



to.plot.upset <- list()
to.plot.upset[["scMPTR"]] <- cor.dt[cor>0.40] %>% head(n=50) %>% .$gene
to.plot.upset[["Liu2020"]] <- cor_liu2020.dt[cor>0.40] %>% head(n=50) %>% .$gene
to.plot.upset[["BJ_CHAMP"]] <- cor_bj_champ.dt[cor>0.40] %>% head(n=50) %>% .$gene

upset(fromList(to.plot.upset), nsets = length(to.plot.upset), order.by = "freq")

library(ggVennDiagram)

p <- ggVennDiagram(to.plot.upset) +
  scale_fill_distiller(palette = "Greys", direction = 1) +
  theme(
    legend.position = "none"
  )


pdf(file.path(io$outdir,"venn_diagram_coexpression_scMPTR.pdf"), width=8, height=5)
print(p)
dev.off()


dput(Reduce("intersect",to.plot.upset))

#################
## Scatterplot ##
#################

library(UpSetR)

# i <- "HERVK-int"
i <- "LTR7Y"

# scMPTR
cor_scmptr.dt <- data.table(
  feature = colnames(correlation.mtx),
  cor = correlation.mtx[i,]
  # dataset = "scMPTR"
) %>%
  .[feature!=i] %>%
  .[feature%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[transcript_type%in%c("Surface markers")]

# Liu2020
cor_liu2020.dt <- data.table(
  feature = colnames(correlation_liu2020.mtx),
  cor = correlation_liu2020.mtx[i,]
  # dataset = "Liu2020"
) %>%
  .[feature!=i] %>%
  .[feature%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[transcript_type%in%c("Surface markers")]

# BJ CHAMP
cor_bj_champ.dt <- data.table(
  feature = colnames(correlation_bj_champ.mtx),
  cor = correlation_bj_champ.mtx[i,]
  # dataset = "bj_champ"
) %>%
  .[feature!=i] %>%
  .[feature%in%surface_markers,transcript_type:="Surface markers"] %>%
  .[transcript_type%in%c("Surface markers")]

# to.plot <- merge(cor_scmptr.dt, cor_liu2020.dt, by=c("feature","transcript_type"), suffixes=c("_A","_B")) %>%
to.plot <- merge(cor_scmptr.dt, cor_bj_champ.dt, by=c("feature","transcript_type"), suffixes=c("_A","_B")) %>%
  .[,tmp:=abs(cor_A)] %>% setorder(-tmp) %>% .[,tmp:=NULL]
ylim.max <- max(abs(to.plot$cor_A))
xlim.max <- max(abs(to.plot$cor_B))

ggplot(to.plot, aes(x=cor_A, y=cor_B)) +
  geom_jitter(size=0.70, width=0.01, height=0.01, color="gray40") +
  stat_cor(method = "pearson") +
  stat_smooth(method="lm", color="gray30", alpha=0.30) +
  # geom_jitter(size=0.15, width=0.01, height=0.01, color="gray30", data=to.plot[!(abs(cor_B)>=0.40 | abs(cor_A)>=0.40)]) +
  # geom_jitter(size=0.45, width=0.01, height=0.01, color="gray30", data=to.plot[abs(cor_B)>=1.5 | abs(cor_A)>=1.5]) +
  geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
  geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
  ggrepel::geom_text_repel(data=head(to.plot[abs(cor_B)>=0.35 & abs(cor_A)>=0.35],n=20), aes(x=cor_A, y=cor_B, label=feature), color="black", size=2, segment.size=0.1, max.overlaps = Inf) +
  # labs(x=sprintf("Correlation with %s (scMPTR)",i), y=sprintf("Correlation with %s (Liu2020)",i), title="") +
  labs(x=sprintf("Correlation with %s (scMPTR)",i), y=sprintf("Correlation with %s (BJ CHAMP)",i), title="") +
  coord_cartesian(xlim=c(-max(xlim.max,ylim.max),max(xlim.max,ylim.max)), ylim=c(-max(xlim.max,ylim.max),max(xlim.max,ylim.max))) +
  facet_wrap(~transcript_type, nrow=1) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    legend.position = "none",
    legend.title = element_blank()
)

pdf(file.path(io$outdir,sprintf("%s_vs_%s_coexpression_scatter_metacells_scMPTR.pdf",i,j)), width=10.5, height=5.75)
print(p)
dev.off()
