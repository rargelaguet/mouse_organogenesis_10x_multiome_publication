source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$feature_metadata <- "/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/processed/atac/archR/PeakCalls/peak_metadata.tsv.gz"
io$diff_results_bj_reprogramming <- "/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/results/atac/archR/differential/metacells/alias/PeakMatrix/parsed/diff_results.txt.gz"
io$diff_results_bj_passaging <- "/Users/rargelaguet/data/Chan2022_senescence/results/atac/limma/groups/control_groups/limma_results_BJ_reprogramming_peaks.rds"
io$outdir <- file.path(io$basedir,"results/atac/compare_passaging_to_reprogramming"); dir.create(io$outdir, showWarnings = F)

# Options

###########################
## Load feature metadata ##
###########################

feature_metadata.dt <- fread(io$feature_metadata, select=c(1,2,3,5,13)) %>%
  setnames(c("chr","start","end","score","peak_type")) %>%
  .[,feature:=sprintf("%s:%s-%s",gsub("chr","",chr),start,end)] %>%
  .[,c("feature","chr","start","end","score","peak_type")] %>%
  .[!feature%in%duplicated.sites] %>% setkey(feature) %>% .[rownames(atac.mtx)]

############################################
## Load differential results BJ passaging ##
############################################

limma_results.list <- readRDS(io$diff_results_bj_passaging)
diff_passage.dt <- limma_results.list[["Passage"]] %>% as.data.table %>%
  merge(feature_metadata.dt, by="feature") %>%
  .[,feature:=sprintf("chr%s",feature)]

################################################
## Load differential results BJ reprogramming ##
################################################

diff_reprogramming.dt <- fread(io$diff_results_bj_reprogramming)
diff_reprogramming.dt$transcript_type <- NULL

#######################################################
## Scatterplots comparing passaging vs reprogramming ##
#######################################################

# # reprogramming.foreground.condition <- "4F_4OFF-4ON"; reprogramming.background.condition <- "4F_8OFF"
# reprogramming.foreground.condition <- "2F_4OFF-4ON"; reprogramming.background.condition <- "2F_8OFF"
# 
# diff_reprogramming_filt.dt <- diff_reprogramming.dt[groupA==reprogramming.background.condition & groupB==reprogramming.foreground.condition,c("feature","logFC")]
# diff_passaging_filt.dt <- diff_passage.dt[,c("feature","logFC")] %>% .[abs(logFC)>=0.05]
# 
# to.plot <- merge(diff_reprogramming_filt.dt, diff_passaging_filt.dt, by=c("feature"), suffixes=c("_A","_B")) %>%
#   .[,tmp:=abs(logFC_A)] %>% setorder(-tmp) %>% .[,tmp:=NULL]
# 
# foo <- to.plot[abs(logFC_A)>=0.5 | abs(logFC_B)>=0.5]
# 
# xlim.max <- max(abs(to.plot$logFC_A))
# ylim.max <- max(abs(to.plot$logFC_B))
# 
# ggplot(to.plot, aes(x=logFC_A, y=logFC_B)) +
#   geom_point(size=0.50, color="gray70") +
#   # geom_jitter(size=0.15, width=0.01, height=0.01, color="gray70", data=to.plot[!(abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5)]) +
#   # geom_jitter(size=0.75, width=0.01, height=0.01, color="black", data=to.plot[abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5]) +
#   geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
#   geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
#   # ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_B)>=1.5 & abs(logFC_A)<=0.5],n=50), aes(x=logFC_A, y=logFC_B, label=feature), color="darkred", size=2, segment.size=0.1, max.overlaps = Inf) +
#   # ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_A)>=1.5 & abs(logFC_B)<=0.5],n=50), aes(x=logFC_A, y=logFC_B, label=feature), color="darkgreen", size=2, segment.size=0.1, max.overlaps = Inf) +
#   labs(x="logFC (Reprogramming 4ON)", y="logFC (Passaging)", title="") +
#   # coord_cartesian(xlim=c(-max.expr,max.expr), ylim=c(-max.expr,max.expr)) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.text.x = element_text(size=rel(1.0), color="black"),
#     axis.text.y = element_text(size=rel(1.0), color="black"),
#     legend.position = "none",
#     legend.title = element_blank()
#   )
# 
# # pdf(file.path(io$outdir.scatterplots,sprintf("scatterplot_acc_%s_vs_%s_2F_vs_4F_BJ_multiome.pdf",i,background_condition)), width=11, height=7)
# # print(p)
# # dev.off()


##############################################################################################################
## Scatterplots comparing BJ passaging vs reprogrmaming, split peaks that increase vs decrease with passage ##
##############################################################################################################

# reprogramming.foreground.condition <- "4F_8ON"; reprogramming.background.condition <- "4F_8OFF"
# reprogramming.foreground.condition <- "4F_4OFF-4ON"; reprogramming.background.condition <- "4F_8OFF"
reprogramming.foreground.condition <- "4F_4ON-4OFF"; reprogramming.background.condition <- "4F_8OFF"
# reprogramming.foreground.condition <- "2F_8ON"; reprogramming.background.condition <- "2F_8OFF"
# reprogramming.foreground.condition <- "2F_4OFF-4ON"; reprogramming.background.condition <- "2F_8OFF"
# reprogramming.foreground.condition <- "2F_4ON-4OFF"; reprogramming.background.condition <- "2F_8OFF"

diff_passaging_filt.dt <- diff_passage.dt[,c("feature","logFC")] %>% .[abs(logFC)>=0.05] %>% .[,sign:=ifelse(logFC>0,"Open passaging","Close passaging")]
diff_reprogramming_filt.dt <- diff_reprogramming.dt[groupA==reprogramming.background.condition & groupB==reprogramming.foreground.condition,c("feature","logFC")]

to.plot <- merge(diff_passaging_filt.dt, diff_reprogramming_filt.dt , by=c("feature"), suffixes=c("_A","_B")) %>%
  .[,tmp:=abs(logFC_A)] %>% setorder(-tmp) %>% .[,tmp:=NULL]

xlim.max <- max(abs(to.plot$logFC_A))
ylim.max <- max(abs(to.plot$logFC_B))

ggplot(to.plot, aes(x=abs(logFC_A), y=logFC_B)) +
  geom_point(size=0.75, color="gray70") +
  # geom_jitter(size=0.15, width=0.01, height=0.01, color="gray70", data=to.plot[!(abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5)]) +
  # geom_jitter(size=0.75, width=0.01, height=0.01, color="black", data=to.plot[abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5]) +
  geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
  geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
  # ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_B)>=1.5 & abs(logFC_A)<=0.5],n=50), aes(x=logFC_A, y=logFC_B, label=feature), color="darkred", size=2, segment.size=0.1, max.overlaps = Inf) +
  # ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_A)>=1.5 & abs(logFC_B)<=0.5],n=50), aes(x=logFC_A, y=logFC_B, label=feature), color="darkgreen", size=2, segment.size=0.1, max.overlaps = Inf) +
  labs(x="logFC (Passaging)", y=sprintf("logFC (%s Reprogramming)",reprogramming.foreground.condition), title="") +
  facet_wrap(~sign, scales="fixed") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  )

