source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$feature_metadata <- "/Users/rargelaguet/data/Chan2022_senescence/processed/atac/jack/peak_calls.bed"
# io$diff_results_bj_reprogramming <- "/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/results/atac/deeptools_pseudobulk/differential/diff_2F.txt"
io$diff_results_bj_reprogramming <- "/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_10x_multiome_CHAMP/results/atac/bigwig/Chan2022_senescence/differential/diff_4F.txt"
io$diff_results_bj_passaging <- "/Users/rargelaguet/data/Chan2022_senescence/results/atac/limma/samples/control_samples/limma_results.rds"
io$outdir <- file.path(io$basedir,"results/atac/compare_passaging_to_reprogramming"); dir.create(io$outdir, showWarnings = F)

# Options

###########################
## Load feature metadata ##
###########################

feature_metadata.dt <- fread(io$feature_metadata, select=c(1,2,3)) %>%
  setnames(c("chr","start","end")) %>%
  .[,feature:=sprintf("%s:%s-%s",chr,start,end)] %>%
  .[,c("feature","chr","start","end")]

############################################
## Load differential results BJ passaging ##
############################################

limma_results.list <- readRDS(io$diff_results_bj_passaging)
diff_passage.dt <- limma_results.list[["Passage"]] %>% as.data.table %>%
  merge(feature_metadata.dt, by="feature") 

stopifnot(unique(diff_passage.dt$feature)%in%feature_metadata.dt$feature)

############################################
## Load differential results BJ + Tg ##
############################################

io$diff_results_bj <- "/Users/rargelaguet/data/Chan2022_senescence/results/atac/limma/samples/all_samples/limma_results_ATAC_peaks.rds"

limma_results.list <- readRDS(io$diff_results_bj)
diff_tg.dt <- limma_results.list[["Treatment"]] %>% as.data.table %>%
  merge(feature_metadata.dt, by="feature") 

stopifnot(unique(diff_passage.dt$feature)%in%feature_metadata.dt$feature)

################################################
## Load differential results BJ reprogramming ##
################################################

diff_reprogramming.dt <- fread(io$diff_results_bj_reprogramming)
stopifnot(unique(diff_reprogramming.dt$feature)%in%feature_metadata.dt$feature)

##########################################################
## Scatterplots comparing Tg treatment vs reprogramming ##
##########################################################

# reprogramming.foreground.condition <- "2F_reprogramming_2ON"; reprogramming.background.condition <- "2F_unresponsive"
# reprogramming.foreground.condition <- "2F_reprogramming_4ON"; reprogramming.background.condition <- "2F_unresponsive"
# reprogramming.foreground.condition <- "2F_unresponsive_dividing"; reprogramming.background.condition <- "2F_unresponsive"
reprogramming.foreground.condition <- "4F_reprogramming_2ON"; reprogramming.background.condition <- "4F_unresponsive"
# reprogramming.foreground.condition <- "4F_reprogramming_4ON"; reprogramming.background.condition <- "4F_unresponsive"
# reprogramming.foreground.condition <- "4F_unresponsive_dividing"; reprogramming.background.condition <- "4F_unresponsive"

diff_tg_filt.dt <- diff_tg.dt[,c("feature","logFC")] %>% .[abs(logFC)>=0.03] %>% .[,sign:=ifelse(logFC>0,"Open TG","Close TG")]
diff_reprogramming_filt.dt <- diff_reprogramming.dt[group_A==reprogramming.background.condition & group_B==reprogramming.foreground.condition,c("feature","diff")] %>% setnames("diff","logFC")

to.plot <- merge(diff_tg_filt.dt, diff_reprogramming_filt.dt , by=c("feature"), suffixes=c("_A","_B")) %>%
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
  # labs(x="logFC (TG)", y=sprintf("logFC (%s Reprogramming)",reprogramming.foreground.condition), title="") +
  labs(x="logFC (TG)", y="logFC (4F Reprogramming)", title="") +
  facet_wrap(~sign, scales="fixed") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  )


ggplot(to.plot, aes(x=logFC_A, y=logFC_B)) +
  geom_point(size=0.55, color="gray70") +
  geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
  geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
  labs(x="logFC (TG)", y="logFC (4F Reprogramming)", title="") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size=rel(1.0), color="black"),
    axis.text.y = element_text(size=rel(1.0), color="black"),
    legend.position = "none",
    legend.title = element_blank()
  )

##############################################################################################################
## Scatterplots comparing BJ passaging vs reprogrmaming, split peaks that increase vs decrease with passage ##
##############################################################################################################

# reprogramming.foreground.condition <- "2F_reprogramming_2ON"; reprogramming.background.condition <- "2F_unresponsive"
# reprogramming.foreground.condition <- "2F_reprogramming_4ON"; reprogramming.background.condition <- "2F_unresponsive"
# reprogramming.foreground.condition <- "2F_unresponsive_dividing"; reprogramming.background.condition <- "2F_unresponsive"
# reprogramming.foreground.condition <- "4F_reprogramming_2ON"; reprogramming.background.condition <- "4F_unresponsive"
# reprogramming.foreground.condition <- "4F_reprogramming_4ON"; reprogramming.background.condition <- "4F_unresponsive"
reprogramming.foreground.condition <- "4F_unresponsive_dividing"; reprogramming.background.condition <- "4F_unresponsive"

diff_passaging_filt.dt <- diff_passage.dt[,c("feature","logFC")] %>% .[abs(logFC)>=0.03] %>% .[,sign:=ifelse(logFC>0,"Open passaging","Close passaging")]
diff_reprogramming_filt.dt <- diff_reprogramming.dt[group_A==reprogramming.background.condition & group_B==reprogramming.foreground.condition,c("feature","diff")] %>% setnames("diff","logFC")

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

