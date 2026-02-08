source(here::here("settings.R"))
source(here::here("utils.R"))

library(enrichplot)

#####################
## Define settings ##
#####################

# I/O
io$diff_results <- file.path(io$basedir,"results/rna/differential/metacells/alias/parsed/diff_expr_results.txt.gz")
io$gsea_go_results <- file.path(io$basedir,"results/rna/gene_signatures/differential/gsea_go_list.rds")
io$gsea_msigdb_results <- file.path(io$basedir,"results/rna/gene_signatures/differential/gsea_mgsidb_list.rds")
io$outdir <- file.path(io$basedir,"results/rna/gene_signatures/differential/pdf"); dir.create(io$outdir, showWarnings = F)

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

#######################
## Load GSEA results ##
#######################

msigdb_results_list <- readRDS(io$gsea_msigdb_results)
go_results_list <- readRDS(io$gsea_go_results)

#############
## Plot GO ##
#############

# i <- "4F_8ON"; j <- "4F_8OFF"
for (i in opts$groups) {
  for (j in opts$groups) {
    if (i!=j) {
      p_list <- list()
      for (k in c("up","down")) {
        tmp <- go_results_list[[sprintf("%s_vs_%s",i,j)]][[k]]
        if (is.null(tmp)) { 
          to.plot <- data.table(Description="NA", Count=1, qvalue=1)
        } else {
          to.plot <- tmp@result %>% as.data.table %>% head(n=25) %>%
            setorder(-qvalue) %>%
            .[,Description:=gsub("HALLMARK_","",Description)] %>%
            .[,Description:=factor(Description, levels=Description)]
        }
        p_list[[k]] <- ggplot(to.plot, aes(x=Count, y=Description)) +
          geom_point(aes(fill=-log10(qvalue)), size=4, shape=21) +
          scale_fill_gradient(low = "gray80", high = "red") + 
          labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",k,j)) +
          theme_bw() +
          theme(
            axis.text.x = element_text(color="black"),
            axis.text.y = element_text(color="black", size=rel(0.80)),
          )
        
      }
      p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
      
      pdf(file.path(io$outdir,sprintf("GO_%s_vs_%s_BJ_multiome.pdf",i,j)), width=13, height=7)
      print(p)
      dev.off()
    }
  }
}


#################
## Plot msigdb ##
#################

# i <- "4F_8ON"; j <- "4F_8OFF"
for (i in opts$groups) {
  for (j in opts$groups) {
    if (i!=j) {
      p_list <- list()
      for (k in c("up","down")) {
        tmp <- msigdb_results_list[[sprintf("%s_vs_%s",i,j)]][[k]]
        if (is.null(tmp)) { 
          to.plot <- data.table(Description="NA", Count=1, qvalue=1)
        } else {
          to.plot <- tmp@result %>% as.data.table %>% head(n=25) %>%
            setorder(-qvalue) %>%
            .[,Description:=gsub("HALLMARK_","",Description)] %>%
            .[,Description:=factor(Description, levels=Description)]
        }
        p_list[[k]] <- ggplot(to.plot, aes(x=Count, y=Description)) +
          geom_point(aes(fill=-log10(qvalue)), size=4, shape=21) +
          scale_fill_gradient(low = "gray80", high = "red") + 
          labs(x="Number of diff. genes in the pathway", y="", title=sprintf("%s-regulated in %s",k,j)) +
          theme_bw() +
          theme(
            axis.text.x = element_text(color="black"),
            axis.text.y = element_text(color="black", size=rel(0.80)),
          )
        
      }
      p <- cowplot::plot_grid(plotlist=p_list, nrow=1)
      
      pdf(file.path(io$outdir,sprintf("msigdb_%s_vs_%s_BJ_multiome.pdf",i,j)), width=13, height=7)
      print(p)
      dev.off()
    }
  }
}
  

#####################################
## Scatterplots comparing 2F vs 4F ##
#####################################

# opts$conditions <- c(
#   # "8OFF",
#   "6OFF-2ON",
#   "4OFF-4ON",
#   "2OFF-2ON-4OFF",
#   "4ON-4OFF",
#   "8ON"
# )
# 
# background_condition <- "8OFF"
# 
# # i <- opts$conditions[1]
# for (i in opts$conditions) {
#   if (i!=j) {
# 
#     diff_1.dt <- diff.dt[conditionA==background_condition & conditionB==i & cell_lineA=="2F" & cell_lineB=="2F",c("feature","transcript_type","logFC")]
#     diff_2.dt <- diff.dt[conditionA==background_condition & conditionB==i & cell_lineA=="4F" & cell_lineB=="4F",c("feature","transcript_type","logFC")]
#     to.plot <- merge(diff_1.dt, diff_2.dt, by=c("feature","transcript_type"), suffixes=c("_A","_B")) %>%
#       .[,tmp:=abs(logFC_A)] %>% setorder(-tmp) %>% .[,tmp:=NULL] %>%
#       .[transcript_type=="gene",transcript_type:="Genes"] %>%
#       .[transcript_type=="te",transcript_type:="Transposable Elements"] %>%
#       .[transcript_type=="tf",transcript_type:="Transcription Factors"]
# 
#     foo <- to.plot[abs(logFC_A)>=0.5 | abs(logFC_B)>=0.5]
# 
#     ylim.max <- max(abs(to.plot$logFC_A))
#     xlim.max <- max(abs(to.plot$logFC_B))
# 
#     p <- ggplot(to.plot, aes(x=logFC_A, y=logFC_B)) +
#       # geom_jitter(size=1, width=0.01, height=0.01) +
#       geom_jitter(size=0.15, width=0.01, height=0.01, color="gray70", data=to.plot[!(abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5)]) +
#       geom_jitter(size=0.75, width=0.01, height=0.01, color="black", data=to.plot[abs(logFC_B)>=1.5 | abs(logFC_A)>=1.5]) +
#       geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", linewidth=0.1) +
#       geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", linewidth=0.1) +
#       ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Genes" & abs(logFC_B)>=2 & abs(logFC_A)>=2],n=50), aes(x=logFC_A, y=logFC_B, label=feature), color="black", size=2, segment.size=0.1, max.overlaps = Inf) +
#       ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Transposable Elements" & abs(logFC_B)>=1 & abs(logFC_A)>=1],n=25), aes(x=logFC_A, y=logFC_B, label=feature), color="black", size=2, segment.size=0.1, max.overlaps = Inf) +
#       ggrepel::geom_text_repel(data=head(to.plot[transcript_type=="Transcription Factors" & abs(logFC_B)>=1 & abs(logFC_A)>=1],n=25), aes(x=logFC_A, y=logFC_B, label=feature), color="black", size=2, segment.size=0.1, max.overlaps = Inf) +
#       ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_B)>=1.5 & abs(logFC_A)<=0.5],n=50), aes(x=logFC_A, y=logFC_B, label=feature), color="darkred", size=2, segment.size=0.1, max.overlaps = Inf) +
#       ggrepel::geom_text_repel(data=head(to.plot[abs(logFC_A)>=1.5 & abs(logFC_B)<=0.5],n=50), aes(x=logFC_A, y=logFC_B, label=feature), color="darkgreen", size=2, segment.size=0.1, max.overlaps = Inf) +
#       # labs(x=sprintf("logFC (%s vs %s, 2F)",i,background_condition), y=sprintf("logFC (%s vs %s, 4F)",i,background_condition), title=sprintf("%s vs %s",i,background_condition)) +
#       labs(x="logFC (2F)", y="logFC (4F)", title=sprintf("%s vs %s",i,background_condition)) +
#       # annotate("text", x=1, y=ylim/1.1, size=4, label="Down in L1-ASO") +
#       # annotate("text", x=1, y=-ylim/1.1, size=4, label="Up in L1-ASO") +
#       # coord_cartesian(xlim=c(-max.expr,max.expr), ylim=c(-max.expr,max.expr)) +
#       facet_wrap(~transcript_type, nrow=1) +
#       theme_classic() +
#       theme(
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(size=rel(1.0), color="black"),
#         axis.text.y = element_text(size=rel(1.0), color="black"),
#         legend.position = "none",
#         legend.title = element_blank()
#       )
# 
#     pdf(file.path(io$outdir.scatterplots,sprintf("scatterplot_%s_vs_%s_2F_vs_4f_BJ_multiome.pdf",i,background_condition)), width=11, height=7)
#     print(p)
#     dev.off()
#   }
# }


#########################################
## Compare to BJ reprogramming (CHAMP) ##
#########################################

# diff_bj_champ_d2_vs_d0.dt <- fread("/Users/rargelaguet/data/BJ_fibroblasts_reprogramming_scrna_altos/results/differential/metacells/group/diff_group_D02_Dox_vs_D00_Dox.txt.gz") %>%
#   setnames("gene","feature")

# features <- intersect(diff_bj_champ_d2_vs_d0.dt$feature, unique(diff.dt$feature))
# to.plot <- rbind(
#   diff.dt[feature%in%features & groupA=="D2" & groupB=="D0",c("feature","transcript_type","logFC")] %>% .[,dataset:="Xing2020"],
#   diff_bj_champ_d2_vs_d0.dt[feature%in%features,c("feature","transcript_type","logFC")] %>% .[,dataset:="CHAMP"]
# ) %>%  .[is.na(logFC),logFC:=0] %>%
#   dcast(feature+transcript_type~dataset, value.var="logFC")

# xlim.max <- max(abs(to.plot$CHAMP), na.rm=T)
# ylim.max <- max(abs(to.plot$Xing2020), na.rm=T)

# p <- ggplot(to.plot, aes(x=CHAMP, y=Xing2020)) +
#   # geom_point(size=0.30, width=0.01, height=0.01) +
#   geom_point(size=1, color="black", data=to.plot[abs(Xing2020)>=1.5 & abs(CHAMP)>=1.5]) +
#   ggrastr::geom_point_rast(size=0.15, color="gray70", data=to.plot[!(abs(Xing2020)>=1.5 & abs(CHAMP)>=1.5)]) +
#   # stat_smooth(method="lm", alpha=0.5, color="gray30") +
#   stat_cor(method = "pearson") +
#   geom_segment(x=-xlim.max, xend=xlim.max, y=0, yend=0, color="orange", size=0.1) +
#   geom_segment(x=0, xend=0, y=-ylim.max, yend=ylim.max, color="orange", size=0.1) +
#   ggrepel::geom_text_repel(data=to.plot[abs(Xing2020)>=1.5 & abs(CHAMP)>=1.5], aes(x=CHAMP, y=Xing2020, label=feature), color="black", size=2, segment.size=0.1, max.overlaps = Inf) +
#   labs(x="logFC (D2 vs D0, CHAMP)", y="logFC (D2 vs D0, Xing2020)", title="") +
#   theme_classic() +
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     axis.text.x = element_text(size=rel(1.0), color="black"),
#     axis.text.y = element_text(size=rel(1.0), color="black"),
#     legend.position = "none",
#     legend.title = element_blank()
#   )

# pdf(file.path(io$outdir,"scatterplot_Xing2020_vs_CHAMP_D0_vs_D2.pdf"), width=8, height=6)
# print(p)
# dev.off()
