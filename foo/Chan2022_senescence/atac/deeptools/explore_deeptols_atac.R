foo <- fread("/Users/rargelaguet/config/conda_environments/genome_browse/test.gz", skip = 1)

foo[1:5,1:15]
dim(foo)

dim(foo)



set.seed(123)
mat = matrix(rnorm(100), 10)
rownames(mat) = paste0("R", 1:10)
colnames(mat) = paste0("C", 1:10)

# column_ha = HeatmapAnnotation(foo1 = runif(10), bar1 = anno_barplot(runif(10)))
column_ha = HeatmapAnnotation(foo = 1:10)

row_ha = rowAnnotation(foo2 = runif(10), bar2 = anno_barplot(runif(10)))

Heatmap(mat, name = "mat", top_annotation = column_ha, right_annotation = row_ha)


atac.dt <- fread("/Users/rargelaguet/data/fibros_reik_lab/results/atac/archR/pseudobulk_exploration/sample_gapdh/test_multiBigwigSummary.tab", skip = 1)
colnames(atac.dt) <- c("chr", "start", "end", "25F_gapdh_high", "25F_gapdh_low", "28F_gapdh_high", "28F_gapdh_low", "52F_gapdh_high", 
                      "52F_gapdh_low", "54F_gapdh_high", "54F_gapdh_low", "63F_gapdh_high", "63F_gapdh_low", "70M_gapdh_high", 
                      "70M_gapdh_low", "FSFA_gapdh_high", "FSFB_gapdh_high")

to.plot <- atac.dt %>% melt(id.vars=c("chr","start","end"), variable.name="sample")

to.plot[,cell_line:=strsplit(as.character(sample),split = "_") %>% map_chr(1)]
to.plot[,population:=strsplit(as.character(sample),split = "_") %>% map_chr(~ paste(.[c(2, 3)], collapse = "_"))]

to.plot[,cell_line:=factor(cell_line,levels=c("FSFA", "FSFB", "25F", "28F", "52F", "54F", "63F", "70M"))]

ggplot(to.plot, aes(x=population, y=value, fill=population)) +
  geom_boxplot(outlier.shape=NA) +
  facet_wrap(~cell_line, scales = "fixed", nrow=1) +
  # scale_fill_manual(values=opts$celltype.colors) +
  coord_cartesian(ylim=c(0,0.005)) +
  labs(x="", y="ATAC signal") +
  theme_classic() +
  theme(
    axis.text.y = element_text(size=rel(0.90), color="black"),
    axis.text.x = element_text(colour="black",size=rel(1.0), angle=90, hjust=1, vjust=0.5),
    # axis.title = element_text(colour="black",size=rel(0.75)),
    # axis.ticks.x = element_blank(),
    legend.position = "right"
  )


tmp <- atac.dt[,c("chr","start","end","63F_gapdh_high","63F_gapdh_low")] %>%
  .[,diff:=`63F_gapdh_high`-`63F_gapdh_low`]


################################################################################################
################################################################################################
################################################################################################
################################################################################################


atac.dt <- fread("/Users/rargelaguet/data/adult_fibroblast_reprogramming_10x_multiome_reik/results/atac/archR/deeptools_pseudobulk/multiBigwigSummary_chromHMM.tab", skip = 1)
colnames(atac.dt) <- c("chr", "start", "end", "25F",  "25F_day0_GAPDH_high",  "25F_day0_GAPDH_low",  "28F",  "28F_day0_GAPDH_high",  "28F_day0_GAPDH_low",  "52F",  "52F_day0_GAPDH_high",  "52F_day0_GAPDH_low",  "54F",  "54F_day0_GAPDH_high",  "54F_day0_GAPDH_low",  "63F",  "63F_day0_GAPDH_high",  "63F_day0_GAPDH_low",  "63F_day13_Maturation",  "63F_day13_Refractory",  "63F_day27_Unsorted_Reversion",  "63F_day2_Initiation",  "63F_day2_Refractory",  "63F_day7_Initiation",  "63F_day7_Maturation",  "63F_day7_Refractory",  "63F_day7_re6_unknown",  "70M",  "70M_day0_GAPDH_high",  "70M_day0_GAPDH_low",  "BJ_day0_NoReprogramming",  "BJ_day13_NoReprogramming",  "BJ_day27_Unsorted_Reversion",  "BJ_day2_NoReprogramming",  "BJ_day2_Reprogramming",  "BJ_day7_NoReprogramming",  "FSFA",  "FSFA_day0_GAPDH_high",  "FSFB",  "FSFB_day0_GAPDH_high",  "FSFB_day0_GAPDH_low",  "FSF_day13_Refractory",  "FSF_day27_FTR_Reversion",  "FSF_day27_TR_Reversion",  "FSF_day27_Unsorted_Reversion",  "FSF_day2_Initiation",  "FSF_day2_Refractory",  "FSF_day7_Initiation",  "FSF_day7_Maturation",  "FSF_day7_Refractory",  "FSF_day7_re6_unknown")
atac.dt[,idx:=sprintf("%s:%s-%s",chr,start,end)]

samples.to.plot <- c("25F", "28F", "52F", "54F", "63F", "70M", "FSFA", "FSFB")

chromHMM.dt <- fread("/Users/rargelaguet/data/hg38_regulation/chromHMM/HDF_6_dense.merged.bed", select=c(1,2,3,4)) %>%
  setnames(c("chr","start","end","annotation"))
chromHMM.dt[,idx:=sprintf("chr%s:%s-%s",chr,start,end)]
chromHMM.dt[,c("chr","start","end"):=NULL]
atac.dt <- atac.dt %>% merge(chromHMM.dt, by="idx")

to.plot <- atac.dt %>%
  .[,c("chr","start","end","annotation",samples.to.plot), with = FALSE] %>% 
  melt(id.vars=c("chr","start","end","annotation"), variable.name="sample")

tmp <- atac.dt[,c("chr","start","end","25F","70M")] %>%
  .[,diff:=`70M`-`25F`]

