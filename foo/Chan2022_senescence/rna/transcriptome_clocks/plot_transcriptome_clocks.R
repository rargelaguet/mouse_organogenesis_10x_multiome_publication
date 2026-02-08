source(here::here("settings.R"))
source(here::here("utils.R"))

library(pheatmap)

#####################
## Define settings ##
#####################

# I/O
# io$metadata <- file.path(io$basedir,"processed/rna/sample_metadata.txt")
io$clock_results <- file.path(file.path(io$basedir,"rna/results/transcriptome_clocks/Senescence_Chan2022_meta_with_ProbPredict_5394genes_Ricard.csv"))
io$outdir <- file.path(io$basedir,"rna/results/transcriptome_clock"); dir.create(io$outdir, showWarnings = F, recursive = T)

# Options
opts$conditions <- c(
  "hTERT",
  "Replicative_Senescence",
  "Cell_Density",
  "Radiation_Senescence"
)

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>%
  .[condition%in%opts$conditions] %>% .[,condition:=factor(condition,levels=opts$conditions)] %>%
  .[time%in%opts$timepoints] %>% .[,time:=factor(time,levels=opts$timepoints)]

########################################
## Load transcriptional clock results ##
########################################

clock_results.dt <- fread(io$clock_results)[,c("alias","ProbPredict(Fleischer+GTEx)")] %>%
  setnames(c("alias","predicted_age"))

##########
## Plot ##
##########


to.plot <- sample_metadata.dt %>% merge(clock_results.dt, by="alias")

to.plot.means <- to.plot[,.(predicted_age=mean(predicted_age), sd=sd(predicted_age)), by=c("condition","time")]

p <- ggplot(to.plot, aes(x=time, y=predicted_age)) +
  geom_line(aes(group=condition), color="gray40", linewidth=1, alpha=0.75, data=to.plot.means) +
  geom_point(color="gray40", size=2.50) +
  geom_errorbar(aes(ymin=predicted_age-sd, ymax=predicted_age+sd), width=0.4, alpha=0.9, linewidth=0.5, data=to.plot.means) +
  facet_wrap(~condition, nrow=2, scales="free_x") +
  labs(x="Time", y="Predicted age", title="") +
  theme_classic() +
  theme(
    strip.text = element_text(size=rel(1.10)),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(colour="black",size=rel(1.0)),
    axis.title = element_text(colour="black",size=rel(1.10)),
    axis.text.y = element_text(colour="black",size=rel(1.10)),
    legend.position = "none",
    legend.title = element_blank(),
  )

pdf(file.path(io$outdir,"predicted_age_probpredict_Chan2022_senescence.pdf"), width=8.5, height=7)
print(p)
dev.off()
