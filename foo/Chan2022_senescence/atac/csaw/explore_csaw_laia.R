library(edgeR)
library(DESeq2)

csaw.se <- readRDS("/Users/rargelaguet/data/Chan2022_senescence/processed/atac/csaw/windows_width5000_spacing2500.rds")

colData(csaw.se)$sample <- basename(colData(csaw.se)$bam.files ) %>% gsub(".bam","",.)

metadata.dt <- data.table(sample=colData(csaw.se)$sample) %>%
  .[,passage:=strsplit(sample,"_") %>% map_chr(1)] %>%
  .[,condition:=strsplit(sample,"_") %>% map_chr(2)]

tmp <- as.data.table(rowRanges(csaw.se))
colnames(assay(csaw.se)) <- basename(colData(csaw.se)$bam.files ) %>% gsub(".bam","",.)
rownames(assay(csaw.se)) <- sprintf("%s:%s-%s",tmp$seqnames,tmp$start,tmp$end)

# Compute average log2 counts-per-million for each row of counts.
abundance <- edgeR::aveLogCPM(csaw.se, normalized.lib.sizes = T, prior.count = 1)

p <- ggplot(as.data.frame(abundance), mapping = aes(x = abundance)) +
  geom_density() +
  geom_vline(xintercept = -1.75, linetype = 1, linewidth = 1, colour = "tomato") +
  theme_classic() + xlim(-3.5,6)

pdf(paste0(work_dir,"window.counts/graphs/win.2kb.abundance.pdf"), width = 4, height = 4)
print(p)
dev.off()

csaw.se <- csaw.se[abundance>(-2),]


#######################################
## Differential analysis with DESeq2 ##
#######################################
  
counts <- csaw::asDGEList(csaw.se)$counts

# colnames(counts) <- unlist(lapply(csaw.se$bam.files, function(x){strsplit(x, split = "/", fixed = T)[[1]][3]}))
colnames(counts) <- colData(csaw.se)$sample

colData <- DataFrame(metadata.dt)
rownames(colData) <- colData$sample

sf <- csaw::asDGEList(csaw.se)$samples$norm.factors

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

dds$sizeFactor <- sf

# dds$time_point <- stats::relevel(x = dds$time_point, ref = "T2")

dds <- DESeq2::DESeq(dds)

res <- as.data.frame(DESeq2::results(dds))

# Add window coordinates information to the res object:
res <- cbind(res, data.frame(csaw.se@rowRanges))

##########################################################################################
## Stratify windows into deciles according to their mean accessibility in control cells ##
##########################################################################################

res$mean.ctrl <- apply(DESeq2::counts(dds, normalized = TRUE)[,colData[colData$condition == "Control",]$sample],1,mean)

breaks <- quantile(res$mean.ctrl, probs = seq(0, 1, 0.1))
breaks[1] <- 0

res$quantile <- cut(x = res$mean.ctrl, breaks)

# saveRDS(object = res, file = paste0(work_dir,"window.counts/res.diff.ATACSeq.win.2kb.T2-T3.RDS"))
# res <- readRDS(paste0(work_dir,"window.counts/res.diff.ATACSeq.win.2kb.T2-T3.RDS")) 


# Plot the mean ATAC-Seq counts in control cells divided into deciles. Subset the res object for representation:

subset.res <- res[sample(x = 1:nrow(res), size = round(nrow(res)*0.10, digits = 0)),]

p <- ggplot(data = subset.res, mapping = aes(x = quantile, y = log2(mean.ctrl + 1), fill = quantile)) +
  geom_boxplot(width = 0.8, notch = T, alpha = 0.7, outliers = F) +
  scale_fill_manual(values = colorRampPalette(colors = c("#FFF5F0","#FB6A4A","#67000D"))(10)) +
  ylim(2.75,9.5) +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(paste0(work_dir,"window.counts/graphs/boxplot.ATAC-Seq.T2-T3.mean.ctrl.in.deciles.pdf"), width = 6, height = 4)
print(p)
dev.off()


# Plot the fold change in accessibility of windows divided into deciles. Subset the res object for representation:
  
p <- ggplot(data = subset.res, mapping = aes(x = quantile, y = log2FoldChange, fill = quantile)) +
  geom_violin() +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", notch = T) +
  scale_fill_manual(values = colorRampPalette(colors = c("#FFF5F0","#FB6A4A","#67000D"))(10)) +
  geom_hline(yintercept = 0, linetype = 3, linewidth = 0.7) +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(paste0(work_dir,"window.counts/graphs/violin.ATAC-Seq.T2-T3.FC.T3.vs.T2.deciles.pdf"), width = 6, height = 4)
print(p)
dev.off()