# i="T_789"; j <- "NMP"
tfs.to.use <- colnames(virtual_chip.mtx)
celltypes.to.use <- "NMP"
# i <- "T"
tf_enrichment_insilico_chip.dt <- tfs.to.use %>% map(function(i) { 
  print(i)
  
  celltypes.to.use %>% map(function(j) {
    
    foreground.peaks <- diff.dt[celltype==j & sign=="Downregulated in KO" & sig==TRUE,feature]
    background.peaks <- diff.dt[celltype==j & sig==FALSE,feature]
    foreground.nmatches <- virtual_chip_logical.mtx[foreground.peaks,i] %>% sum
    background.nmatches <- virtual_chip_logical.mtx[background.peaks,i] %>% sum
    
    p.value <- phyper(foreground.nmatches-1, background.nmatches, length(background.peaks)-background.nmatches,length(foreground.peaks), lower.tail = F)

    data.table(tf=i, celltype=j, pval=p.value)
  }) %>% rbindlist }) %>% rbindlist


to.plot <- tf_enrichment_insilico_chip.dt[pval<=0.10] %>% 
  .[,pval:=as.numeric(pval)] %>%
  .[,log_pval:=-log10(pval)] %>%
  .[,celltype:=factor(celltype,levels=celltypes.to.use)]

to.plot[,dot_size:=minmax.normalisation(abs(log_pval))]

to.plot.text <- to.plot[pval<=0.01]

ggplot(to.plot[pval<=0.01], aes_string(x="log_pval", y="tf", size="dot_size")) +
  geom_point(shape=21) +
  # scale_x_discrete(drop=F) +
  # scale_size_continuous(range = c(0.25,2)) +
  # guides(x = guide_axis(angle = 90)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color="black", size=rel(0.75)),
    axis.text.y = element_text(color="black")
  )
