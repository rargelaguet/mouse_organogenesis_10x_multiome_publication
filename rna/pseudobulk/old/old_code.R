
###########################################################################
## Calculate average expression as the average of log-transformed values ##
###########################################################################

# expr.dt <- unique(sce$celltype.mapped) %>% map(function(i) {
#   dt <- logcounts(sce[,sce$celltype.mapped==i]) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
#     melt(id.vars="rn") %>% setnames(c("symbol","cell","value")) %>%
#     .[,.(mean_expr=round(mean(value),3)),by="symbol"] %>%
#     .[,celltype:=i]
#   return(dt)
# }) %>% rbindlist
# 
# length(unique(expr.dt$symbol))
# length(unique(expr.dt$celltype))

#################################################################
## Calculate average expression using the average count values ##
#################################################################

# NOTE: NOT WORKING
# expr.dt <- unique(sce$celltype.mapped) %>% map(function(i) {
#   dt <- counts(sce[,sce$celltype.mapped==i]) %>% as.matrix %>% as.data.table(keep.rownames = T) %>%
#     melt(id.vars="rn") %>% setnames(c("symbol","cell","value")) %>%
#     .[,.(counts=sum(value), mean_counts=round(mean(value),3)),by="symbol"] %>%
#     .[,celltype:=i]
#   return(dt)
# }) %>% rbindlist
# 
# foo <-expr.dt %>%
#   .[,sum_counts:=sum(counts),by="celltype"] %>%
#   .[,.(mean_counts=unique(mean_counts), mean_counts2=counts/unique(sum_counts)),by="symbol"]
# 
# length(unique(expr.dt$symbol))
# length(unique(expr.dt$celltype))

##########
## Save ##
##########

# to.save <- expr.dt %>%
#   merge(gene_metadata[,c("symbol","ens_id")], all.x=T) %>%
#   setnames("symbol","gene")
# fwrite(to.save, paste0(io$outdir,"/avg_expr_per_celltype_and_gene.txt.gz"), sep="\t")
