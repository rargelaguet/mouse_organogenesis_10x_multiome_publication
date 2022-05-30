
plot.dimred <- function(plot_df, query.label, atlas.label = "Atlas") {
  
  # Define dot size  
  size.values <- c(opts$size.mapped, opts$size.nomapped)
  names(size.values) <- c(query.label, atlas.label)
  
  # Define dot alpha  
  alpha.values <- c(opts$alpha.mapped, opts$alpha.nomapped)
  names(alpha.values) <- c(query.label, atlas.label)
  
  # Define dot colours  
  colour.values <- c("red", "lightgrey")
  names(colour.values) <- c(query.label, atlas.label)
  
  # Plot
  ggplot(plot_df, aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = size.values) +
    scale_alpha_manual(values = alpha.values) +
    scale_colour_manual(values = colour.values) +
    # labs(x="UMAP Dimension 1", y="UMAP Dimension 2") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}

plot.dimred.wtko <- function(plot_df, wt.label = "WT", ko.label = "KO", nomapped.label = "-") {
  
  # Define dot size  
  size.values <- c(opts$size.mapped, opts$size.mapped, opts$size.nomapped)
  names(size.values) <- c(wt.label, ko.label, nomapped.label)
  
  # Define dot alpha  
  alpha.values <- c(opts$alpha.mapped, opts$alpha.mapped, opts$alpha.nomapped)
  names(alpha.values) <- c(wt.label, ko.label, nomapped.label)
  
  # Define dot colours  
  colour.values <- c("red", "blue", "lightgrey")
  names(colour.values) <- c(wt.label, ko.label, nomapped.label)
  
  # Plot
  ggplot(plot_df, aes(x=V1, y=V2)) +
    ggrastr::geom_point_rast(aes(size=mapped, alpha=mapped, colour=mapped)) +
    scale_size_manual(values = size.values) +
    scale_alpha_manual(values = alpha.values) +
    scale_colour_manual(values = colour.values) +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    theme_classic() +
    theme(
      legend.position = "top", 
      legend.title = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )
}
