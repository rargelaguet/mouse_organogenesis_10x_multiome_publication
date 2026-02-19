###############
## libraries ##
###############

# general stuff
library(SingleCellExperiment)
library(scater)
library(markdown)
library(R.utils)
library(data.table)
library(purrr)
library(DT)

# general viz
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)
require(patchwork)
require(pheatmap)
require(RColorBrewer)

# shiny
library(shiny)
library(shinythemes)
library(ggiraph)

#
library(irlba)

# font management
library(gdtools)
register_gfont(family="Open Sans")

#############################
## Define global variables ##
#############################

groups <- c(
  "Ctrl",
  # "Ctrl_OFF_None",
  # "Ctrl_12ON-14OFF_None",
  "Ctrl_4ON_CD13", 
  "Ctrl_8ON_CD13", 
  "Ctrl_12ON_CD13", 
  "Ctrl_4ON-14OFF_CD13",
  "Ctrl_8ON-14OFF_CD13", 
  "Ctrl_12ON-14OFF_CD13", 
  "Dox_4ON_CD13", 
  "Dox_4ON_CD70", 
  "Dox_8ON_CD13", 
  "Dox_8ON_CD70",
  "Dox_12ON_CD13", 
  "Dox_12ON_CD70", 
  "Dox_12ON_SSEA4", 
  "Dox_22ON_CD13", 
  # "Dox_22ON_CD70", 
  "Dox_22ON_SSEA4", 
  "Dox_4ON-14OFF_CD13", 
  "Dox_4ON-14OFF_CD70", 
  "Dox_8ON-14OFF_CD13", 
  "Dox_8ON-14OFF_CD70", 
  "Dox_12ON-14OFF_CD13", 
  "Dox_12ON-14OFF_CD70", 
  "Dox_12ON-14OFF_SSEA4"
)

groups_no_ctrl <- c(
  "Ctrl",
  # "Ctrl_OFF_None",
  # "Ctrl_12ON-14OFF_None",
  # "Ctrl_4ON_CD13", 
  # "Ctrl_8ON_CD13", 
  # "Ctrl_12ON_CD13", 
  # "Ctrl_4ON-14OFF_CD13",
  # "Ctrl_8ON-14OFF_CD13", 
  # "Ctrl_12ON-14OFF_CD13", 
  "Dox_4ON_CD13", 
  "Dox_4ON_CD70", 
  "Dox_8ON_CD13", 
  "Dox_8ON_CD70",
  "Dox_12ON_CD13", 
  "Dox_12ON_CD70", 
  "Dox_12ON_SSEA4", 
  "Dox_22ON_CD13", 
  # "Dox_22ON_CD70", 
  "Dox_22ON_SSEA4", 
  "Dox_4ON-14OFF_CD13", 
  "Dox_4ON-14OFF_CD70", 
  "Dox_8ON-14OFF_CD13", 
  "Dox_8ON-14OFF_CD70", 
  "Dox_12ON-14OFF_CD13", 
  "Dox_12ON-14OFF_CD70", 
  "Dox_12ON-14OFF_SSEA4"
)


databases_gsea <- c("GO", "MSigDB_Hallmarks", "MitoCarta", "EnrichR_GO_Molecular_Function_2023", "EnrichR_GO_Cellular_Component_2023", "EnrichR_GO_Biological_Process_2023", "EnrichR_Reactome_2022", "EnrichR_MSigDB_Hallmark_2020", "EnrichR_BioCarta_2016", "EnrichR_KEGG_2021_Human")

#####################
## Colour palettes ##
#####################

# gene_expr_palette_color <- scale_color_gradient(low = "gray90", high = "red")
# gene_expr_palette_fill <- scale_fill_gradient(low = "gray90", high = "red")

# gene_signature_palette_color <- scale_color_gradient(low = "gray80", high = "darkgreen")
# gene_signature_palette_fill <- scale_fill_gradient(low = "gray80", high = "darkgreen")


timepoints <- c("OFF", "4ON", "8ON", "12ON", "22ON", "4ON-14OFF", "8ON-14OFF", "12ON-14OFF")

cell_lines <- donors <- c("25F", "63F", "70M")

surface_markers <- c("None", "CD13", "CD70","SSEA4")

timepoints_colours <- c(
  "OFF" = "#C1CDCD",
  "4ON" = "#9ECAE1",
  "8ON" = "#6BAED6",
  "12ON" = "#3182BD",
  "22ON" = "#08519C",
  "4ON-14OFF" = "#BAE4B3",
  "8ON-14OFF" = "#31A354",
  "12ON-14OFF" = "#006D2C"
)

conditions_colours <- c(
  "Ctrl" = "#A1A1A1", 
  "Ctrl_CD13" = "#CD5555", 
  "Dox_CD13" = "#8B2323", 
  "Dox_CD70" = "#4F94CD", 
  "Dox_SSEA4" = "#4169E1"
)


###############
## Functions ##
###############

minmax.normalisation <- function(x) {
  return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

ggplot_theme_NoAxes <- function() {
  theme(
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
}

matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]

give.n <- function(x) { return(c(y = max(x), label = length(x))) }

sidebarPanel2 <- function (..., out = NULL, width = 4)  {
  div(class = paste0("col-sm-", width), 
      tags$form(class = "well", ...),
      out
  )
}
