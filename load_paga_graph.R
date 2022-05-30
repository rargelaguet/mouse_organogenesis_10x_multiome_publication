library(GGally)
library(network)
library(sna)
library(ggraph)
library(igraph)
library(tidygraph)

connectivity.mtx <- fread(io$paga.connectivity) %>%
  matrix.please %>% .[opts$celltypes,opts$celltypes]

df.coordinates <- fread(io$paga.coordinates) %>% 
  matrix.please %>% .[opts$celltypes,]

# Parse data
connectivity.mtx[connectivity.mtx<0.20] <- 0
connectivity.mtx[connectivity.mtx>=0.20] <- 1

# Create igraph object
igraph.paga <- graph_from_adjacency_matrix(connectivity.mtx, mode = "undirected")

# Create tbl_graph object
igraph.paga.tbl <- as_tbl_graph(igraph.paga) %>%
  activate(nodes) %>%
  mutate(celltype=rownames(connectivity.mtx)) %>%
  mutate(x=df.coordinates[,1]) %>% mutate(y=df.coordinates[,2])

# Create network object
net.paga = network(connectivity.mtx)
net.paga %v% "x" = connectivity.mtx[, 1]
net.paga %v% "y" = connectivity.mtx[, 2]

##########
## TEST ##
##########

# sum(connectivity.mtx==1)
# connectivity.mtx["Epiblast","Rostral_neurectoderm"]
# igraph.paga.tbl %>% activate(edges) %>% as.data.table()  %>% nrow
# filter(celltype=="Epiblast")
