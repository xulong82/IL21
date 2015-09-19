# Copyright: Xulong Wang (xulong.wang@jax.org)

library(xlsx)
library(KEGG.db)
library(org.Mm.eg.db)
library(igraph)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

# PARSE THE IREGULON OUTPUT
load(file = "./data/profile12.rdt")
load(file = "./data/myTpm.rdt")
emean <- data.frame("NN"=rowMeans(myTpm[1:2]), "NP"=rowMeans(myTpm[3:4]), "PP"=rowMeans(myTpm[5:6]))
maxId <- apply(emean, 1, which.max)

# *** OPTIONAL: N
file <- "Pathway/iregulon_naive.tsv"
univ <- rownames(myTpm)[emean$NN > 30 & maxId == 1]
# *** OPTIONAL: ACT
file <- "Pathway/iregulon_act.tsv"
univ <- rownames(myTpm)[emean$NP > 30 & maxId == 2]
# *** OPTIONAL: ACT IL21
file <- "Pathway/iregulon_il21.tsv"
univ <- rownames(myTpm)[emean$PP > 30 & maxId == 3]

ireg <- read.delim(file, comment.char = ";", stringsAsFactors = F)
factor <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Transcription.factor[x], split = ",")))
target <- lapply(1:nrow(ireg), function(x) unlist(strsplit(ireg$Target.genes[x], split = ",")))

factor <- lapply(factor, function(x) x[x %in% univ])
idx <- sapply(factor, length) > 0

factor <- factor[idx]
target <- target[idx]

edges <- lapply(1:length(factor), function(x) expand.grid(factor[[x]], target[[x]], stringsAsFactors = F))
edges <- do.call(rbind, edges)
edges <- edges[! duplicated(edges), ]

# VISUALIZATION: IGRAPH
igraph.dt <- graph.data.frame(edges)

igraph.dt$layout <- layout.sphere
igraph.dt$layout <- layout.circle
igraph.dt$layout <- layout.fruchterman.reingold 

V(igraph.dt)$color = rep("chartreuse3", length(V(igraph.dt)$name))
V(igraph.dt)$color[V(igraph.dt)$name %in% unlist(factor)] <- "gold"
V(igraph.dt)$size = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 10 + 2
V(igraph.dt)$label.cex = (degree(igraph.dt) - 1) / (max(degree(igraph.dt)) - 1) * 2 + 0.5
V(igraph.dt)$label.color = rep("dodgerblue3", length(V(igraph.dt)$name))
V(igraph.dt)$label.color[V(igraph.dt)$color == "gold"] = "firebrick1"

plot.igraph(igraph.dt, vertex.frame.color = "white", edge.arrow.size = 0.3)

igraphList <- list()
igraphList$NN <- igraph.dt
igraphList$NP <- igraph.dt
igraphList$PP <- igraph.dt
save(igraphList, file = "data/igraph.rdt")
