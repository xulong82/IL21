# Copyright: Xulong Wang (xulong.wang@jax.org)

library(ape)
library(amap)
library(ggplot2)
library(lattice)
library(ggdendro)
library(grid)
library(xtable)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")
load("./data/myTpm.rdt")

expr <- myTpm[apply(myTpm, 1, function(x) max(x) > 50), ]
group <- gsub("[12]", "", colnames(expr))
pval <- apply(expr, 1, function (x) min(summary(aov(x ~ group))[[1]][["Pr(>F)"]], na.rm = T))
expr <- expr[pval < 0.05, ]
expr <- log2(expr + 1)

emean <- data.frame("NN"=rowMeans(expr[1:2]), "NP"=rowMeans(expr[3:4]), "PP"=rowMeans(expr[5:6]))
idx.sample = apply(emean, 1, which.max)

profile1 <- profile2 <- list()
group <- c("1", "1", "2", "2")
for (sample in c("NN", "NP", "PP")) {
  cat(sample, "\n")
  
  if (sample == "NN") {
    expr0 <- expr[idx.sample == 1, ]
    expr01 <- expr0[, 1:4]
    expr02 <- expr0[, c(1:2, 5:6)]
    expr03 <- expr0[, 3:6]
  } else if (sample == "NP") {
    expr0 <- expr[idx.sample == 2, ]
    expr01 <- expr0[, 1:4]
    expr03 <- expr0[, c(1:2, 5:6)]
    expr02 <- expr0[, 3:6]
  } else {
    expr0 <- expr[idx.sample == 3, ]
    expr03 <- expr0[, 1:4]
    expr01 <- expr0[, c(1:2, 5:6)]
    expr02 <- expr0[, 3:6]
  }
  
  lm1 <- apply(expr01, 1, function(x) lm(x ~ group))
  lm2 <- apply(expr02, 1, function(x) lm(x ~ group))
  lm3 <- apply(expr03, 1, function(x) lm(x ~ group))
  
  idx1 <- lapply(lm1, function(x) summary(x)$coefficients[2, "Pr(>|t|)"] < 0.05)
  idx2 <- lapply(lm2, function(x) summary(x)$coefficients[2, "Pr(>|t|)"] < 0.05)
  idx3 <- lapply(lm3, function(x) summary(x)$coefficients[2, "Pr(>|t|)"] > 0.05)
  
  idx.gene1 <- unlist(idx1) & unlist(idx2)
  idx.gene2 <- unlist(idx1) & unlist(idx2) & unlist(idx3)
  profile1[[sample]] = rownames(expr0)[idx.gene1]
  profile2[[sample]] = rownames(expr0)[idx.gene2]
}

save(profile1, profile2, file = "./data/profile12.rdt")
geneId <- unlist(profile1)

# --- SCATTERPLOT ON PROFILIC GENES ---
scatter <- myTpm[apply(myTpm, 1, function(x) max(x) > 5), ]
scatter <- log2(scatter + 1)
scatter <- data.frame("NN" = rowMeans(scatter[1:2]), "NP" = rowMeans(scatter[3:4]), "PP" = rowMeans(scatter[5:6]))

scatter$NP.NN <- scatter$NP - scatter$NN
scatter$PP.NN <- scatter$PP - scatter$NN
scatter$group <- rep("None", nrow(scatter))
scatter$group[rownames(scatter) %in% profile1[[1]]] <- "N"
scatter$group[rownames(scatter) %in% profile1[[2]]] <- "ACT"
scatter$group[rownames(scatter) %in% profile1[[3]]] <- "ACT IL21"

scatter1 <- scatter[scatter$group == "None", ]
scatter2 <- scatter[! scatter$group == "None", ]
scatter2$group <- factor(scatter2$group, levels = c("N", "ACT", "ACT IL21"))

hlId <- c("Foxo1", "Il4ra", "Sell", "Bcl6", "Btla", "Cd28", "Cxcr5", "Il21", "Maf", "Foxp3", "Il2ra", "Ccl5")
myLabel <- scatter[hlId, ]

pdf("Results/scatterplot.pdf", width = 5, height = 4)
ggplot() + 
  geom_point(data = scatter1, aes(x = NP.NN, y = PP.NN), color = "grey70", shape = 4, size = 1.5) +
  geom_point(data = scatter2, aes(x = NP.NN, y = PP.NN, color = group, shape = group), size = 1.5) +
  theme_bw() + xlab("ACT vs N") + ylab("ACT IL21 vs N") +
  scale_color_manual(values = c("darkorchid2", "chartreuse3", "dodgerblue3")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12, face = "bold"),
        legend.key = element_blank()) +
  geom_point(data = myLabel, aes(x = NP.NN, y = PP.NN), shape = 2, size = 3, color = "red") 
#   geom_text(data = myLabel, aes(x = NP.NN, y = PP.NN, label = rownames(myLabel)), 
#             color = "black", size = 3, vjust = -1.0)
dev.off()

# --- PHYLO GRAPH ON SAMPLE HC ---
dt <- expr[geneId, ]
colnames(dt) <- c("N S1", "N S2", "ACT S1", "ACT S2", "ACT IL21 S1", "ACT IL21 S2")
hc1 <- hcluster(t(dt), method = "pearson", link = "average")
hc2 <- hcluster(dt, method = "correlation", link = "centroid")
hc1$labels <- gsub(" S[12]", "", hc1$labels)

pdf("Results/phylo_profile1.pdf", fonts = "Helvetica", width = 5, height = 4)
par(mar = c(5, 2, 4, 2))
col.manual <- c("darkorchid2", "chartreuse3", "dodgerblue3", "grey20")
clusts = cutree(hc1, 3)
plot(as.phylo(hc1), tip.color = col.manual[clusts], 
     edge.color = col.manual[c(4, 1, 1, 4, 4, 2, 2, 4, 3, 3)],
     edge.width = 5, srt = 90, adj = 0.8, cex = 1, font = 2, 
     direction = "downward", label.offset = -5e-4)
dev.off()

# --- HEATMAP ON PROLIFIC GENES ---
tile1 <- dt[hc2$order, hc1$order]
tile2 <- t(apply(tile1, 1, scale))
colnames(tile2) <- colnames(tile1)
tile2 <- data.frame(value = c(tile2), gene = rep(rownames(tile1), 6), 
                    sample = rep(colnames(tile1), each = nrow(tile1)))
tile2$sample <- factor(tile2$sample, levels = colnames(tile1))
tile2$gene <- factor(tile2$gene, levels = rownames(tile1))

pdf("Results/heatmap.pdf", width = 15, height = 3)
ggplot(tile2, aes(x = sample, y = gene, fill = value)) + 
  geom_tile() + guides(alpha = F) + coord_flip() +
  scale_fill_gradient2(low = "blue", high = "red") +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.title = element_blank())
dev.off()

# --- TABLE ---
scatter <- log2(scatter + 1)
scatter <- data.frame("NN" = rowMeans(scatter[1:2]), "NP" = rowMeans(scatter[3:4]), "PP" = rowMeans(scatter[5:6]))

scatter$NP.NN <- scatter$NP - scatter$NN
scatter$PP.NN <- scatter$PP - scatter$NN
scatter$group <- rep("None", nrow(scatter))
scatter$group[rownames(scatter) %in% profile1[[1]]] <- "N"
scatter$group[rownames(scatter) %in% profile1[[2]]] <- "ACT"
scatter$group[rownames(scatter) %in% profile1[[3]]] <- "ACT IL21"

group <- c(rep("N", length(profile1[[1]])), rep("ACT", length(profile1[[2]])), rep("ACT IL21", length(profile1[[3]])))
table <- cbind(group, geneId, myTpm[geneId, ], ANOVA = pval[geneId]) 
table$"FC(NP/NN)" <- rowMeans(table[, c("NP1", "NP2")]) / rowMeans(table[, c("NN1", "NN2")])
table$"FC(PP/NN)" <- rowMeans(table[, c("PP1", "PP2")]) / rowMeans(table[, c("NN1", "NN2")])
rownames(table) <- NULL
write.xlsx(table, file = "Results/signatures.xlsx", sheetName = "genes", row.names = F, append = T)
