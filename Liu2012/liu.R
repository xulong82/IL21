# Copyright: Xulong Wang (xulong.wang@jax.org)

library(ape)
library(amap)
library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")
load("./data/myTpm.rdt")
load("./data/ebseq.rdt")
geneId <- rownames(sig0.95)
load("./data/liu.rdt")

rnaseq <- myTpm[geneId, ]
rnaseq <- log2(rnaseq + 1)
rnaseq <- rnaseq[apply(rnaseq, 1, max) > 5, ]

idx <- intersect(rownames(rnaseq), rownames(liu))  
data1 <- cbind(rnaseq[idx, ], liu[idx, ])

# --- BATCH CORRECTION ---
dt.ctr <- data1 - apply(data1, 1, mean)
dt.svd <- svd(dt.ctr)  # SVD
barplot(dt.svd$v[, 1])

dt2 <- t(apply(dt.ctr, 1, function (x) {lm(x ~ dt.svd$v[, 1])$res}))  # PC1
dt2 <- as.data.frame(dt2)

hc1 <- hcluster(t(dt2), method = "pearson", link = "average")
plot(as.phylo(hc1), edge.width = 2, font = 2, label.offset = .03)

col.manual <- c("darkblue", "darkorchid2", "darkgreen")
clusts = cutree(hc1, 3)

pdf("Results/liuPhylo1.pdf", fonts = "Helvetica")
plot(as.phylo(hc1), edge.width = 3, font = 2, label.offset = .03)
# plot(as.phylo(hc1), tip.color = col.manual[clusts], edge.width = 2, font = 2, label.offset = .03)
dev.off()

hc2 <- hcluster(dt2, method = "correlation", link = "centroid")
tile1 <- dt2[hc2$order, hc1$order]
tile2 <- t(apply(tile1, 1, scale))
colnames(tile2) <- colnames(tile1)
tile2 <- data.frame(value = c(tile2), gene = rep(rownames(tile1), 6), 
                    sample = rep(colnames(tile1), each = nrow(tile1)))
tile2$sample <- factor(tile2$sample, levels = colnames(tile1))
tile2$gene <- factor(tile2$gene, levels = rownames(tile1))

pdf("Results/heatmapLiu.pdf", width = 15, height = 5)
  ggplot(tile2, aes(x = sample, y = gene, fill = value)) + 
  geom_tile() + guides(alpha = F) + coord_flip() +
  scale_fill_gradient2(low = "blue", high = "red") +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(size = 2, angle = 90),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.title = element_blank())
dev.off()

# --- PROFILIC GENE LIST ---
# aov.pval <- NULL
# group <- c("PP", "PP", "NN", "NN", "NP", "NP")
# for (i in 1:nrow(liu)) {
#   if (i %% 1e3 == 0) cat(i, "\n")
#   dat1 <- data.frame(tpm = as.matrix(liu)[i, ], group)
#   aov1 = aov(tpm ~ group, data = dat1)
#   aov.pval <- c(aov.pval, summary(aov1)[[1]][["Pr(>F)"]][1])
# }
# 
# expr.aov <- liu[aov.pval < 0.05, ]  # Select genes
# expr.mean <- data.frame("NN" = rowMeans(expr.aov[3:4]), "NP" = rowMeans(expr.aov[5:6]), "PP" = rowMeans(expr.aov[1:2]))
# idx.sample = apply(expr.mean, 1, which.max)
# 
# profile3 <- list()
# group <- c("1", "1", "2", "2")
# for (sample in c("NN", "NP", "PP")) {
#   cat(sample, "\n")
#   
#   if (sample == "NN") {
#     expr0 <- expr.aov[idx.sample == 1, ]
#     expr01 <- expr0[, 1:4]
#     expr02 <- expr0[, 3:6]
#   } else if (sample == "NP") {
#     expr0 <- expr.aov[idx.sample == 2, ]
#     expr01 <- expr0[, c(1:2, 5:6)]
#     expr02 <- expr0[, 3:6]
#   } else {
#     expr0 <- expr.aov[idx.sample == 3, ]
#     expr01 <- expr0[, 1:4]
#     expr02 <- expr0[, c(1:2, 5:6)]
#   }
#   
#   lm1 <- apply(expr01, 1, function(x) lm(x ~ group))
#   lm2 <- apply(expr02, 1, function(x) lm(x ~ group))
#   idx1 <- lapply(lm1, function(x) summary(x)$coefficients[2, "Pr(>|t|)"] < 0.05)
#   idx2 <- lapply(lm2, function(x) summary(x)$coefficients[2, "Pr(>|t|)"] < 0.05)
#   # idx3 <- lapply(lm3, function(x) summary(x)$coefficients[2, "Pr(>|t|)"] > 0.05)
#   # idx.gene <- unlist(idx1) & unlist(idx2) & unlist(idx3)
#   idx.gene <- unlist(idx1) & unlist(idx2)
#   
#   profile[[sample]] = rownames(expr0)[idx.gene]
# }
