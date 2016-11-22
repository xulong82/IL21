# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: TCR transcript analysis

rm(list = ls())
library(ggplot2)
library(pheatmap)

setwd("~/Dropbox/GitHub/Il21")
load("./data/myTpm.rdt")

myTpm <- data.frame(row.names = rownames(myTpm), 
  NN = rowMeans(myTpm[1:2]), NP = rowMeans(myTpm[3:4]), PP = rowMeans(myTpm[5:6]))

tx2type <- read.delim("./data/tx2type.txt", stringsAsFactors = F)
tr.v <- tx2type[grep("Tr[ab]v", tx2type$Associated.Gene.Name), ]
tr.j <- tx2type[grep("Tr[ab]j", tx2type$Associated.Gene.Name), ]
tr.d <- tx2type[grep("Trbd", tx2type$Associated.Gene.Name), ]

tx.tr.v <- myTpm[unique(tr.v$Associated.Gene.Name), ]
tx.tr.j <- myTpm[unique(tr.j$Associated.Gene.Name), ]
tx.tr.d <- myTpm[unique(tr.d$Associated.Gene.Name), ]

tx.tr.v <- tx.tr.v[apply(tx.tr.v, 1, max) > 5, ]
tx.tr.v <- tx.tr.v[apply(tx.tr.v, 1, sd) > 5, ]
colnames(tx.tr.v) <- c("N", "ACT", "ACT IL21")

gdt <- data.frame(TPM = c(as.matrix(tx.tr.v)),
  sample = rep(c("N", "ACT", "ACT IL21"), each = nrow(tx.tr.v)),
  gene = rep(rownames(tx.tr.v), ncol(tx.tr.v)))
gdt$sample = factor(gdt$sample, levels = c("N", "ACT", "ACT IL21"))

pdf(file = "Results/TCR.pdf", family = "Arial", width = 12, height = 5)
ggplot(gdt, aes(x = gene, y = TPM, fill = gene)) + 
  geom_bar(stat = "identity") + facet_grid(sample ~ .) +
  theme_bw() + xlab("") +
  theme(panel.border = element_rect(size = 1, color = "grey30")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold", vjust = 1.5),
        legend.position = "none")
dev.off()

pdf("Figure/TCR_pheatmap.pdf", width = 2.5, height = 3)
pheatmap(cor(tx.tr.v, method = "pearson"), 
  display_number = T, treeheight_row = 0, border_color = "white", legend = F)
dev.off()

# --- STATISTICS ---
load("./data/myTpm.rdt")
myTpm <- myTpm[rownames(tx.tr.v), ]
group <- gsub("[12]", "", colnames(myTpm))
pval <- apply(myTpm, 1, function (x) min(summary(aov(x ~ group))[[1]][["Pr(>F)"]], na.rm = T))
pval[pval < 0.05]
              