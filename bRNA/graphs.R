library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/Il21")
load("./data/myTpm.rdt")

derry <- read.delim("./data/derry.txt", stringsAsFactors = F, header = F)
geneId <- gsub(" ", "", derry$V1)
geneId <- gsub("\\(GITR\\)", "", geneId)

chosen <- myTpm[geneId, ]
group <- gsub("[12]", "", colnames(chosen))
pval <- apply(chosen, 1, function (x) min(summary(aov(x ~ group))[[1]][["Pr(>F)"]], na.rm = T))
table <- cbind(geneId, chosen, ANOVA = pval[geneId]) 
table$NN <- rowMeans(table[, c("NN1", "NN2")] + 1)
table$NP <- rowMeans(table[, c("NP1", "NP2")] + 1)
table$PP <- rowMeans(table[, c("PP1", "PP2")] + 1)
table$"FCNP" <- table$NP / table$NN
table$"FCPP" <- table$PP / table$NN
write.xlsx(table, file = "Results/chosen.xlsx", sheetName = "genes", row.names = F, append = T)

geneId.tfh <- geneId[1:34]
gdt <- table[geneId.tfh, ]

pdf("./Figure/TFH.pdf", width = 10, height = 5)
# pdf("./Figure/THC.pdf", width = 5, height = 7)

ggplot() + 
  geom_point(data = gdt, aes(x = geneId, y = log2(FCNP), size = NP), color = "chartreuse3") + 
  geom_point(data = gdt, aes(x = geneId, y = log2(FCPP), size = PP), color = "dodgerblue3") +
  theme_bw() + xlab("") + ylab("Log2 fold change") + # coord_flip() +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "grey30"),
#       axis.text.x = element_text(size = 10, face = 2, vjust = 0.5),
        axis.text.x = element_text(size = 10, face = 2, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10, face = 2),
        axis.title = element_text(size = 12, face = 2),
        legend.position = "none")

dev.off()

gdt <- myTpm[geneId.tfh, ]
geneId.thc <- geneId[35:61]
gdt <- myTpm[geneId.thc, ]

gdt <- data.frame("NN" = rowMeans(gdt[1:2]), "NP" = rowMeans(gdt[3:4]), "PP" = rowMeans(gdt[5:6]))
gdt <- log2(gdt + 1)
gdt <- data.frame(value = c(as.matrix(gdt)), 
                  gene = factor(rep(rownames(gdt), time = ncol(gdt)), levels = rownames(gdt)),
                  sample = factor(rep(colnames(gdt), each = nrow(gdt)), levels = colnames(gdt)))

pdf("./Figure/TFH.pdf", width = 5, height = 8)
pdf("./Figure/TFH.pdf", width = 8, height = 5)

ggplot() + 
  geom_line(data = gdt, aes(x = gene, y = value), color = "grey80", size = 2) +
  geom_point(data = gdt, aes(x = gene, y = value, color = sample), size = 3) +
  theme_bw() + xlab("") + ylab("") + # coord_flip() +
  scale_color_manual(values = c("darkorchid2", "chartreuse3", "dodgerblue3"),
                     breaks = c("NN", "NP", "PP"),
                     labels = c("N", "ACT", "ACT IL21")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "grey30"),
#       axis.text.x = element_text(size = 10, face = 2, vjust = 0.5),
        axis.text.x = element_text(size = 10, face = 2, angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 10, face = 2),
        legend.position = "none")

dev.off()

bar.dt <- data.frame("NN" = rowMeans(expr[1:2]), "NP" = rowMeans(expr[3:4]), "PP" = rowMeans(expr[5:6]))
bar.dt <- bar.dt[geneId, 1:3]
colnames(bar.dt) <- c("Naive", "Act", "Act Il21")
bar.dt <- data.frame(value = c(as.matrix(bar.dt)), 
                     gene = factor(rep(rownames(bar.dt), time = ncol(bar.dt)), levels = geneId),
                     sample = factor(rep(colnames(bar.dt), each = nrow(bar.dt)), levels = colnames(bar.dt)))

pdf("./Public/barplot.pdf", width = 12, height = 7)

# ggplot(bar.dt[bar.dt$gene %in% geneId.tfh, ], aes(x = gene, y = value)) +
ggplot(bar.dt, aes(x = gene, y = value)) +
  # geom_bar(aes(fill = sample), stat = "identity", position = "dodge", width = 0.65) +
  geom_bar(aes(fill = sample), stat = "identity", width = 0.65) +
  scale_fill_manual(values = c("darkorchid2", "chartreuse3", "dodgerblue3")) +
  theme_bw() + xlab("") + ylab("") + # coord_flip() +
  scale_color_manual(values = c("darkorchid2", "chartreuse3", "dodgerblue3")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = 'grey30'),
        axis.text.x = element_text(size = 10, angle = -90, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "none", 
        legend.direction = "horizontal",
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.key = element_blank()) 

dev.off()

pdf("./Results/myplot3.pdf")

par(mar = c(5, 4, 4, 2))
ggplot(data.tfh4, aes(x = as.numeric(sample), y = Log2.Ratio, group = gene, colour = gene)) + 
  geom_line(size = 1.5) + geom_point(size = 5) +
  theme_bw() +
  xlab("Samples") + ylab("Log2.Ratio") +
  scale_x_continuous(breaks=c(1,2,3), labels = c("-/-", "-/+", "+/+")) +
  scale_colour_brewer(palette="Set1") +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(vjust = -0.5)) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())

dev.off()

pdf("./Results/myplot2.pdf")

par(mar = c(5, 4, 4, 2))
ggplot(data.tfh2, aes(x = gene, y = Log2.Ratio, fill = sample)) + geom_boxplot() +
  theme_bw() +
  theme(panel.border = element_rect(size = 2, color = "black")) +
  theme(axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 

dev.off()

# Jan 25, 2016, heatmap on Derry's. HEATMAP ON CHOSEN GENES

gdt <- data.frame("N" = rowMeans(myTpm[1:2]), "ACT" = rowMeans(myTpm[3:4]), "IL21-ACT" = rowMeans(myTpm[5:6]))

name = c("Bcl6", "E2f2", "Id3", "Fosb", "Tox", "Tox2", "Egr2", "Maf", "Nfatc1", "Pou2af1")
name = c("Tcf7", "Lef1", "Prdm1")
name = c("Il21", "Tnsf8", "Tgfb3", "Angpotl2")

dt = gdt[name, ]

# dt <- log2(dt + 1)
# dt <- t(apply(dt, 1, scale))

tile <- data.frame(value = c(as.matrix(dt)), gene = rep(rownames(dt), 3), sample = rep(names(gdt), each = nrow(dt)))
tile$sample = factor(tile$sample, levels = names(gdt))

tile$value <- log2(tile$value + 1)
tile$value <- scale(tile$value)

pdf("./Results/heatmap_tf.pdf", width = 6, height = 2)

ggplot(tile, aes(x = sample, y = gene, fill = value)) + 
  geom_tile() + guides(alpha = F) + coord_flip() +
  scale_fill_gradient2(low = "blue", high = "red") +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(size = 12, face = "bold", angle = 90),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_blank())

dev.off()

# Jan 18, 2016, graph request from Liz and Sandy

sandy = c("Bcl6", "Maf", "Nfatc1", "Prdm1", "Ascl2", "Irf4", "Il21", "Il10", "Ifng", "Cxcr4", "Cxcr5", "Cd28", "Cd200", "Il4", "Il5", "Il21r", "Foxp3", "Il12rb1", "Il12rb2", "Slamf6")
dt = myTpm[sandy, ]

names(dt) <- c("N S1", "N S2", "ACT S1", "ACT S2", "ACT IL21 S1", "ACT IL21 S2")
hc1 <- hcluster(t(dt))
hc2 <- hcluster(dt)
hc2 <- hcluster(dt, method = "correlation", link = "centroid")

# --- HEATMAP ON PROLIFIC GENES ---

tile1 <- dt[hc2$order, hc1$order]
gene = rownames(tile1)
sample = names(tile1)

# tile1 <- log2(tile1 + 1)
# tile2 <- apply(tile1, 2, function(x) scale(x, center = T, scale = T))
# tile2 <- apply(tile1, 2, function(x) scale(x, center = F, scale = T))

# tile2 <- t(apply(tile1, 1, scale))
# colnames(tile2) <- colnames(tile1)

tile2 <- data.frame(value = c(as.matrix(tile1)), gene = rep(gene, 6), sample = rep(sample, each = nrow(tile1)))

tile2$sample <- factor(tile2$sample, levels = sample)
tile2$gene <- factor(tile2$gene, levels = gene)

tile2$value <- log2(tile2$value + 1)
tile2$value <- scale(tile2$value)

pdf("./heatmap.pdf", width = 12, height = 3)
ggplot(tile2, aes(x = sample, y = gene, fill = value)) + 
  geom_tile() + guides(alpha = F) + coord_flip() +
  scale_fill_gradient2(low = "blue", high = "red") +
  theme_bw() + xlab("") + ylab("") + 
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.title = element_blank())
dev.off()
