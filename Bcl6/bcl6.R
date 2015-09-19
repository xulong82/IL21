# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: Bcl6

library(ggplot2)
library(circlize)

rm(list = ls())
cat("--- Bcl6 binding sites --- \n")
bcl6.bed <- read.delim("~/Dropbox/Lupus/Bcl6/MA0463.1.bed", header = F, stringsAsFactors = F)[, 1:3]
bcl6.pm.bed <- read.delim("~/Dropbox/Lupus/Bcl6/bcl6.promoter.bed", header = F, stringsAsFactors = F)[, 1:3]
bcl6.bed$type <- "full"
bcl6.pm.bed$type <- "promoter"
bcl6.bed <- rbind(bcl6.bed, bcl6.pm.bed)
bcl6.bed$chr <- with(bcl6.bed, factor(V1, levels = paste("chr",c(1:19,"X","Y"),sep=""), ordered = TRUE))
bcl6.bed$pos <- (bcl6.bed$V2 + bcl6.bed$V3) / 2 * 1e-6

pdf("~/Dropbox/Lupus/Bcl6/stacked.pdf", width = 12)
ggplot(bcl6.bed, aes(x = pos, y = 1)) +
  geom_point(color = "black", size = .1) +
  geom_vline(aes(xintercept = pos, colour = type)) +
  facet_grid(chr ~ .) +
  theme_bw() +
  xlab("Genome position (Mb)") + ylab("") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_blank())
dev.off()

pdf("~/Dropbox/Lupus/Bcl6/circlize.pdf")
circos.initializeWithIdeogram(species = "mm10", plotType = NULL)
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function (x, y) {
  chr = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.rect(xlim[1], 0, xlim[2], 0.2, col = rgb(runif(1), runif(1), runif(1)))
  circos.text(mean(xlim), 0.9, chr, cex = 0.5, facing = "clockwise", niceFacing = TRUE)
}, bg.border = NA)
circos.genomicDensity(bcl6.bed)
# circos.genomicRainfall(bcl6.pm.bed)
})  
dev.off()

cat("--- Genes with Bcl6 peak in the promoter regions --- \n")
load("~/Dropbox/Lupus/R/data2.rdt")
bcl6.pm.bed <- read.delim("~/Dropbox/Lupus/Bcl6/bcl6.promoter.bed", header = F, stringsAsFactors = F)
bcl6.pm.dt <- data2[unique(bcl6.pm.bed$V4), ] 
intersect(unique(bcl6.pm.bed$V4), rownames(data2.sig))

bcl6.pm.dt <- data.frame(id <- rep(rownames(bcl6.pm.dt), ncol(bcl6.pm.dt)),
                       tpm <- c(as.matrix(bcl6.pm.dt)),
                       cond <- rep(colnames(bcl6.pm.dt), each = nrow(bcl6.pm.dt)))

pdf("~/Dropbox/Lupus/Bcl6/bcl6genes.pdf", width = 10)
ggplot(bcl6.pm.dt, aes(x = id, y = tpm, fill = cond)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  xlab("") + ylab("TPM") +
  theme(panel.border = element_rect(size = 1, color = "black")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_discrete(breaks=c("vfp.neg.icos.neg", "vfp.neg.icos.pos", "vfp.pos.icos.pos"),
                      labels=c("-/-", "-/+", "+/+")) +
  theme(legend.position = "top", legend.direction = "horizontal", 
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank())
dev.off()
