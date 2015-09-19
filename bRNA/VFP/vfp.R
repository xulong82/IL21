# Copyright: Xulong Wang (xulong.wang@jax.org)
# Purpose: VFP transcript analysis

rm(list = ls())
library(ggplot2)

setwd("~/Dropbox/GitHub/Lupus")
rsem.tx <- read.delim("./VFP/rsem_genes.txt")
rsem.tx <- data.frame(row.names = rsem.tx$transcript_id.s., rsem.tx[, 2:7])
head(rsem.tx)

load("./data/myTpm.rdt")

vfp.tx <- rsem.tx["YFP_Algosome", ]
colnames(vfp.tx) <- colnames(myTpm)

ctrls <- c("Cd44", "Icos", "Il21")
tpm <- rbind(myTpm[ctrls, ], vfp.tx)

# il21.tx <- c("ENSMUST00000161015", "ENSMUST00000029273")
# il21.tx <- rsem.tx[il21.tx, ]
# il21.tx <- colSums(il21.tx)
# tpm <- rbind(il21.tx, vfp.tx)[, grep("TPM", colnames(rsem.tx))]

gdt <- data.frame(TPM = c(as.matrix(tpm)), 
                  Gene = rep(c(ctrls, "VFP"), 6),  
                  Sample = c(rep("N", 8), rep("ACT", 8), rep("ACT IL21", 8)))
gdt$Sample <- factor(gdt$Sample, levels = c("N", "ACT", "ACT IL21"))

pdf("Results/vfp.pdf", width = 6, height = 5)
ggplot(gdt, aes(x = Gene, y = TPM)) +
  geom_bar(aes(fill = Sample), stat = "identity", position = "dodge", width = 0.65) +
  scale_fill_manual(values = c("darkorchid2", "chartreuse3", "dodgerblue3")) +
  theme_bw() + xlab("") + ylab("") + # coord_flip() +
  scale_color_manual(values = c("darkorchid2", "chartreuse3", "dodgerblue3")) +
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 1, color = "black"),
        axis.text.x = element_text(size = 10, vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
dev.off()

pdf("~/Dropbox/Lupus/vfp.pdf")
ggplot(gdt, aes(x = Sample, y = TPM, colour = Gene, shape = Gene)) + 
  geom_point(size = 5) + theme_bw() + xlab("") +
  scale_colour_manual(values = c("chartreuse3", "firebrick1")) +
  theme(panel.border = element_rect(size = .5, color = "grey30")) +
  theme(axis.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold", vjust = 1.5),
        legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

