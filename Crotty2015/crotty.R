library(affy)
library(oligo)
library(ape)
library(amap)
library(xlsx)
library(genefilter)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(preprocessCore)
library(quantro)
library(qvalue)
library(dplyr)
library(tidyr)
library(reshape)
library(Biobase)

rm(list = ls())
load("~/GitHub/X/summary.rdt")
source("~/GitHub/X/function.R")
setwd("~/GitHub/Il21/Crotty2015")

load("../data/myTpm.rdt")
matboxplot(myTpm, groupFactor = gsub("[12]", "", names(myTpm)))
myTpm[rowMax(as.matrix(myTpm)) > 1e4, ]
myTpm <- myTpm[! rowMax(as.matrix(myTpm)) > 1e4, ]
myTpm <- sweep(myTpm, 2, colSums(myTpm), "/") * 1e6 
myTpm <- log2(myTpm + 1)
myTpm <- myTpm[rowMax(as.matrix(myTpm)) > 5, ]
group <- gsub("[12]", "", names(myTpm))
aov.pval <- apply(myTpm, 1, function (x) summary(aov(x ~ group))[[1]][["group", "Pr(>F)"]])
myTpm <- myTpm[aov.pval < 0.05, ]
myTpm_mean <- data.frame(NN = rowMeans(myTpm[1:2]), NP = rowMeans(myTpm[3:4]), PP = rowMeans(myTpm[5:6]))
colSums(myTpm_mean)

load("../data/myCount.rdt") # note: do not cluster well 
matboxplot(myCount, groupFactor = gsub("[12]", "", names(myCount)))
myCount[rowMax(as.matrix(myCount)) > 1e5, ]
myCount <- myCount[! rowMax(as.matrix(myCount)) > 1e5, ]
myCount <- sweep(myCount, 2, colSums(myCount), "/") * 1e7 
myCount <- myCount[rowMax(as.matrix(myCount)) > 10, ]
group <- gsub("[12]", "", names(myCount))
aov.pval <- apply(myCount, 1, function (x) summary(aov(x ~ group))[[1]][["group", "Pr(>F)"]])
myCount <- myCount[aov.pval < 0.05, ]
myCount <- data.frame(NN = rowMeans(myCount[1:2]), NP = rowMeans(myCount[3:4]), PP = rowMeans(myCount[5:6]))

cel.files <- read.celfiles(list.celfiles("GSE21380_CEL/", full.name = T))
eight <- exprs(rma(cel.files))

file <- "Mouse430_2.na35.annot.csv.gz"
mouse430 <- read.table(file, header = T, comment.char = "#", sep = ",", stringsAsFactors = F)
mouse430 <- mouse430[c("Probe.Set.ID", "Gene.Symbol")]
mouse430$Gene.Symbol <- gsub(" .*", "", mouse430$Gene.Symbol)

eight <- eight[mouse430$Probe.Set.ID, ]
eight <- apply(eight, 2, function(x) tapply(x, mouse430$Gene.Symbol, max))
colnames(eight) <- c("NonTFH_Rep1", "NonTFH_Rep2", "TFH_Rep1", "TFH_Rep2", "GCTFH_Rep1", "GCTFH_Rep2", "Naive_CD4_T")
eight <- eight[-1, ] %>% as.data.frame
matboxplot(eight, groupFactor = gsub("_.*", "", names(eight)))
eight_mean <- data.frame(TH1 = rowMeans(eight[1:2]), TFH = rowMeans(eight[3:4]), GCTFH = rowMeans(eight[5:6]), NV = eight[7])

file <- "GSE67334_TFH.Vs.TH1.DESeq_Analysis_DE_Genes.txt"
three <- read.table(file, header = T, row.names = 1, sep = "\t") 
three_expr <- three[grep("Rep", names(three))] # read count
matboxplot(three_expr, groupFactor = gsub("_.*", "", names(three_expr)))
three_mean <- three[grep("^Mean", names(three))]
names(three_mean) <- gsub("Mean.", "", names(three_mean))

txdb <- keepStandardChromosomes(TxDb.Mmusculus.UCSC.mm10.knownGene)
entrez = AnnotationDbi::select(org.Mm.eg.db, rownames(three), columns=c("ENTREZID"), keytype="SYMBOL")
nrow(entrez[is.na(entrez$ENTREZID), ])
entrez.keep = entrez$ENTREZID[! is.na(entrez$ENTREZID)]
tx <- transcriptsBy(txdb, "gene")
tx <- tx[names(tx) %in% entrez.keep]
tx.len <- mean(width(tx)) # caveat
names(tx.len) <- entrez$SYMBOL[match(names(tx.len), entrez$ENTREZID)]

three_mean <- three_mean[names(tx.len), ]
three_mean <- sweep(three_mean, 1, tx.len, "/")
three_mean <- sweep(three_mean, 2, colSums(three_mean), "/") * 1e6 
three_mean <- log2(three_mean + 1)

three_expr <- three_expr[names(tx.len), ]
three_expr <- sweep(three_expr, 1, tx.len, "/")
three_expr <- sweep(three_expr, 2, colSums(three_expr), "/") * 1e6 
three_expr <- log2(three_expr + 1)

# BATCH CORRECTION
data <- inner_join(add_rownames(myTpm, "gene"), add_rownames(three_expr, "gene"))
data <- inner_join(add_rownames(myTpm_mean, "gene"), add_rownames(three_mean, "gene"))
data <- data.frame(row.names = data$gene, data[-1])
data <- normalize.quantiles(as.matrix(data))
heatmap(cor(data)) # normalize to transcript length is a must do
colSums(data)

names(data) <- gsub("^T.*\\.", "", names(data))

data <- inner_join(add_rownames(myTpm_mean, "gene"), add_rownames(eight, "gene"))
data <- data.frame(row.names = data$gene, data[-1])
heatmap(cor(data)) # two naive samples clustered without norm.

dt.ctr <- data - rowMeans(data)
dt.svd <- svd(dt.ctr) # SVD
barplot(dt.svd$v[, 1]) # PC1
dt2 <- t(apply(dt.ctr, 1, function (x) {lm(x ~ dt.svd$v[, 1])$res})) 
dt3 <- dt2 + rowMeans(data)

heatmap(cor(dt2))
heatmap(cor(dt3))
hc1 <- hcluster(t(dt3), method = "pearson", link = "average")

hc1 <- hcluster(t(data), method = "pearson", link = "average")

pdf("phylo3.pdf", width = 6, height = 5)
plot(as.phylo(hc1), edge.width = 2, font = 2, label.offset = 3e-4)
dev.off()

tfh <- dt3[, c("PP1", "PP2", "TFH_Rep1", "TFH_Rep2")]
colSums(tfh)

group <- factor(c("PP", "PP", "TFH", "TFH"), levels = c("PP", "TFH"))
ttest <- rowttests(tfh, group) # ttests

tfh <- as.data.frame(tfh)
tfh$PP_avg <- rowMeans(tfh[group == "PP"])
tfh$TFH_avg <- rowMeans(tfh[group == "TFH"])
tfh <- cbind(tfh, ttest)

tfh_select <- tfh[tfh$p.value < 0.05 & abs(tfh$dm) > 0.5, ]
gk <- mmGK(rownames(tfh_select))
data.frame(KEGG = gk$KEGG$Term[1:20], BP = gk$GO$BP$Term[1:20], MF = gk$GO$BP$Term[1:20])

gdt <- tfh
gdt$select = "N"
gdt[rownames(tfh_select), "select"] = "Y"
gdt$select[gdt$select == "Y" & gdt$dm < 0] = "TFH"
gdt$select[gdt$select == "Y" & gdt$dm > 0] = "PP"

gdt$select = factor(gdt$select, levels = c("N", "TFH", "PP"))

derry <- read.table("select_derry")$V1 %>% as.vector
gdt$label = rownames(gdt)
gdt$label[! gdt$label %in% derry] = ""

pdf("scatterplot7.pdf", width = 7, height = 6)

ggplot(gdt, aes(x = PP_avg, y = TFH_avg, label = label)) + 
  geom_point(aes(colour = select)) +
  geom_abline(intercept = 0, slope = 1) + 
# geom_text(size = 3) +
  scale_color_manual(values = c("grey30", "dodgerblue3", "firebrick1")) + 
  theme_bw() + theme(legend.title = element_blank())

dev.off()

x = tfh_select[order(tfh_select$dm), ]
write.xlsx(x, file = "tfh.xlsx", sheetName = "TFH_2")

y = dt2[c("PP", "TFH")]
y = y[abs(y$PP - y$TFH) > 0.1, ]
plot(y)
abline(0, 1)

(select = unlist(read.table(file = "select")) %>% as.vector)

setdiff(select, rownames(myTpm))

gdt = myTpm[intersect(select, rownames(myTpm)), ]
gdt = reshape2::melt(add_rownames(gdt, "gene"))

pdf("select.pdf", width = 9, height = 5)
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
        axis.text.y = element_text(size = 10, face = 2))
dev.off()

#
tfh = read.xlsx("tfh.xlsx", sheetName = "TFH_2", row.names = 1, stringsAsFactors = F)
tfh_up = rownames(tfh)[tfh$dm < 0] # High in TFH
tfh_down = rownames(tfh)[tfh$dm > 0] # High in PP

tfh_up_go = mmGK(tfh_up)
tfh_down_go = mmGK(tfh_down)

write.xlsx(tfh_up_go$KEGG, file = "TFH_KEGG.xlsx", sheetName = "up", append = T)
write.xlsx(tfh_down_go$KEGG, file = "TFH_KEGG.xlsx", sheetName = "down", append = T)

write.xlsx(tfh_up_go$GO$BP, file = "TFH_GO.xlsx", sheetName = "up_BP", append = T)
write.xlsx(tfh_up_go$GO$MF, file = "TFH_GO.xlsx", sheetName = "up_MF", append = T)
write.xlsx(tfh_up_go$GO$CC, file = "TFH_GO.xlsx", sheetName = "up_CC", append = T)

write.xlsx(tfh_down_go$GO$BP, file = "TFH_GO.xlsx", sheetName = "down_BP", append = T)
write.xlsx(tfh_down_go$GO$MF, file = "TFH_GO.xlsx", sheetName = "down_MF", append = T)
write.xlsx(tfh_down_go$GO$CC, file = "TFH_GO.xlsx", sheetName = "down_CC", append = T)
