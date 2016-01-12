library(xlsx)
library(dplyr)
library(VennDiagram)
library(Biobase)

rm(list = ls())
setwd("~/GitHub/Il21")

options(stringsAsFactors = F)

# p.adjust

tfh = read.xlsx("./Crotty2015/tfh.xlsx", sheetName = "TFH_2", row.names = 1)

tfh_up = rownames(tfh)[tfh$dm < 0] # High in TFH
tfh_down = rownames(tfh)[tfh$dm > 0] # High in PP

filter(tfh, p.value < 0.05 / 1e3)

rownames(tfh)[tfh$p.value < 0.05 / 1e2]
rownames(tfh)[tfh$p.value < 0.05 / 1e3]

# Crotty's

file <- "./Crotty2015/GSE67334_TFH.Vs.TH1.DESeq_Analysis_DE_Genes.txt"
three <- read.table(file, header = T, row.names = 1, sep = "\t") 

three$gene <- rownames(three)

three_exp <- three[grep("Rep", names(three))] # read count
table(rowMax(as.matrix(three_exp)) > 20)

three = filter(three, rowMax(as.matrix(three_exp)) > 20)

table(three$padj < 0.05)

tfh = filter(three, padj < 0.05 & log2FoldChange > 0.5)
th1 = filter(three, padj < 0.05 & log2FoldChange < -0.5)

# IL21's

nn = read.xlsx("./Results/new.xlsx", sheetName = "Naive", row.names = 1)
np = read.xlsx("./Results/new.xlsx", sheetName = "ACT", row.names = 1)
pp = read.xlsx("./Results/new.xlsx", sheetName = "ACT-IL21", row.names = 1)

col.a = c("red", "white", "white")
col.b = c("white", "white", "white", "red")

venn.diagram(list(PP = rownames(pp), TFH = tfh$gene, TH1 = th1$gene), fill = col.a, file = "./venn.png", imagetype = "png")
venn.diagram(list(NP = rownames(np), TFH = tfh$gene, TH1 = th1$gene), fill = col.a, file = "./venn.png", imagetype = "png")
venn.diagram(list(NN = rownames(nn), TFH = tfh$gene, TH1 = th1$gene), fill = col.a, file = "./venn.png", imagetype = "png")

venn.diagram(list(NN = rownames(nn), NP = rownames(np), PP = rownames(pp), TFH = tfh$gene), fill = col1, 
             file = "./venn.png", imagetype = "png")

