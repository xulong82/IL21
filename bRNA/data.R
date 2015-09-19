# Copyright: Xulong Wang (xulong.wang@jax.org)

library(mygene)
library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

#--- Roopenian Lab RnaSeq TPM
myTpm <- read.delim("./RSEM/genes.TPM", sep = "", header = F, stringsAsFactors = F)
myTpm <- aggregate(. ~ V1, data = myTpm, sum)
myTpm <- data.frame(row.names = myTpm$V1, myTpm[, 2:ncol(myTpm)])
colnames(myTpm) <- c("NN1", "NN2", "NP1", "NP2", "PP1", "PP2")
save(myTpm, file = "./data/myTpm.rdt")

myCount <- read.delim("./RSEM/genes.count", sep = "", header = F, stringsAsFactors = F)
myCount <- aggregate(. ~ V1, data = myCount, sum)
myCount <- data.frame(row.names = myCount$V1, myCount[, 2:ncol(myCount)])
colnames(myCount) <- c("NN1", "NN2", "NP1", "NP2", "PP1", "PP2")
save(myCount, file = "./data/myCount.rdt")

# --- Gene summary
myTpm <- myTpm[apply(myTpm, 1, function (x) max(x) > 5), ]
geneId <- rownames(myTpm)
mus2hg <- read.delim("~/Dropbox/X/hg2mus.map", header = F, stringsAsFactors = F)
table <- queryMany(geneId, scopes="symbol", species="mouse", fields = c("name", "summary"))
table <- table[!duplicated(table$query), ]
table$hg <- mus2hg$V1[match(geneId, mus2hg$V4)]
query.hg <- queryMany(table$hg, scopes="symbol", species="human", fields = c("name", "summary"))
query.hg <- query.hg[!duplicated(query.hg$query), ]
summary_hg <- query.hg$summary[match(table$hg, query.hg$query)]
table$summary <- paste("Mus:", table$summary, "Hg:", summary_hg)
table <- table[c("query", "name", "hg", "summary")]
rownames(table) <- NULL
table <- as.data.frame(table)
save(table, file = "./data/gene_summary.rdt")
# table <- queryMany(geneId, scopes="symbol", species="mouse", fields = c("name", "summary", "go", "kegg"))
# table <- annotation[c("query", "name", "summary")]
# table$go.MF <- sapply(annotation$go.MF, function(x) paste(x$term, collapse = ","))
# table$go.BP <- sapply(annotation$go.BP, function(x) paste(x$term, collapse = ","))
# table$go.CC <- sapply(annotation$go.CC, function(x) paste(x$term, collapse = ","))
# --- Shiny
save(myTpm, table, file = "./Shiny/gene_expression_and_summary_10922.rdt")

# --- ImmGen dataset
setwd("~/Dropbox/GitHub/Lupus/Immgen")
immgen <- read.csv("ImmgenTotal.csv.gz", header = T)  # 214 samples
# immgen <- read.csv("gdTCells.csv", header = T)  
# immgen <- read.csv("abTCells.csv", header = T)  
# immgen <- read.csv("BCells.csv", header = T)  
# immgen <- read.csv("DendriticCells.csv", header = T)  
# immgen <- read.csv("Macrophages.csv", header = T)  
# immgen <- read.csv("MonoCytes.csv", header = T)  
# immgen <- read.csv("NKCells.csv", header = T)  
# immgen <- read.csv("Neutrophils.csv", header = T)  
# immgen <- read.csv("StemCells.csv", header = T)  
# immgen <- read.csv("StromalCells.csv", header = T)  
immgen <- immgen[-c(1, 3)]  # Drop the Description column 
immgen <- aggregate(. ~ GeneSymbol, data = immgen, max)  # Aggregate multiple probesets to its maximal
immgen <- data.frame(row.names = gsub("\\s", "", immgen$GeneSymbol), immgen[-1])
colnames(immgen) <- gsub("^X.", "", colnames(immgen))
colnames(immgen) <- gsub("\\.", "", colnames(immgen))
immgen <- log2(immgen + 1)

setwd("~/Dropbox/GitHub/Lupus")
save(immgen, file = "./data/immgen.rdt")

# --- Liu et al 2013 JEM data
setwd("~/Dropbox/Lupus/Liu2012")
liu <- read.delim("liu.txt", row.names = 1, sep = "\t", header = T)
id2gene.map <- read.delim("mart_export.txt", sep = "\t", header = T, stringsAsFactors = F)
# Affy MoGene ID to gene name
liu <- merge(id2gene.map, liu, by.x = "Affy.MoGene.probeset", by.y = "row.names", all.y = TRUE)
liu <- aggregate(. ~ Associated.Gene.Name, data = liu[-1], max)  # Aggregate multiple records 
rownames(liu) <- liu$Associated.Gene.Name
liu <- liu[, -1]
colnames(liu) <- c("CPBP1", "CPBP2", "CNBN1", "CNBN2", "CPBN1", "CPBN2")
liu <- liu[apply(liu, 1, function(x) max(x) > 8), ]

setwd("~/Dropbox/GitHub/Lupus")
save(liu, profile, file = "./data/liu.rdt")
