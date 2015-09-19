#---------------------------------------------------------------------------------------------------
# method <- 2
# 
# if (method == 1) {
# # Method 1: Expressed genes in Immgen data
#   myunion <- merge(roopenian, liu, by = "row.names", all.x = TRUE)  # Union genes
#   myunion <- merge(myunion, immgen.expression, by.x = "Row.names", by.y = "row.names", all.y = TRUE)
# } else if (method == 2) {
# # Method 2: Signature genes in RNA-seq data (DE genes in EBSeq)
#   myunion <- merge(roopenian.signatures, liu, by = "row.names", all.x = TRUE)  # Union genes
#   myunion <- merge(myunion, immgen, by.x = "Row.names", by.y = "row.names", all.x = TRUE)
# } else if (method == 3) {
# # Method 3: Signature genes in RNA-seq data (Log2 FC greater than 2 genes in EBSeq)
#   myunion <- merge(roopenian.signatures.log2fc, liu, by = "row.names", all.x = TRUE)  # Union genes
#   myunion <- merge(myunion, immgen, by.x = "Row.names", by.y = "row.names", all.x = TRUE)
# } else if (method == 4) {
# # Method 4: [Roopenian + Liu]
#   myunion <- merge(roopenian.signatures, liu, by = "row.names", all.x = TRUE)  # Union genes
# } else if (method == 5) {
# # Method 5: [Roopenian + Immgen]
#   myunion <- merge(roopenian.signatures, immgen, by = "row.names", all.x = TRUE)  # Union genes
# } else {
#   stop("Invalid method ID!")
# }
# myunion[is.na(myunion)] <- 0
# myunion <- data.frame(row.names = myunion$Row.names, myunion[-1])


########################### Select and Z-Zanking top correlated samples ############################
# mycor <- cor(myunion, method = "spearman")  # Correlation analysis (method = ["pearson"|"spearman"])
# 
# mycor.vfp.neg.icos.neg.top30 <- 
#   data.frame(Sample = names(sort(mycor[, "vfp.neg.icos.neg"], decreasing = T)[1:30]))
# mycor.vfp.neg.icos.pos.top30 <- 
#   data.frame(Sample = names(sort(mycor[, "vfp.neg.icos.pos"], decreasing = T)[1:30]))
# mycor.vfp.pos.icos.pos.top30 <- 
#   data.frame(Sample = names(sort(mycor[, "vfp.pos.icos.pos"], decreasing = T)[1:30]))
# top30s <- merge(merge(mycor.vfp.neg.icos.neg.top30, mycor.vfp.neg.icos.pos.top30, all = T),
#                 mycor.vfp.pos.icos.pos.top30, all = T)
# myunion.top30s <- myunion[c(as.character(top30s$Sample))]  # Select top 30 correlated samples
# 
# myunion.top30s.rank <- apply(myunion.top30s, 2, rank)  # Ranking
# myunion.top30s.z <- t(apply(myunion.top30s.rank, 1, scale))  # Z-score transformation
# colnames(myunion.top30s.z) <- colnames(myunion.top30s.rank)
# 
# mycor2 <- cor(myunion.top30s.z, method = "pearson")  # Pearson's correlation of Z-transformed data
# mycor2.vfp.neg.icos.neg.top3 <- 
#             data.frame(Sample = names(sort(mycor2[, "vfp.neg.icos.neg"])[49:54]))
# mycor2.vfp.neg.icos.pos.top3 <- 
#             data.frame(Sample = names(sort(mycor2[, "vfp.neg.icos.pos"])[c(46, 47, 48)]))
# mycor2.vfp.pos.icos.pos.top3 <- 
#             data.frame(Sample = names(sort(mycor2[, "vfp.pos.icos.pos"])[c(46, 47, 48)]))
# top3s <- rbind(mycor2.vfp.neg.icos.neg.top3, 
#                mycor2.vfp.neg.icos.pos.top3, 
#                mycor2.vfp.pos.icos.pos.top3)
# myunion.top3s <- myunion[c(as.character(top3s$Sample))]
# 
# myunion.top3s.rank <- apply(myunion.top3s, 2, rank)  # Ranking
# myunion.top3s.z <- t(apply(myunion.top3s.rank, 1, scale))  # Z-score transformation
# colnames(myunion.top3s.z) <- colnames(myunion.top3s.rank)
# 
#---------------------------------------------------------------------------------------------------
# MyCutoff <- function(p1, p2) {
#   #  Return: mean(p1) + sd(p1) | p2 = 1
#   #          mean(p1) - sd(p1) | p2 = 2
#   if (p2 == 1) {
#     return(mean(mypca2$x[, p1]) + sd(mypca2$x[, p1]))
#   } else if (p2 == 2) {
#     return(mean(mypca2$x[, p1]) - sd(mypca2$x[, p1]))
#   } else {
#     stop("Invalid value of p2!")
#   }
#   
# } 
# 
# for (i in 1:5) {  # Top genes in PC 1-5
#   name1 <- paste("MyPCA.1.PC", i, sep = "")  # Positively correlated
#   name2 <- paste("MyPCA.2.PC", i, sep = "")  # Negatively correlated
#   # Extract the expression data
#   assign(name1, myunion.top3s.z[names(mypca2$x[mypca2$x[, i] > MyCutoff(i, 1), i]), ])
#   assign(name2, myunion.top3s.z[names(mypca2$x[mypca2$x[, i] < MyCutoff(i, 2), i]), ])
#   # Save the data
#   write.table(get(name1), file = paste("~/Dropbox/Lupus/", name1, ".txt", sep = ""))
#   write.table(get(name2), file = paste("~/Dropbox/Lupus/", name2, ".txt", sep = ""))
# }
# 
# par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
# plot(apply(MyPCA.2.PC1, 2, mean), type = "b", ylim = c(-3, 3))
# apply(MyPCA.2.PC1, 1, function(x) {lines(x, type = "p", col = colors)})
# 
# plot3d(mypca$x[, 2:4], col = colors, size = 10)  
# # Copyright: Xulong Wang (xulong.wang@jax.org)

library(xlsx)
library(KEGG.db)
library(org.Mm.eg.db)
library(GOstats)  # Conflicts with "igraph"

# --- GO AND KEGG ENRICHMENT ---

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")

load("./data/profile12.rdt")
goAnn <- get("org.Mm.egGO")
universe <- Lkeys(goAnn)

myGO <- list()
for (sample in c("NN", "NP", "PP")) {
  cat(sample, "\n")
  geneId <- profile1[[sample]]
  entrezId <- mget(geneId, org.Mm.egSYMBOL2EG, ifnotfound = NA)
  entrezId <- entrezId[! is.na(entrezId)]
  entrezId <- as.character(entrezId)
  for (category in c("BP", "MF", "CC")) {
    params <- new("GOHyperGParams", geneIds = entrezId, universeGeneIds = universe, annotation = "org.Mm.eg.db", 
                  ontology = category, pvalueCutoff = 0.001, testDirection = "over")  
    over = hyperGTest(params)
    myGO[[sample]][[category]] <- cbind(group = rep(sample, nrow(summary(over))), summary(over))
  }
}

write.xlsx(rbind(myGO$NN$BP, myGO$NP$BP, myGO$PP$BP), 
           file = "Results/tables.xlsx", sheetName = "GO_BP", row.names = F, append = T)
write.xlsx(rbind(myGO$NN$MF, myGO$NP$MF, myGO$PP$MF), 
           file = "Results/tables.xlsx", sheetName = "GO_MF", row.names = F, append = T)
write.xlsx(rbind(myGO$NN$CC, myGO$NP$CC, myGO$PP$CC), 
           file = "Results/tables.xlsx", sheetName = "GO_CC", row.names = F, append = T)

keggAnn <- get("org.Mm.egPATH")
universe <- Lkeys(keggAnn)

myKegg <- list()
for (sample in c("NN", "NP", "PP")) {
  cat(sample, "\n")
  geneId <- profile1[[sample]]
  entrezId <- mget(geneId, org.Mm.egSYMBOL2EG, ifnotfound = NA)
  entrezId <- entrezId[! is.na(entrezId)]
  entrezId <- as.character(entrezId)
  params <- new("KEGGHyperGParams", 
                geneIds=entrezId, universeGeneIds=universe, annotation="org.Mm.eg.db", 
                categoryName="KEGG", pvalueCutoff=0.01, testDirection="over")
  over <- hyperGTest(params)
  kegg <- cbind(group = rep(sample, nrow(summary(over))), summary(over))
  glist <- geneIdsByCategory(over)
  glist <- sapply(glist, function(x) {y <- mget(x, envir=org.Mm.egSYMBOL); paste(y, collapse=";")})
  kegg$Symbols <- glist[as.character(kegg$KEGGID)]
  myKegg[[sample]] <- kegg
}

write.xlsx(rbind(myKegg$NN, myKegg$NP, myKegg$PP), 
           file = "Results/tables.xlsx", sheetName = "KEGG", row.names = F, append = T)
