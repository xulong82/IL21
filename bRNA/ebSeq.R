# Copyright: Xulong Wang (xulong.wang@jax.org)

library(EBSeq)

rm(list = ls())
setwd("~/Dropbox/GitHub/Lupus")
load(file = "./data/myCount.rdt")

expr <- myCount[apply(myCount, 1, function(x) max(x) > 6), ]
expr <- as.matrix(expr)
size <- MedianNorm((expr))
condition <- gsub("[12]", "", colnames(expr))
patterns <- GetPatterns(condition)
PlotPattern(patterns)
ebout <- EBMultiTest(expr, Conditions = condition, sizeFactors = size, maxround = 5)
QQP(ebout)

pp <- GetMultiPP(ebout)$PP
fc <- GetMultiFC(ebout)$Log2FCMat
average <- GetMultiFC(ebout)$CondMeans
complete <- cbind(average, fc, pp)

sig <- NULL
cut <- 0.75
for (i in 2:5) {
  name <- paste("Pattern", i, sep = "")
  assign(name, subset(complete, complete[, name] >= cut))
  sig <- rbind(sig, get(name))
} 

geneId <- rownames(sig)

write.table(complete, file = "data/ebseq.txt")
save(ebout, sig, file = "./data/ebseq.rdt")
load("./data/ebseq.rdt")

# --- PROLIFIC GENES FROM EBSEQ ---
profile1 <- profile2 <- list()

profile1$NN <- geneId[sig[, "NN"] > 50 & sig[, "NNOverNP"] > 1 & sig[, "NNOverPP"] > 1]
profile1$NP <- geneId[sig[, "NP"] > 50 & sig[, "NNOverNP"] < -1 & sig[, "NPOverPP"] > 1]
profile1$PP <- geneId[sig[, "PP"] > 50 & sig[, "NNOverPP"] < -1 & sig[, "NPOverPP"] < -1]

profile2$NN <- profile1$NN[sig[profile1$NN, "Pattern4"] > 0.75]
profile2$NP <- profile1$NP[sig[profile1$NP, "Pattern3"] > 0.75]
profile2$PP <- profile1$PP[sig[profile1$PP, "Pattern2"] > 0.75]
