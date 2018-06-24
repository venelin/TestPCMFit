library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(TestPCMFit)
library(surface)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  id <- as.integer(args[1])
} else {
  id <- 1
}


tree <- testData$tree[[id]]
tree$tip.label <- paste0("t", tree$tip.label)
tree$node.label <- NULL
tree <- nameNodes(tree)
values <- testData$X[[id]][, 1:PCMTreeNumTips(tree)]
values <- data.frame(t(values))
rownames(values) <- tree$tip.label


fitSurface <- runSurface(tree, values, plotaic = TRUE, verbose = TRUE)

#save(fitMappings, file = paste0("fitMappings_", id, "_t1_.RData"))


