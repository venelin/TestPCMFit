library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(data.table)
library(phytools)
library(ggtree)


source("./R-local/GeneratePCMModels.R")
set.seed(2)

# number of regimes
R <- 2
# number of traits
k <- 2

# this will generate a non-ultrametric tree with about 300 tips
treeFossil <- pbtree(n=200, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossil)
# normalize branch lengths, so that the length of the tree equals 1
treeFossil$edge.length <- treeFossil$edge.length / max(PCMTreeNodeTimes(treeFossil))

# this will generate an ultrametric tree with exactly 300 tips
treeExtant <- pbtree(n=318, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtant)
# normalize branch lengths, so that the length of the tree equals 1
treeExtant$edge.length <- treeExtant$edge.length / max(PCMTreeNodeTimes(treeExtant))

# PCMTreePlot(treeFossil) + geom_nodelab()
# PCMTreePlot(treeExtant) + geom_nodelab()

# number of different modelMappings per per clustering
nRandomModelMappingsPerClustering <- 4

# number of generated random models per model-mapping  of a tree
nRandomParamsPerMapping <- 3

# number of simulations per random model
nSimulationsPerRandomParam <- 2

# number of traits
k <- 2

testData <- rbindlist(
  list(
    data.table(
      treeType = "extant",
      tree = list(treeExtant),
      numClusters = 2,
      clusterNodes = list(as.character(c(319, 499))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "extant",
      tree = list(treeExtant),
      numClusters = 8,
      clusterNodes = list(as.character(c(319, 499, 438, 360,
                                         320, 486, 376, 583))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    ),
    data.table(
      treeType = "fossil",
      tree = list(treeFossil),
      numClusters = 2,
      clusterNodes = list(as.character(c(319, 393))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "fossil",
      tree = list(treeFossil),
      numClusters = 8,
      clusterNodes = list(as.character(c(319, 393, 430, 484,
                                         517, 486, 600, 343))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    )

    ))


testData[, treeWithRegimes:=lapply(1:.N, function(i) {
  PCMTreeSetRegimes(tree[[i]], nodes = clusterNodes[[i]], inplace = FALSE)
})]

# this will reorder the nodes in clusterNodes column according to the preorder traversal of
# the tree
testData[, clusterNodes:=lapply(1:.N, function(i) {
  PCMTreeGetLabels(treeWithRegimes[[i]])[PCMTreeGetStartingNodesRegimes(treeWithRegimes[[i]])]
})]

# see how the treees with regimes look like
# PCMTreePlot(testData$treeWithRegimes[[2]]) + geom_nodelab()
# PCMTreePlot(testData$treeWithRegimes[[4]]) + geom_nodelab()


# replicate each row in testData nRandomParamsPerMapping times
# replicate each row in testData nSimulationsPerRandomParam times
testData <- testData[sort(rep_len(1:.N, length.out = .N * nRandomParamsPerMapping))]

options(PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal = 0.2)
options(PCMBase.ParamValue.LowerLimit = -2)
options(PCMBase.ParamValue.UpperLimit = 10)



testData[, model:=lapply(1:.N, function(i) {
  model <- do.call(
    MixedGaussian,
    c(list(k = k, modelTypes = simulatedModels, mapping = mapping[[i]]),
      argsMixedGaussian_SimulatedModels))

  vecParams <- PCMParamRandomVecParams(model, n = 1)
  vecParams <- round(vecParams, digits = 1)
  PCMParamLoadOrStore(model, vecParams, 0, load = TRUE)
  model
})]

# replicate each row in testData nSimulationsPerRandomParam times
testData <- testData[sort(rep_len(1:.N, length.out = .N * nSimulationsPerRandomParam))]

# simulate the trait values
testData[, X:=lapply(1:.N, function(i) {
  cat(i, ", ")
  PCMSim(treeWithRegimes[[i]], model[[i]], X0 = model[[i]]$X0)
})]

# log-likelihood calculated at the model used to generate the data
testData[, logLik:=lapply(1:.N, function(i) {
  cat(i, ", ")
  PCMLik(X[[i]], treeWithRegimes[[i]], model[[i]], PCMInfoCpp(X[[i]], treeWithRegimes[[i]], model[[i]]))
})]

# AIC calculated at the model used to generate the data
testData[, AIC:=lapply(1:.N, function(i) {
  cat(i, ", ")
  model <- model[[i]]
  attr(model, "tree") <- treeWithRegimes[[i]]
  attr(model, "X") <- X[[i]]
  attr(model, "PCMInfoFun") <- PCMInfoCpp
  AIC(model)
})]

testData[, nobs:=sapply(treeWithRegimes, PCMTreeNumTips)]
testData[, df:=sapply(model, PCMParamCount, TRUE,TRUE)]

devtools::use_data(testData, overwrite = TRUE)
devtools::use_data(treeExtant, overwrite = TRUE)
devtools::use_data(treeFossil, overwrite = TRUE)

