library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(data.table)
library(phytools)
library(ggtree)
library(TestPCMFit)

GeneratePCMModels()

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

simulatedModels <- SimulatedModels()
argsMixedGaussian_SimulatedModels <- ArgsMixedGaussian_SimulatedModels()


# We overwrite the limits for the model parameters, setting them to narrower sub-region
# of the limits specified in TestPCMFit/GeneratePCMModels. This prevents generating
# data for which some of the transition covariance matrices are singular due to
# too extreme random parameters. Note that these limits have been chosent to be adequate with
# the time-scale of the trees (all trees are of depth 166.2)

# these options will affect the H matrix, so that we don't have to specify it explicitly
# in the functions below
options(PCMBase.ParamValue.LowerLimit = -4,
        PCMBase.ParamValue.UpperLimit = 4,
        PCMBase.ParamValue.LowerLimit.NonNegativeDiagonal = .1)

PCMParamLowerLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[1, 1] <- .05
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- .0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[1, 1, r] <- .05
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .0
      }
    }
  }
  o
}

PCMParamUpperLimit.BM <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .5
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- .2
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .5
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .2
      }
    }
  }
  o
}

PCMParamLowerLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Theta)) {
    o$Theta[1] <- 3.0
    o$Theta[2] <- 2.0
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 3.0
      o$Theta[2, r] <- 2.0
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .05
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- -.0
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .05
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .0
      }
    }
  }
  o
}

PCMParamUpperLimit.OU <- function(o, k, R, ...) {
  o <- NextMethod()
  k <- attr(o, "k", exact = TRUE)
  R <- length(attr(o, "regimes", exact = TRUE))

  if(is.Global(o$Theta)) {
    o$Theta[1] <- 6.0
    o$Theta[2] <- 4.0
  } else {
    for(r in seq_len(R)) {
      o$Theta[1, r] <- 6.0
      o$Theta[2, r] <- 4.0
    }
  }
  if(is.Global(o$Sigma_x)) {
    o$Sigma_x[1, 1] <- o$Sigma_x[2, 2] <- .5
    if(!is.Diagonal(o$Sigma_x)) {
      o$Sigma_x[1, 2] <- .2
    }
  } else {
    for(r in seq_len(R)) {
      o$Sigma_x[1, 1, r] <- o$Sigma_x[2, 2, r] <- .5
      if(!is.Diagonal(o$Sigma_x)) {
        o$Sigma_x[1, 2, r] <- .2
      }
    }
  }
  o
}



set.seed(2)

# number of regimes
R <- 2
# number of traits
k <- 2

# this will generate a non-ultrametric tree with 318 tips
treeFossilSmall <- pbtree(n=200, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossilSmall)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossilSmall$edge.length <- treeFossilSmall$edge.length / max(PCMTreeNodeTimes(treeFossilSmall)) * 166.2
PCMTreeSetDefaultRegime(treeFossilSmall, regime = 1)

# this will generate an ultrametric tree with exactly 300 tips
treeExtantSmall <- pbtree(n=318, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtantSmall)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtantSmall$edge.length <- treeExtantSmall$edge.length / max(PCMTreeNodeTimes(treeExtantSmall)) * 166.2
PCMTreeSetDefaultRegime(treeExtantSmall, regime = 1)

# this will generate a non-ultrametric tree with 638 tips
treeFossilBig <- pbtree(n=374, scale=1, b = 1, d = 0.4)
PCMTreeSetLabels(treeFossilBig)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeFossilBig$edge.length <- treeFossilBig$edge.length / max(PCMTreeNodeTimes(treeFossilBig)) * 166.2
PCMTreeSetDefaultRegime(treeFossilBig, regime = 1)

# this will generate an ultrametric tree with exactly 300 tips
treeExtantBig <- pbtree(n=638, scale=1, b = 1, d = 0.4, extant.only = TRUE)
PCMTreeSetLabels(treeExtantBig)
# set the depth of the tree to the depth of the mammal tree (166.2)
treeExtantBig$edge.length <- treeExtantBig$edge.length / max(PCMTreeNodeTimes(treeExtantBig)) * 166.2
PCMTreeSetDefaultRegime(treeExtantBig, regime = 1)


# PCMTreePlot(treeFossilSmall, layout="fan") + geom_nodelab()
# PCMTreePlot(treeExtantSmall, layout="fan") + geom_nodelab()
# PCMTreePlot(treeFossilBig, layout="fan") + geom_nodelab(size = 3)
# PCMTreePlot(treeExtantBig, layout="fan") + geom_nodelab(size = 3)


# number of different modelMappings per per clustering
nRandomModelMappingsPerClustering <- 4

# number of generated random models per model-mapping  of a tree
nRandomParamsPerMapping <- 3

# number of simulations per random model
nSimulationsPerRandomParam <- 2

# number of traits
k <- 2

testData_t2 <- rbindlist(
  list(
    data.table(
      treeType = "extant-small",
      tree = list(treeExtantSmall),
      numClusters = 2,
      clusterNodes = list(as.character(c(319, 499))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "extant-small",
      tree = list(treeExtantSmall),
      numClusters = 8,
      clusterNodes = list(as.character(c(319, 499, 438, 360,
                                         320, 486, 376, 583))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    ),
    data.table(
      treeType = "fossil-small",
      tree = list(treeFossilSmall),
      numClusters = 2,
      clusterNodes = list(as.character(c(319, 393))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "fossil-small",
      tree = list(treeFossilSmall),
      numClusters = 8,
      clusterNodes = list(as.character(c(319, 393, 430, 484,
                                         517, 486, 600, 343))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    ),

    # big trees
    data.table(
      treeType = "extant-big",
      tree = list(treeExtantBig),
      numClusters = 2,
      clusterNodes = list(as.character(c(639, 1007))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "extant-big",
      tree = list(treeExtantBig),
      numClusters = 8,
      clusterNodes = list(as.character(c(639, 1105, 1007, 817,
                                         867, 700, 1236, 1177))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    ),
    data.table(
      treeType = "fossil-big",
      tree = list(treeFossilBig),
      numClusters = 2,
      clusterNodes = list(as.character(c(639, 883))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 2, replace = TRUE))
    ),
    data.table(
      treeType = "fossil-big",
      tree = list(treeFossilBig),
      numClusters = 8,
      clusterNodes = list(as.character(c(639, 685, 644, 641,
                                         914, 971, 1136, 975))),
      mapping = lapply(1:nRandomModelMappingsPerClustering, function(i) sample(1:length(simulatedModels), size = 8, replace = TRUE))
    )

    ))


testData_t2[, treeWithRegimes:=lapply(1:.N, function(i) {
  PCMTreeSetRegimes(tree[[i]], nodes = clusterNodes[[i]], inplace = FALSE)
})]

# this will reorder the nodes in clusterNodes column according to the preorder traversal of
# the tree
testData_t2[, clusterNodes:=lapply(1:.N, function(i) {
  PCMTreeGetLabels(treeWithRegimes[[i]])[PCMTreeGetStartingNodesRegimes(treeWithRegimes[[i]])]
})]

# see how the treees with regimes look like
# PCMTreePlot(testData_t2$treeWithRegimes[[2]]) + geom_nodelab()
# PCMTreePlot(testData_t2$treeWithRegimes[[4]]) + geom_nodelab()


# replicate each row in testData_t2 nRandomParamsPerMapping times
# replicate each row in testData_t2 nSimulationsPerRandomParam times
testData_t2 <- testData_t2[sort(rep_len(1:.N, length.out = .N * nRandomParamsPerMapping))]

# generate random models
testData_t2[, model:=lapply(1:.N, function(i) {
  model <- do.call(
    MixedGaussian,
    c(list(k = k, modelTypes = simulatedModels, mapping = mapping[[i]]),
      argsMixedGaussian_SimulatedModels))

  vecParams <- PCMParamRandomVecParams(model, n = 1)
  vecParams <- round(vecParams, digits = 1)
  PCMParamLoadOrStore(model, vecParams, 0, load = TRUE)
  # use a fixed value for X0 in all simulations
  model$X0[] <- c(1.0, -1.0)
  model
})]

# replicate each row in testData_t2 nSimulationsPerRandomParam times
testData_t2 <- testData_t2[sort(rep_len(1:.N, length.out = .N * nSimulationsPerRandomParam))]

# simulate the trait values
testData_t2[, X:=lapply(1:.N, function(i) {
  cat(i, ", ")
  PCMSim(treeWithRegimes[[i]], model[[i]], X0 = model[[i]]$X0)
})]

# for 10 simulations there was an eigenvalue below the threshold of 1e-5 for the
# variance covariance transition matrix on 1 node. To calculate the likelihood
# and AIC for these 10 simulations we set the threshold to a more tolerant value.
# In the inference, the threshold will be left as per defaul (i.e. 1e-5).
options(PCMBase.Threshold.EV = 1e-7)

# log-likelihood calculated at the model used to generate the data
testData_t2[, logLik:=lapply(1:.N, function(i) {
  cat(i, ", ")
  PCMLik(X[[i]], treeWithRegimes[[i]], model[[i]], PCMInfoCpp(X[[i]], treeWithRegimes[[i]], model[[i]]))
})]

# AIC calculated at the model used to generate the data
testData_t2[, AIC:=lapply(1:.N, function(i) {
  cat(i, ", ")
  model <- model[[i]]
  attr(model, "tree") <- treeWithRegimes[[i]]
  attr(model, "X") <- X[[i]]
  attr(model, "PCMInfoFun") <- PCMInfoCpp
  AIC(model)
})]

testData_t2[, nobs:=sapply(treeWithRegimes, PCMTreeNumTips)]
testData_t2[, df:=sapply(model, PCMParamCount, TRUE,TRUE)]

# store the data within the package TestPCMFit
devtools::use_data(testData_t2, overwrite = TRUE)

