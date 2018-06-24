library(ape)
library(PCMBase)
library(PCMBaseCpp)
library(PCMFit)
library(data.table)
library(TestPCMFit)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  id <- as.integer(args[1])
} else {
  id <- 1
}

if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a cluster)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else {
    cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE),
                                     outfile = paste0("log_", id, ".txt"))
    doParallel::registerDoParallel(cluster)
  }

}

tree <- treeMammals
PCMTreeSetLabels(tree)
tree$edge.length <- tree$edge.length / max(PCMTreeNodeTimes(tree))
values <- valuesMammals[, 1:PCMTreeNumTips(tree)]

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

tableFits <- NULL

#uncomment if re-running
load("tableFits_Mammals.RData")

tableFits2 <- rbindlist(lapply(list.files(".", pattern=paste0("fits_Mammals_t2_*")), function(f)  {
  load(f)
  data
}))

tableFits4 <- rbindlist(lapply(list.files(".", pattern=paste0("fits_Mammals_t4_*")), function(f)  {
  load(f)
  data
}))

tableFits5 <- rbindlist(lapply(list.files(".", pattern=paste0("fits_Mammals_t5_*")), function(f)  {
  load(f)
  data
}))

tableFits <- rbindlist(list(tableFits, tableFits2, tableFits4, tableFits5))
setkey(tableFits, hashCodeTree, hashCodeStartingNodesRegimesLabels, hashCodeMapping)

print(nrow(tableFits))
tableFits <- unique(tableFits, by = key(tableFits))
tableFits[, duplicated:=FALSE]
print(tableFits[treeEDExpression=="tree"][order(aic)])

# setnames(tableFits, "AIC", "aic")
# tableFits[, mapping:=lapply(mapping, function(m) modelTypes[m])]
# tableFits[, hashCodeMapping:=sapply(1:.N, function(i) {
#  if(i%%100==0) cat("i=",i,"\n")
#  modelMapping <- match(mapping[[i]], modelTypes)
#  #tree, modelTypes, startingNodesRegimesLabels, modelMapping
#  treei <- PCMTreeEvalNestedEDxOnTree(treeEDExpression[[i]], tree)
#  PCMTreeSetRegimes(treei, startingNodesRegimesLabels[[i]])
#  hashCodes <- HashCodes(treei, modelTypes, startingNodesRegimesLabels[[i]], modelMapping)
#  hashCodes$hashCodeMapping
#})]
#
#save(tableFits, file="tableFits_Mammals.RData")
#
# print(tableFits)
# tableFits[, duplicated:=FALSE]
# if(nrow(tableFits) == 0) {
#   tableFits <- NULL
# }

argsMRG <- list(
  X0 = list(default = rep(0, 2), type = c("gvector", "full")),
  Sigmae_x = list(default = matrix(0, 2, 2), type = c("gmatrix", "fixed"),
                  description = "Fixed upper triangular Choleski factor of the variance-covariance matrix for the non-phylogenetic trait component"))

argsPCMLowerBound <- list(lowerBoundValue = -10, lowerBoundValuePositiveDiag = 0)
argsPCMUpperBound <- list(upperBoundValue = 10, upperBoundValuePositiveDiag = 10)

fitMappings <- PCMFitModelMappings(
  values, tree, modelTypes = modelTypes,
  metaIFun = PCMInfoCpp, positiveValueGuard = 1000,

  tableFitsPrev = tableFits,

  prefixFiles = paste0("fits_Mammals_t6_"),
  maxCladePartitionLevel = 4, minCladeSizes = 20,

  argsMRG = argsMRG,
  argsPCMLowerBound = argsPCMLowerBound,
  argsPCMUpperBound = argsPCMUpperBound,
  argsConfigOptimAndMCMC1 = list(nCallsOptim = 200, genInitNumEvals = 100000, genInitVerbose = FALSE),
  argsConfigOptimAndMCMC2 = list(nCallsOptim = 10, genInitNumEvals = 2000, genInitVerbose = FALSE),

  numJitterAllRegimeFits = 1000, numJitterRootRegimeFit = 1000,

  printFitVectorsToConsole = FALSE,
  doParallel = TRUE,
  verbose = TRUE,
  verbosePCMFit = FALSE)

save(fitMappings, file = paste0("fitMappings_Mammals_t6_.RData"))

if(exists("cluster") && !is.null(cluster)) {
  parallel::stopCluster(cluster)
  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.

  cluster <- NULL
}

