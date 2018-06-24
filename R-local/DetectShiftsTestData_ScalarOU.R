# This script can be run locally or on a cluster using a command like:
#
# bsub -M 100000 -n 100 -W 23:59 -R ib sh R --vanilla --slave -f DetectShiftsTestData_XXX.R --args 96
#
# In the above command the argument 96 after --args denotes the row-number in testData data.table from the
# package TestPCMFit. The blank XXX stays for the type of model to be inferred, e.g. ScalarOU or MixedGaussian
#
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

prefixFiles = paste0("ScalarOU_testData_id_", id, "_t1_")

if(!exists("cluster") || is.null(cluster)) {
  if(require(doMPI)) {
    # using MPI cluster as distributed node cluster (possibly running on a cluster)
    # Get the number of cores. Assume this is run in a batch job.
    p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
    cluster <- startMPIcluster(count = p-1, verbose = TRUE)
    doMPI::registerDoMPI(cluster)
  } else {
    cluster <- parallel::makeCluster(parallel::detectCores(logical = TRUE),
                                     outfile = paste0("log_", prefixFiles, ".txt"))
    doParallel::registerDoParallel(cluster)
  }
}

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

tableFits <- NULL
# try using a previously stored tableFits from a previous run that was interupted
fileCurrentResults <- paste0("CurrentResults_", prefixFiles, ".RData")
if(file.exists(fileCurrentResults)) {
  cat("Loading previously stored tableFits from file", fileCurrentResults, "...\n")
  load(fileCurrentResults)


  tempFiles <- list.files(pattern = paste0("^", prefixFiles, ".*.RData"))
  if(length(tempFiles) > 0) {
    cat("Loading previously stored tableFits from temporary files (", toString(tempFiles), ")...\n")
    tableFitsTempFiles <- rbindlist(
      lapply(tempFiles, function(file) {
        load(file)
        data
      }))
    tableFits <- rbindlist(list(tableFits, tableFitsTempFiles))
  }

  setkey(tableFits, hashCodeTree,hashCodeStartingNodesRegimesLabels,hashCodeMapping)

  tableFits <- unique(tableFits, by = key(tableFits))
  tableFits[, duplicated:=FALSE]

  print(tableFits)
}

modelTypes <- InferredModel_ScalarOU()
argsMixedGaussian <- ArgsMixedGaussian_SurfaceOU()
argsPCMParamLowerLimit <- list()
argsPCMParamUpperLimit <- list()

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

tree <- testData$tree[[id]]
values <- testData$X[[id]][, 1:PCMTreeNumTips(tree)]


fitMappings <- PCMFitModelMappings(
  values, tree, modelTypes = modelTypes,
  generatePCMModelsFun = GeneratePCMModels,
  metaIFun = PCMInfoCpp, positiveValueGuard = 1000,

  tableFits = tableFits,

  prefixFiles = prefixFiles,

  maxCladePartitionLevel = 100, maxNumNodesPerCladePartition = 1, minCladeSizes = 20,

  argsMixedGaussian = argsMixedGaussian,
  argsPCMParamLowerLimit = argsPCMParamLowerLimit,
  argsPCMParamUpperLimit = argsPCMParamUpperLimit,
  argsConfigOptimAndMCMC1 = list(nCallsOptim = 200, genInitNumEvals = 100000, genInitVerbose = FALSE),
  argsConfigOptimAndMCMC2 = list(nCallsOptim = 10, genInitNumEvals = 2000, genInitVerbose = FALSE),

  numJitterAllRegimeFits = 1000, numJitterRootRegimeFit = 1000,

  printFitVectorsToConsole = TRUE,
  doParallel = TRUE,
  verbose = TRUE)

save(fitMappings, file = paste0("FinalResult_", prefixFiles, ".RData"))

if(exists("cluster") && !is.null(cluster)) {
  parallel::stopCluster(cluster)
  # Don't forget to destroy the parallel cluster to avoid leaving zombie worker-processes.

  cluster <- NULL
}

