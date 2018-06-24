library(PCMFit)
library(TestPCMFit)
library(PCMBase)
library(PCMBaseCpp)
library(data.table)

GeneratePCMModels()

options(PCMBase.Value.NA = -1e20)
options(PCMBase.Lmr.mode = 11)

load("FinalResult_MixedGaussian_MammalData_id_1_t3_.RData")


tree <- treeMammals

# segment long branches
while(TRUE) {
  points <- PCMTreeLocateMidpointsOnBranches(tree, 16)
  if(length(points$nodes) == 0) {
    break
  } else {
    tree <- PCMTreeInsertSingletons(tree, points$nodes, points$positions)
  }
}

inferred <- RetrieveBestFitAIC(fitMappings)
inferredModel <- inferred$inferredModel
inferredTree <- attr(inferredModel, "tree")

PCMTreeSetLabels(inferredTree, PCMTreeGetLabels(tree))

colnames(valuesMammals) <- inferredTree$tip.label


AIC(inferredModel)
logLik(inferredModel)

argsPCMParamLowerLimit <- list()
argsPCMParamUpperLimit <- list()

bestFitRep <- PCMFit(valuesMammals, inferredTree, inferredModel,
                     metaIFun = PCMInfoCpp, positiveValueGuard = 1000,
                     argsPCMParamLowerLimit = argsPCMParamLowerLimit,
                     argsPCMParamUpperLimit = argsPCMParamUpperLimit,
                     argsConfigOptimAndMCMC = list(
                       nCallsOptim = 1000, genInitNumEvals = 2000000, genInitVerbose = FALSE),
                     verbose = TRUE)

save(bestFitRep, file = "RepeatBestFitOnMammalsData_t3.RData")
