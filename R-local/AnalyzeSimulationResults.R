library(PCMFit)
library(TestPCMFit)
library(ggtree)
library(ggplot2)
library(ggimage)
library(cowplot)
library(PCMBase)
library(data.table)
library(ROCR)
library(fpc)

GeneratePCMModels()

testResultData <- rbindlist(
  lapply(1:nrow(testData), function(i) {
    resultFile <- paste0("../local-data/FinalResult_MixedGaussian_testData_id_", i, "_t1_.RData")

    if(file.exists(resultFile)) {
      cat("Loading ", resultFile, "\n")
      load(resultFile)
      bestFitAIC <- RetrieveBestFitAIC(fitMappings)
      trueFromTestData <- RetrieveTrueFromTestData(i)

      matPredicted_Cluster <- PCMTreeMatrixNodesInSameRegime(bestFitAIC$tree)
      mode(matPredicted_Cluster) <- "double"
      matTRUE_Cluster <- PCMTreeMatrixNodesInSameRegime(trueFromTestData$tree)
      mode(matTRUE_Cluster) <- "double"

      perf_Cluster <- performance(prediction(matPredicted_Cluster, matTRUE_Cluster), "tpr", "fpr")

      # work-around because ROCR doesn't support all true classes equaling the same class
      perf_BM <- performance(prediction(
        predictions = c(is.BM(bestFitAIC$inferredMappedModels), c(0,1)),
        labels = c(is.BM(trueFromTestData$trueMappedModels), c(0,1))
      ), "tpr", "fpr")

      perf_Uncorrelated <- performance(prediction(
        predictions = c(is.Uncorrelated(bestFitAIC$inferredMappedModels), c(0, 1)),
        labels = c(is.Uncorrelated(trueFromTestData$trueMappedModels), c(0, 1))
      ), "tpr", "fpr")

      perf_NonDiagonalH <- performance(prediction(
        predictions = c(is.NonDiagonalH(bestFitAIC$inferredMappedModels), c(0, 1)),
        labels = c(is.NonDiagonalH(trueFromTestData$trueMappedModels), c(0, 1))
      ), "tpr", "fpr")

      perf_AsymetricH <- performance(prediction(
        predictions = c(is.NonSymmetricH(bestFitAIC$inferredMappedModels), c(0, 1)),
        labels = c(is.NonSymmetricH(trueFromTestData$trueMappedModels), c(0, 1))
      ), "tpr", "fpr")


      data.table(i = i,
                 treeType=testData[i, treeType],
                 numClusters = testData[i, numClusters],
                 crit = factor(c("Cluster",
                                 "BM process",
                                 "Uncorrelated traits",
                                 "NonDiagonal H",
                                 "Asymetric H")
                               ),
                 fpr = c(
                   perf_Cluster@x.values[[1]][2],
                   perf_BM@x.values[[1]][2],
                   perf_Uncorrelated@x.values[[1]][2],
                   perf_NonDiagonalH@x.values[[1]][2],
                   perf_AsymetricH@x.values[[1]][2]
                 ),
                 tpr = c(
                   perf_Cluster@y.values[[1]][2],
                   perf_BM@y.values[[1]][2],
                   perf_Uncorrelated@y.values[[1]][2],
                   perf_NonDiagonalH@y.values[[1]][2],
                   perf_AsymetricH@y.values[[1]][2]
                 ),
                 AIC_Final = AIC(bestFitAIC$inferredModel),
                 AIC_True = AIC(trueFromTestData$trueModel)
                 )
    } else {
      NULL
    }
  })
)

testResultData[, crit2:=factor(crit, levels = c("Cluster",
                                                "BM process",
                                                "Uncorrelated traits",
                                                "NonDiagonal H",
                                                "Asymetric H"))]
#ggplot(testResultData[AIC_Final<AIC_True]) +
ggplot(testResultData) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.2, color="red") +
  geom_label(aes(label = i, x = fpr, y = tpr,
                 size = 0.5*apply(cbind((0.2+1-tpr), (0.2+fpr)), 1, max),
                 color = 0.5*apply(cbind((0.2+1-tpr), (0.2+fpr)), 1, max),
                 fill = AIC_Final<AIC_True ),
             position = position_jitter(width=0.1, height = 0.1),
             label.padding = unit(0.1, "lines"), fontface = "bold") +
  scale_color_continuous(low="green", high="red") +
  scale_fill_manual(values = c("TRUE"="white", "FALSE"="black")) +
  scale_size_continuous(range = c(1.5, 3)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1)) +
  facet_grid(numClusters*treeType~crit2) +
  theme_bw() +
  theme(legend.position = "none")

i <- 54
resultFile <- paste0("local-data/FinalResult_MixedGaussian_testData_id_", i, "_t1_.RData")
load(resultFile)
#fit <- RetrieveBestFitAIC(fitMappings)

#plList <- PlotCompareFitWithTrueTestData(fitMappings, i, doPlot2DData = FALSE)
#plot_grid(plotlist = plList)

bestFitAIC <- RetrieveBestFitAIC(fitMappings)
logLik(bestFitAIC$inferredModel)

MeanVec <- PCMMean(bestFitAIC$tree, bestFitAIC$inferredModel, bestFitAIC$inferredModel$X0)
VarMat <- PCMVar(bestFitAIC$tree, bestFitAIC$inferredModel)
VarMat <- 0.5*(VarMat + t(VarMat))

mvtnorm::dmvnorm(as.vector(bestFitAIC$X), as.vector(MeanVec), VarMat, log=TRUE)


tree3 <- attr(bestFitAIC$inferredModel, "tree")
#for(epoch in seq(0.01, max(PCMTreeNodeTimes(tree3)), by = 0.02)) {
  #points <- PCMTreeLocateEpochOnBranches(tree3, epoch)
  points <- PCMTreeLocateMidpointsOnBranches(tree3, 0.2)
  tree3 <- PCMTreeInsertSingletons(tree3, points$nodes, points$positions)
#}

inferredModel3 <- bestFitAIC$inferredModel
attr(inferredModel3, "tree") <- tree3
logLik(inferredModel3)

MeanVec3 <- PCMMean(tree3, PCMApplyTransformation(inferredModel3), inferredModel3$X0)
VarMat3 <- PCMVar(tree3, PCMApplyTransformation(inferredModel3))
VarMat3 <- 0.5*(VarMat3 + t(VarMat3))
eigen(VarMat3)$values[-(1:600)]

mvtnorm::dmvnorm(as.vector(bestFitAIC$X), as.vector(MeanVec3), VarMat3, log=TRUE)

true <- RetrieveTrueFromTestData(i)

MeanVecTrue <- PCMMean(true$tree, true$trueModel, true$trueModel$X0)
VarMatTrue <- PCMVar(true$tree, true$trueModel)
VarMatTrue <- 0.5*(VarMatTrue+t(VarMatTrue))

logLik(true$trueModel)

mvtnorm::dmvnorm(as.vector(bestFitAIC$X), as.vector(MeanVecTrue), VarMatTrue, log=TRUE)


# bhattacharyya.dist(as.vector(MeanVecTrue), as.vector(MeanVec), VarMatTrue, VarMat)
# bhattacharyya.dist(as.vector(MeanVecTrue), as.vector(MeanVecTrue), VarMatTrue, VarMatTrue)




resultFile <- paste0("local-data/FinalResult_MixedGaussian_MammalData_id_1_t1_.RData")
load(resultFile)
#pl <- PlotSearchHistory(fitMappings)

bestFitAIC <- RetrieveBestFitAIC(fitMappings)
#options(PCMBase.Skip.Singular=TRUE, PCMBase.Threshold.SV=1e-6, PCMBase.Threshold.EV=1e-6, PCMBase.Threshold.Skip.Singular=0.2)
logLik(bestFitAIC$inferredModel)

library(PCMBaseCpp)
PCMLik(bestFitAIC$X, bestFitAIC$tree, bestFitAIC$inferredModel, PCMInfoCpp(bestFitAIC$X, bestFitAIC$tree, bestFitAIC$inferredModel))

treeMammals3 <- attr(bestFitAIC$inferredModel, "tree")

points <- PCMTreeLocateMidpointsOnBranches(treeMammals3, 15)
treeMammals3 <- PCMTreeInsertSingletons(treeMammals3, points$nodes, points$positions)


inferredModel3 <- bestFitAIC$inferredModel
attr(inferredModel3, "tree") <- treeMammals3


logLik(inferredModel3)

MeanVecMammals3 <- PCMMean(treeMammals3, PCMApplyTransformation(inferredModel3), inferredModel3$X0)
VarMatMammals3 <- PCMVar(treeMammals3, PCMApplyTransformation(inferredModel3))
VarMatMammals3 <- 0.5*(VarMatMammals3 + t(VarMatMammals3))
eigen(VarMatMammals3)$values[-(1:1200)]

#options(PCMBase.Skip.Singular=TRUE, PCMBase.Threshold.SV=1e-4, PCMBase.Threshold.Skip.Singular=0.2)
#logLik(bestFitAIC$inferredModel)

MeanVecMammals <- PCMMean(bestFitAIC$tree, PCMApplyTransformation(bestFitAIC$inferredModel), bestFitAIC$inferredModel$X0)
VarMatMammals <- PCMVar(bestFitAIC$tree, PCMApplyTransformation(bestFitAIC$inferredModel))
VarMatMammals <- 0.5*(VarMatMammals + t(VarMatMammals))

mvtnorm::dmvnorm(as.vector(valuesMammals), as.vector(MeanVecMammals), VarMatMammals, log=TRUE)

library(LaplacesDemon)

VarMatMammals2 <- Matrix::nearPD(VarMatMammals)
dmvn(as.vector(valuesMammals), as.vector(MeanVecMammals), as.matrix(VarMatMammals2$mat), log=TRUE)


eigen(VarMatMammals)$values[-(1:1000)]
library(rARPACK)
lambda <- eigs_sym(VarMatMammals, k = 5, which = "SM")



