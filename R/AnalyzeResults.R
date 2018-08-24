

RetrieveTrueFromTestData <- function(testData, i) {
  res <- list(
    tree = testData$tree[[i]],
    X = testData$X[[i]][, 1:PCMTreeNumTips(testData$tree[[i]])],
    modelTypes = SimulatedModels(),
    trueRegimeNodes = testData$clusterNodes[[i]],
    trueMapping = SimulatedModels()[testData$mapping[[i]]],
    trueMappingIdx = testData$mapping[[i]],
    trueModel = testData$model[[i]]
  )
  PCMTreeSetLabels(res$tree)
  PCMTreeSetRegimes(res$tree, res$trueRegimeNodes)

  res[["trueMappedModels"]] <- attr(res$trueModel, "mapping")[res$tree$edge.regime]
  attr(res$trueModel, "X") <- res$X
  attr(res$trueModel, "tree") <- res$tree
  res
}

is.BM <- function(x) as.double(x <= 2)
is.OU <- function(x) as.double(x >= 3)
is.Uncorrelated <- function(x) as.double(x %in% c(1, 3))
is.Correlated <- function(x) as.double(x %in% c(2, 4, 5, 6))
is.NonDiagonalH <- function(x) as.double(x > 4)
is.NonSymmetricH <- function(x) as.double(x == 6)

PlotCompareFitWithTrueTestData <- function(fitMappings, i, doPlotSearchHistory = TRUE, doPlot2DData = TRUE) {
  tree <- testData$tree[[i]]
  N <- PCMTreeNumTips(tree)
  trueRegimeNodes <- testData$clusterNodes[[i]]
  trueMappingIdx <- testData$mapping[[i]]
  plTrue <- PlotTreeRegimesAndMapping(tree, regimeNodes = trueRegimeNodes, mappingIdx = trueMappingIdx)

  trueModel <- testData$model[[i]]
  trueTree <- PCMTreeSetRegimes(tree, trueRegimeNodes, inplace = FALSE)
  attr(trueModel, "tree") <- trueTree
  attr(trueModel, "X") <- testData$X[[i]]

  plTrue <- plTrue + ggtitle(paste0("TRUE: AIC=", round(AIC(trueModel)), ", lL=", round(logLik(trueModel)), ", p=", PCMParamCount(trueModel, TRUE, TRUE)))

  if(doPlotSearchHistory) {
    plSearchHistory <- PlotSearchHistory(fitMappings)
  } else {
    plSearchHistory <- NULL
  }


  bestFit <- RetrieveBestFitAIC(fitMappings)

  plInferred <- PlotTreeRegimesAndMapping(
    tree,
    regimeNodes = bestFit$inferredRegimeNodes,
    mappingIdx = bestFit$inferredMappingIdx)
  plInferred <- plInferred +
    ggtitle(paste0("FINAL: AIC=", round(AIC(bestFit$inferredModel)),
                   ", lL=", round(logLik(bestFit$inferredModel)),
                   ", p=", PCMParamCount(bestFit$inferredModel, TRUE, TRUE)))

  if(PCMNumTraits(trueModel) == 2 && doPlot2DData) {
    pl2DTrue <- PCMPlotTraitData2D(
      testData$X[[i]][, 1:N],
      trueTree)

    pl2DInferred <- PCMPlotTraitData2D(
      testData$X[[i]][, 1:N],
      PCMTreeSetRegimes(tree, bestFit$inferredRegimeNodes, inplace = FALSE))
  }

  plotlist = c(list(),
    if(doPlotSearchHistory) {
      plSearchHistory[!sapply(plSearchHistory, is.null)]
    } else {
      NULL
    },
    list(inferred=plInferred),
    if(doPlot2DData) {
      list(pl2DInferred = pl2DInferred)
    } else {
      NULL
    },
    list(final=plTrue),
    if(doPlot2DData) {
      list(pl2DTrue = pl2DTrue)
    } else {
      NULL
    })
  plotlist
}

tpr <- function(pred, true) {
  pred <- as.logical(pred)
  true <- as.logical(true)
  res <- sum(pred & true)/sum(true)
  if(is.finite(res)) {
    res
  } else {
    as.double(NA)
  }
}

fpr <- function(pred, true) {
  pred <- as.logical(pred)
  true <- as.logical(true)
  res <- sum(pred & !true)/sum(!true)
  if(is.finite(res)) {
    res
  } else {
    as.double(NA)
  }
}
