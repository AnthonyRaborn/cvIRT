#' The Standard Bootstrap
#'
#' @param responseData The initial, full sample data as a matrix of item responses.
#' @param modelTypes A character vector specifying the model types to be compared. Uses the `TAM` package format, so must be one of the following: "1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", and "2PL.groups".
#' @param bootSize An integer value greater than 0 that indicates the number of bootstrap samples to draw.
#' @param replications The number of replications of the bootstrap procedure to perform. Not suggested to use as for the standard bootstrap changing the `bootSize` is more efficient.
#' @param leaveOneOut A logical value indicating whether to use `leaveOneOut` bootstrap. Not suggested for use.
#' @param indicator A logical value that controls the progress printing.
#' @param ... Further arguments to be passed to the `tam` function.
#' @param seed Either a positive integer setting the random seed, or `NULL`.
#' @param type A character vector specifying whether the validation treats the "person" or the "item" as the unit of observation. Default is "person".
#'
#' @return A list object of class "cvIRT" with the following values:
#' \item{call}{The original function calls.}
#' \item{seed}{The random seed that produced the results.}
#' \item{bootstrapSamples}{A list of each bootstrap samples' training data.}
#' \item{testLik}{A matrix of the loglikelihood values estimated on the testing data for each model within each holdout replication.}
#' \item{nModelParams}{A matrix of the number model parameters estimated on the training data within each holdout replication.}
#' \item{AIC}{A matrix of the AIC value for each holdout replication, as well as the mean value across each replication.}
#' \item{AICc}{A matrix of the AICc value for each holdout replication, as well as the mean value across each replication.}
#' \item{BIC}{A matrix of the BIC value for each holdout replication, as well as the mean value across each replication.}
#' \item{-2 log-Likelihood Ratio Test}{A list of the log-likelihood ratio test statistics, degrees of freedom, and p-values for each model comparison and each replication, as well as the test using the mean test statistics and degrees of freedom. If only one model is used, returns `NULL`.}
#' \item{warnings}{A character vector of any warnings incurred while the method runs.}
#' \item{time}{A vector of the start and end times of the function.}
#' @export
#' @family bootstrap
#'
#' @examples #None.
#'
simpleBootstrap <- function(responseData, modelTypes, bootSize = 50, replications = 1, leaveOneOut = F, indicator = TRUE, ..., seed = NULL, type = "person") {

  startTime <- Sys.time()

  if (!is.matrix(responseData)&!is.data.frame(responseData)) {
    stop("The responseData needs to be a matrix or data.frame with individuals on the rows and items on the columns.")
  }
  if (is.null(seed)) {
    seed <- sample(1:1e8, size = 1)
    set.seed(seed)
  } else if (!is.numeric(seed)) {
    warning("The input seed value needs to be NULL or a numeric value! Defaulting to a randomly generated seed value.",
            immediate. = T)
    seed <- sample(1:1e8, size = 1)
    set.seed(seed)
  } else {
    set.seed(seed)
  }

  if (!is.integer(bootSize)&&(bootSize < 0)) {
    warning("The `bootSize` argument needs to be positive! Defaulting to `bootSize = 50`.",
            immediate. = T)
    bootSize = 10
  }

  if (!is.numeric(replications)|replications < 0) {
    warning("The `replications` argument needs to be a positive number! Defaulting to `replications = 1`.",
            immediate. = T)
    replications = 1
  }

  if (!is.logical(indicator)) {
    warning("The `indicator` argument was not a logical value. Defaulting to `indicator = TRUE`.",
            immediate. = T)
    indicator <- TRUE
  }
  warningIndicator <- c()

  bootSamples <- vector('list', length = replications)
  testLik <- matrix(nrow = bootSize, ncol = length(modelTypes))
  nParamTrain <- matrix(nrow = bootSize, ncol = length(modelTypes))
  colnames(testLik) = modelTypes
  testLikList <- replicate(n = replications, testLik, simplify = FALSE)
  nParamTrainList <- replicate(n = replications, nParamTrain, simplify = FALSE)

  AICval <- AICCval <- BICval <- LLRtest <- loobLLRtest <- resubLLRtest <- boot632LLRtest <- vector('list', length = replications)

  if (leaveOneOut) boot632AIC <- boot632AICc <- boot632BIC <- boot632LLRtest <- vector('list', length = replications)

  for (i in 1:replications) {

    testModels <- vector('list', length = bootSize)
    for (b in 1:bootSize) {

      # create the bootstrap samples
      bootSamples[[i]][[b]] <- sample(x = 1:nrow(responseData), size = nrow(responseData), replace = TRUE)

    }

    if (leaveOneOut) {

      leftOutIndex <- vector('list', length = nrow(responseData))

      for (n in 1:nrow(responseData)) {

        for (b in 1:bootSize) {

          leftOutIndex[[n]][b] <-  n %in% unlist(bootSamples[[i]][[b]])
          }
      }

      usableBootSamples <- t(sapply(lapply(leftOutIndex, cbind), cbind))

    }

    for (b in 1:bootSize) {

      testModels[[b]] <- vector('list', length = length(modelTypes))

      for (j in 1:length(modelTypes)) {

        if (indicator) {
        cat(paste0("\r Replication: ", i, " of ", replications,
                   ". Model type: ", modelTypes[j], ". Bootstrap sample: ", b, " of ", bootSize, ".   "))

      }

      # estimate model on bootstrap sample
        if (modelTypes[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
      testModels[[b]][[j]] <- TAM::tam.mml(responseData[bootSamples[[i]][[b]],], irtmodel = modelTypes[j], verbose = F, ...)
        } else if (modelTypes[j] %in% c("2PL", "GPCM", "2PL.groups")) {
          testModels[[b]][[j]] <- TAM::tam.mml.2pl(responseData[bootSamples[[i]][[b]],], irtmodel = modelTypes[j], verbose = F, ...)

        }

        # extract CV likelihood and number of parameters
        nParamTrainList[[i]][b,j] <- testModels[[b]][[j]]$ic$np
        if (!leaveOneOut) {
        testLikList[[i]][b,j] <- testModels[[b]][[j]]$ic$loglike
        }
      }
    }

    if (!leaveOneOut) {
    # calculate model selection criteria for each replication and model type
    AICval[[i]] <- cvAIC(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], method = "bootstrap")
    if (!is.null(attr(AICval[[i]], "warning"))) warningIndicator <- c(warningIndicator, i)
    AICCval[[i]] <- cvAICc(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], n = nrow(responseData), method = "bootstrap")
    BICval[[i]] <- cvBIC(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], n = nrow(responseData), method = "bootstrap")
    if (length(modelTypes) > 1){

      LLRtest[[i]] <- cvLogLikRatio(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], models = modelTypes, method = "bootstrap")

    } else {
      LLRtest[[i]] <- NULL
    }

  } else if (leaveOneOut) {

    loobLogLik <- loobLogLikEst(fittedModels = testModels, loobMatrix = usableBootSamples,
                                models = modelTypes, responses = responseData, bootstrapSamples = bootSamples)
    resubLogLik <- resubLogLikEst(models = modelTypes, responses = responseData)
    boot632logLik <- loob632logLikEst(loobEst = loobLogLik, resubEst = resubLogLik)

    AICval[[i]] <- cvAIC(loglikelihood = loobLogLik, numParams = nParamTrainList[[i]], method = "loob")
    AICCval[[i]] <- cvAICc(loglikelihood = loobLogLik, numParams = nParamTrainList[[i]], n = nrow(responseData), method = "loob")
    BICval[[i]] <- cvBIC(loglikelihood = loobLogLik, numParams = nParamTrainList[[i]], n = nrow(responseData), method = "loob")

    resubAIC[[i]] <- cvAIC(loglikelihood = resubLogLik, numParams = nParamTrainList[[i]], method = "resub")
    resubAICc[[i]] <- cvAICc(loglikelihood = resubLogLik, numParams = nParamTrainList[[i]], n = nrow(responseData), method = "resub")
    resubBIC[[i]] <- cvBIC(loglikelihood = resubLogLik, numParams = nParamTrainList[[i]], n = nrow(responseData), method = "resub")

    boot632AIC[[i]] <- .632*AICval[[i]]+ .368*resubAIC[[i]]
    rownames(boot632AIC[[i]]) <- ".632 Bootstrap AIC:"
    boot632AICc[[i]] <- .632*AICCval[[i]]+ .368*resubAICc[[i]]
    rownames(boot632AICc[[i]]) <- ".632 Bootstrap AICc:"
    boot632BIC[[i]] <- .632*BICval[[i]]+ .368*resubBIC[[i]]
    rownames(boot632BIC[[i]]) <- ".632 Bootstrap BIC:"

    if (length(modelTypes) > 1) {
      loobLLRtest[[i]] <- cvLogLikRatio(loglikelihood = loobLogLik, numParams = nParamTrainList[[i]], models = modelTypes, method = 'loob')
      resubLLRtest[[i]] <- cvLogLikRatio(loglikelihood = resubLogLik, numParams = nParamTrainList[[i]], models = modelTypes, method = 'loob')
      boot632LLRtStatistic <- .632*selectRow(sapply(loobLLRtest[[i]], selectRow, "Leave One Out Bootstrap:"), "Test Statistic") + .368*selectRow(sapply(resubLLRtest[[i]], selectRow, "Leave One Out Bootstrap:"), "Test Statistic")
      boot632LLRtDF <- selectRow(sapply(loobLLRtest[[i]], selectRow, "Leave One Out Bootstrap:"), "degrees of freedom")
      boot632LLRtest[[i]] <- cbind("Test Statistic" = boot632LLRtStatistic, "degrees of freedom" = boot632LLRtDF,  "p-value" = pchisq(boot632LLRtStatistic, df = boot632LLRtDF, lower.tail = F))
    } else {
      loobLLRtest <- resubLLRtest <- boot632LLRtest <- NULL
    }

    }


  }


  if (leaveOneOut) {
    resub <- list("AIC resubstitution" = resubAIC, "AICc resubstitution" = resubAICc, "BIC resubstitution" = resubBIC, "-2 log-Likelihood Ratio Test resubstitution" = resubLLRtest)
  boot632 <- list("AIC .632 Bootstrap" = boot632AIC, "AICc .632 Bootstrap" = boot632AICc, "BIC .632 Bootstrap" = boot632BIC, "-2 log-Likelihood Ratio Test .632 Bootstrap" = boot632LLRtest)
  }
  endTime <- Sys.time()

  # create the results

  results <- list()
  results$call <- match.call()
  results$seed <- seed
  results$bootstrapSamples <- bootSamples
  results$testLik <- testLikList
  results$nModelParams <- nParamTrainList
  results$AIC <- AICval
  results$AICc <- AICCval
  results$BIC <- BICval
  results$`-2 log-Likelihood Ratio Test` <- ifelse(leaveOneOut, loobLLRtest, LLRtest)
  results$warnings <- warningAttr(AICval, warningIndicator)
  results$time <- c(startTime, endTime)
  if (leaveOneOut) {
    results$resubstitution <- resub
    results$`.632 Bootstrap Results` <- boot632
  }

  if (leaveOneOut) {
    class(results) <- c('cvIRT', 'cvIRTloob')
  } else {
  class(results) <- c("cvIRT", "cvIRTbootstrap")
  }

  return(results)

}

#' (Repeated) \emph{k}-Fold Bootstrap
#'
#' @param responseData The initial, full sample data as a matrix of item responses.
#' @param modelTypes A character vector specifying the model types to be compared. Uses the `TAM` package format, so must be one of the following: "1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", and "2PL.groups".
#' @param bootSize An integer value greater than 0 that indicates the number of bootstrap samples to draw.
#' @param folds An integer value indicating the number of cross-validation folds to split the data into during the cross-validation process.
#' @param replications The number of replications of the bootstrap procedure to perform. Not suggested to use as for the standard bootstrap changing the `bootSize` is more efficient.
#' @param leaveOneOut A logical value indicating whether to use `leaveOneOut` bootstrap. Not suggested for use.
#' @param indicator A logical value that controls the progress printing.
#' @param ... Further arguments to be passed to the `tam` function.
#' @param seed Either a positive integer setting the random seed, or `NULL`.
#' @param type A character vector specifying whether the validation treats the "person" or the "item" as the unit of observation. Default is "person".
#'
#' @return A list object of class "cvIRT" with the following values:
#' \item{call}{The original function calls.}
#' \item{seed}{The random seed that produced the results.}
#' \item{bootstrapSamples}{A list of each bootstrap samples' training data.}
#' \item{testLik}{A matrix of the loglikelihood values estimated on the testing data for each model within each holdout replication.}
#' \item{nModelParams}{A matrix of the number model parameters estimated on the training data within each holdout replication.}
#' \item{AIC}{A matrix of the AIC value for each holdout replication, as well as the mean value across each replication.}
#' \item{AICc}{A matrix of the AICc value for each holdout replication, as well as the mean value across each replication.}
#' \item{BIC}{A matrix of the BIC value for each holdout replication, as well as the mean value across each replication.}
#' \item{-2 log-Likelihood Ratio Test}{A list of the log-likelihood ratio test statistics, degrees of freedom, and p-values for each model comparison and each replication, as well as the test using the mean test statistics and degrees of freedom. If only one model is used, returns `NULL`.}
#' \item{warnings}{A character vector of any warnings incurred while the method runs.}
#' \item{time}{A vector of the start and end times of the function.}
#' @export
#' @family bootstrap
#'
#' @examples #None.
#'
kfoldBootstrap <- function(responseData, modelTypes, bootSize = 50, folds = 10, replications = 1, indicator = TRUE, ..., seed = NULL, type = "person") {

  startTime <- Sys.time()

  if (!is.matrix(responseData)&!is.data.frame(responseData)) {
    stop("The responseData needs to be a matrix or data.frame with individuals on the rows and items on the columns.")
  }
  if (is.null(seed)) {
    seed <- sample(1:1e8, size = 1)
    set.seed(seed)
  } else if (!is.numeric(seed)) {
    warning("The input seed value needs to be NULL or a numeric value! Defaulting to a randomly generated seed value.",
            immediate. = T)
    seed <- sample(1:1e8, size = 1)
    set.seed(seed)
  } else {
    set.seed(seed)
  }

  if (!is.integer(bootSize)&&(bootSize < 0)) {
    warning("The `bootSize` argument needs to be positive! Defaulting to `bootSize = 50`.",
            immediate. = T)
    bootSize = 10
  }

  if (!is.numeric(folds)&&(folds < 0 | folds > nrow(responseData))) {
    warning("The `folds` argument needs to be positive and less than the number of observations in the data! Defaulting to `folds = 10`.",
            immediate. = T)
    folds = 10
  }

  if (!is.numeric(replications)|replications < 0) {
    warning("The `replications` argument needs to be a positive number! Defaulting to `replications = 1`.",
            immediate. = T)
    replications = 1
  }

  if (!is.logical(indicator)) {
    warning("The `indicator` argument was not a logical value. Defaulting to `indicator = TRUE`.",
            immediate. = T)
    indicator <- TRUE
  }
  warningIndicator <- c()

  bootSamples <- foldAssignment <- replicate(n = replications, vector('list', length = folds), simplify = FALSE)

  testLik <- matrix(nrow = bootSize, ncol = length(modelTypes))
  nParamTrain <- matrix(nrow = bootSize, ncol = length(modelTypes))
  colnames(testLik) = modelTypes
  testLikList <- replicate(n = replications, replicate(n = folds, testLik, simplify = FALSE), simplify = FALSE)
  nParamTrainList <- replicate(n = replications, replicate(n = folds, nParamTrain, simplify = FALSE), simplify = FALSE)

  AICval <- AICCval <- BICval <- LLRtest <- loobLLRtest <- resubLLRtest <- boot632LLRtest <- vector('list', length = replications)

  for (i in 1:replications) {

    testModels <- vector('list', length = folds)
    foldAssignment[[i]] <- sample(x = rep(1:folds, length.out = nrow(responseData)), size = nrow(responseData))

    for (k in 1:folds){

      for (b in 1:bootSize) {

        # create the bootstrap samples
        bootSamples[[i]][[k]][[b]] <- sample(x = 1:nrow(responseData[foldAssignment[[i]]!=k,]),
                                             size = nrow(responseData[foldAssignment[[i]]!=k,]), replace = TRUE)

      }

    for (b in 1:bootSize) {

      testModels[[k]][[b]] <- vector('list', length = length(modelTypes))

      for (j in 1:length(modelTypes)) {

        if (indicator) {
          cat(paste0("\r Replication: ", i, " of ", replications,
                     ". Fold: ", k, " of ", folds,
                     ". Model type: ", modelTypes[j],
                     ". Bootstrap sample: ", b, " of ", bootSize, ".   "))

        }

        # estimate model on bootstrap sample
        if (modelTypes[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
          testModels[[k]][[b]][[j]] <- TAM::tam.mml(responseData[bootSamples[[i]][[k]][[b]],], irtmodel = modelTypes[j], verbose = F)
        } else if (modelTypes[j] %in% c("2PL", "GPCM", "2PL.groups")) {
          testModels[[k]][[b]][[j]] <- TAM::tam.mml.2pl(responseData[bootSamples[[i]][[k]][[b]],], irtmodel = modelTypes[j], verbose = F)

        }

        # extract CV likelihood and number of parameters
        nParamTrainList[[i]][[k]][b,j] <- testModels[[k]][[b]][[j]]$ic$np
        testLikList[[i]][[k]][b,j] <- testModels[[k]][[b]][[j]]$ic$loglike

      }
    }
    }

      # calculate model selection criteria for each replication and model type
      AICval[[i]] <- cvAIC(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], method = "kfold bootstrap")
      if (!is.null(attr(AICval[[i]], "warning"))) warningIndicator <- c(warningIndicator, i)
      AICCval[[i]] <- cvAICc(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], n = nrow(responseData), method = "kfold bootstrap")
      BICval[[i]] <- cvBIC(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], n = nrow(responseData), method = "kfold bootstrap")
      if (length(modelTypes) > 1){

        LLRtest[[i]] <- cvLogLikRatio(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], models = modelTypes, method = "kfold bootstrap")

      } else {
        LLRtest[[i]] <- NULL
      }

  }

  endTime <- Sys.time()

  # create the results

  results <- list()
  results$call <- match.call()
  results$seed <- seed
  results$bootstrapSamples <- bootSamples
  results$testLik <- testLikList
  results$nModelParams <- nParamTrainList
  results$AIC <- AICval
  results$AICc <- AICCval
  results$BIC <- BICval
  results$`-2 log-Likelihood Ratio Test` <- LLRtest
  results$warnings <- warningAttr(AICval, warningIndicator)
  results$time <- c(startTime, endTime)
  class(results) <- c("cvIRT", "cvIRTbootstrap")

  return(results)

}
