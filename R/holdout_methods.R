#' Holdout Validation
#'
#' @param responseData The initial, full sample data as a matrix of item responses.
#' @param modelTypes A character vector specifying the model types to be compared. Uses the `TAM` package format, so must be one of the following: "1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", and "2PL.groups".
#' @param proportion The proportion of the data to hold out for validation. Must be a numeric value between 0 and 1, exclusive.
#' @param replications The number of times to replicate the holdout process using the same data, models, and proportion.
#' @param indicator A logical value that controls the progress printing.
#' @param ... Further arguments to be passed to the `tam` function.
#' @param seed Either a positive integer setting the random seed, or `NULL`.
#' @param type A character vector specifying whether the validation treats the "person" or the "item" as the unit of observation. Default is "person".
#'
#' @return An object of class "cvIRT" with the following values:
#' \item{seed}{The random seed that produced the results.}
#' \item{trainData}{A list of each replications' training data.}
#' \item{testData}{A list of each replications' testing data.}
#' \item{testLik}{A matrix of the loglikelihood values estimated on the testing data for each model within each holdout replication.}
#' \item{nModelParams}{A matrix of the number model parameters estimated on the training data within each holdout replication.}
#' \item{AIC}{A matrix of the AIC value for each holdout replication, as well as the mean value across each replication.}
#' \item{AICc}{A matrix of the AICc value for each holdout replication, as well as the mean value across each replication.}
#' \item{BIC}{A matrix of the BIC value for each holdout replication, as well as the mean value across each replication.}
#' \item{-2 log-Likelihood Ratio Test}{A list of the log-likelihood ratio test statistics, degrees of freedom, and p-values for each model comparison and each replication, as well as the test using the mean test statistics and degrees of freedom. If only one model is used, returns `NULL`.}
#' \item{warnings}{A character vector of any warnings incurred while the method runs.}
#' \item{time}{A vector of the start and end times of the function.}
#'
#' @export
#' @family holdout
#' @importFrom TAM tam.mml tam.mml.2pl
#'
#' @examples #None.
holdout = function(responseData, modelTypes, proportion, replications, indicator = TRUE, ..., seed = NULL, type = "person") {

  startTime <- Sys.time()

  if (!is.matrix(responseData)&!is.data.frame(responseData)) {
    stop("The responseData needs to be a matrix or data.frame with individuals on the rows and items on the columns.",
         immediate. = T)
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

  if (!is.numeric(proportion)||!(0 < proportion && proportion <1)) {
    warning("The input proportion value needs to be a numeric value between 0 and 1. Defaulting to `proportion=.33`",
            immediate. = T)
    proportion = .33
  }

  if (!is.numeric(replications) | replications < 0) {
    warning("The `replications` argument needs to be a positive integer! Defaulting to `replications = 1`.",
            immediate. = T)
    replications = 1
  }

  if (!is.logical(indicator)) {
    warning("The `indicator` argument was not a logical value. Defaulting to `indicator = TRUE`.",
            immediate. = T)
    indicator <- TRUE
  }

  warningIndicator <- c()


  testData <- trainData <- vector("list", length = replications)
  trainModel <- testModel <- c(rep(list(vector("list", length = length(modelTypes))), times = replications))
  testLik <- matrix(nrow = replications, ncol = length(modelTypes))
  nParamTrain <- matrix(nrow = replications, ncol = length(modelTypes))
  colnames(testLik) = modelTypes

  for (i in 1:replications) {

    # split data into training and testing sets

    testSampleRows <- sample(1:nrow(responseData), size = proportion*nrow(responseData), replace = F)
    testData[[i]] <- responseData[testSampleRows,]
    trainData[[i]] <- responseData[-testSampleRows,]

    for (j in 1:length(modelTypes)) {

      if (indicator) cat(paste0("\r Replication: ", i, " of ", replications, ". Model: ", modelTypes[j], ".   "))

      if (modelTypes[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
        # estimate each model on training set
    trainModel[[i]][[j]] <- TAM::tam.mml(trainData[[i]], irtmodel = modelTypes[j], verbose = F, ...)

    # estimate each model on testing set

    testModel[[i]][[j]] <- TAM::tam.mml(testData[[i]], irtmodel = modelTypes[j], xsi.fixed = trainModel[[i]][[j]]$xsi.fixed.estimated, verbose = F, ...)
      } else if (modelTypes[j] %in% c("2PL", "GPCM", "2PL.groups")) {
        # estimate each model on training set
          trainModel[[i]][[j]] <- TAM::tam.mml.2pl(trainData[[i]], irtmodel = modelTypes[j], verbose = F, ...)

          # estimate each model on testing set

          testModel[[i]][[j]] <- TAM::tam.mml.2pl(testData[[i]], irtmodel = modelTypes[j], xsi.fixed = trainModel[[i]][[j]]$xsi.fixed.estimated, verbose = F, ...)

      }

    # extract CV likelihood and number of parameters

    testLik[i,j] <- testModel[[i]][[j]]$ic$loglike
    nParamTrain[i,j] <- trainModel[[i]][[j]]$ic$np

    }
  }

  # calculate model selection criteria for each replication and model type
  AICval <- cvAIC(loglikelihood = testLik, numParams = nParamTrain, method = "holdout")
  AICCval <- cvAICc(loglikelihood = testLik, numParams = nParamTrain, n = proportion*nrow(responseData), method = "holdout")
  BICval <- cvBIC(loglikelihood = testLik, numParams = nParamTrain, n = proportion*nrow(responseData), method = "holdout")
  if (length(modelTypes) > 1){

    LLRtest <- cvLogLikRatio(loglikelihood = testLik, numParams = nParamTrain, models = modelTypes, method = "holdout")

  } else {
    LLRtest <- NULL
  }

  endTime <- Sys.time()

  # create the results

  results <- list()
  results$seed <- seed
  results$trainData <- trainData
  results$testData <- testData
  results$testLik <- testLik
  results$nModelParams <- nParamTrain
  results$AIC <- AICval
  results$AICc <- AICCval
  results$BIC <- BICval
  results$`-2 log-Likelihood Ratio Test` <- LLRtest
  results$warnings <- warningAttr(AICval, warningIndicator)
  results$time <- c(startTime, endTime)

  class(results) <- c("cvIRT", "cvIRTholdout")

  return(results)

  }
