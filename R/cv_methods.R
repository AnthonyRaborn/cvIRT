
#' *k*-Fold Cross-Validation
#'
#' @param responseData The initial, full sample data as a matrix of item responses.
#' @param modelTypes A character vector specifying the model types to be compared. Uses the `TAM` package format, so must be one of the following: "1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", and "2PL.groups".
#' @param folds An integer value indicating the number of cross-validation folds to split the data into during the cross-validation process. If specified as `folds=nrow(responseData)` (or some other way equivalently), results in the leave-one-out cross-validation procedure (a.k.a., when *k*=*n*, *k*-fold CV == *n*-fold CV == LOOCV).
#' @param replications The number of replications of the *k*-fold CV procedure to perform, with data being split into the *k*-fold groups randomly each time.
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
#' @family cross-validation
#'
#' @examples #None.
#'
crossValidation <- function(responseData, modelTypes, folds = 10, replications = 1, indicator = TRUE, ..., seed = NULL, type = "person") {

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
  if (folds==nrow(responseData) && replications!=1){
    warning("Since the number of folds was equal to the number of rows in responseData, `replications` is being forced to equal 1. Check that this was intended.",
            immediate. = T)
    replications = 1
  }

  if (!is.logical(indicator)) {
    warning("The `indicator` argument was not a logical value. Defaulting to `indicator = TRUE`.",
            immediate. = T)
    indicator <- TRUE
  }
  warningIndicator <- c()

  foldAssignment <- vector('list', length = replications)
  testLik <- matrix(nrow = folds, ncol = length(modelTypes))
  nParamTrain <- matrix(nrow = folds, ncol = length(modelTypes))
  colnames(testLik) = modelTypes
  testLikList <- replicate(n = replications, testLik, simplify = FALSE)
  nParamTrainList <- replicate(n = replications, nParamTrain, simplify = FALSE)

  AICval <- AICCval <- BICval <- LLRtest <- vector('list', length = replications)

  if (folds == nrow(responseData)) loocv = T else loocv = F

  for (i in 1:replications) {

    # assign each row to a fold for this replication; notice it stays constant across models j
    if (!loocv){

         foldAssignment[[i]] <- sample(x = rep(1:folds, length.out = nrow(responseData)), size = nrow(responseData))

    } else {

           foldAssignment[[i]] <- 1:folds

    }

    trainModels <- testModels <- vector('list', length = length(modelTypes))

    for (k in 1:folds){

      # for each modelType, create an empty list for each fold
      trainModels <- vector('list', length = length(modelTypes))
      testModels <- vector('list', length = length(modelTypes))

      for (j in 1:length(modelTypes)) {

        if (indicator) cat(paste0("\r Replication: ", i, " of ", replications,
                                  ". Model type: ", modelTypes[j], ". Fold: ", k, " of ", folds, ".   "))


        # separate the data into inSample and outSample prior to model estimation
        inSample <- responseData[foldAssignment[[i]]!=k,]
        outSample <- matrix(responseData[foldAssignment[[i]]==k,], ncol = ncol(responseData))

        if (modelTypes[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
        trainModels[[j]] <- TAM::tam.mml(inSample, irtmodel = modelTypes[j], verbose = F, ...)
        } else if (modelTypes[j] %in% c("2PL", "GPCM", "2PL.groups")) {
          trainModels[[j]] <- TAM::tam.mml.2pl(inSample, irtmodel = modelTypes[j], verbose = F, ...)

        }

        if (loocv) {
          if (rowSums(outSample)==0) {
            # if the person has all 0's, ignore them.
            # will cause problems with the logLike-based methods, but allows for the
            # run to finish
            testModels[[j]] <- NULL
            testModels[[j]]$ic$loglike <- NA
          } else {

            testModels[[j]] <- tam.mml.loocv(outSample, irtmodel = modelTypes[j], maxKiInput = rep(max(responseData), times = ncol(responseData)), xsi.fixed = trainModels[[j]]$xsi.fixed.estimated, xsi.inits = trainModels[[j]]$xsi.fixed.estimated, verbose = F, ...)
          }
        } else {

          if (modelTypes[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
          testModels[[j]] <- tryCatch(TAM::tam.mml(outSample, irtmodel = modelTypes[j], xsi.fixed = trainModels[[j]]$xsi.fixed.estimated, xsi.inits = trainModels[[j]]$xsi.fixed.estimated, verbose = F, ...),
                                      error = function(e) return(NA))
          } else if (modelTypes[j] %in% c("2PL", "GPCM", "2PL.groups")) {
            testModels[[j]] <- tryCatch(TAM::tam.mml.2pl(outSample, irtmodel = modelTypes[j], xsi.fixed = trainModels[[j]]$xsi.fixed.estimated, xsi.inits = trainModels[[j]]$xsi.fixed.estimated, verbose = F, ...),
                                        error = function(e) return(NA))

          }
        }
        if (is.na(testModels[[j]])) {
          testLikList[[i]][k,j] = NA
          nParamTrainList[[i]][k,j] = trainModels[[j]]$ic$np
        } else {
        # extract CV likelihood and number of parameters
        testLikList[[i]][k,j] <- testModels[[j]]$ic$loglike
        nParamTrainList[[i]][k,j] <- trainModels[[j]]$ic$np
        }

      }
    }

    # calculate model selection criteria for each replication and model type
    AICval[[i]] <- cvAIC(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], method = "kfold")
    if (!is.null(attr(AICval[[i]], "warning"))) warningIndicator <- c(warningIndicator, i)
    AICCval[[i]] <- cvAICc(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], n = (1-1/folds)*nrow(responseData), method = "kfold")
    BICval[[i]] <- cvBIC(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], n = (1-1/folds)*nrow(responseData), method = "kfold")
    if (length(modelTypes) > 1){

      LLRtest[[i]] <- cvLogLikRatio(loglikelihood = testLikList[[i]], numParams = nParamTrainList[[i]], models = modelTypes, method = "kfold")

    } else {
      LLRtest[[i]] <- NULL
    }


  }

  endTime <- Sys.time()

  # create the results

  results <- list()
  results$call <- match.call()
  results$seed <- seed
  results$foldAssignment <- foldAssignment
  results$testLik <- testLikList
  results$nModelParams <- nParamTrainList
  results$AIC <- AICval
  results$AICc <- AICCval
  results$BIC <- BICval
  results$`-2 log-Likelihood Ratio Test` <- LLRtest
  results$warnings <- warningAttr(AICval, warningIndicator)
  results$time <- c(startTime, endTime)

  class(results) <- c("cvIRT", "cvIRTkfold")

  return(results)



}
