#' Resubstitution
#'
#' @param responseData The initial, full sample data as a matrix of item responses.
#' @param modelTypes A character vector specifying the model types to be compared. Valid values: "Rasch", "1PL", "2PL", "3PL", "PCM", "PCM2", "RSM", and "GPCM".
#' @param indicator A logical value that controls the progress printing.
#' @param ... Further arguments to be passed to the `mirt` function.
#' @param type A character vector specifying whether the validation treats the "person" or the "item" as the unit of observation. Default is "person".
#'
#' @return An object of class "cvIRT" with the following values:
#' \item{data}{The data used in the function.}
#' \item{testLik}{A matrix of the loglikelihood values estimated on the responseData for each model.}
#' \item{nModelParams}{A matrix of the number model parameters estimated on the data.}
#' \item{AIC}{A matrix of the AIC value for each of the modelTypes.}
#' \item{AICc}{A matrix of the AICc value for each of the modelTypes.}
#' \item{BIC}{A matrix of the BIC value for each of the modelTypes.}
#' \item{-2 log-Likelihood Ratio Test}{A list of the log-likelihood ratio test statistics, degrees of freedom, and p-values for each model comparison. If only one model is used, returns `NULL`.}
#' \item{warnings}{A character vector of any warnings incurred while the method runs.}
#' \item{time}{A vector of the start and end times of the function.}
#'
#' @export
#' @family resubstitution
#'
#' @examples #None.
resubstitution = function(responseData, modelTypes, indicator = TRUE, ..., type = "person") {

  startTime <- Sys.time()
  dots <- list(...)
  dots[c("verbose", "itemtype")] <- NULL

  if (!is.matrix(responseData)&!is.data.frame(responseData)) {
    stop("The responseData needs to be a matrix or data.frame with individuals on the rows and items on the columns.")
  }
  invalid <- setdiff(modelTypes, cvirt_valid_models)
  if (length(invalid) > 0) {
    stop("Unrecognized modelTypes: ", paste(invalid, collapse = ", "),
         ". Must be one of: ", paste(cvirt_valid_models, collapse = ", "))
  }

  if (!is.logical(indicator)) {
    warning("The `indicator` argument was not a logical value. Defaulting to `indicator = TRUE`.",
            immediate. = T)
    indicator <- TRUE
  }

  warningIndicator <- c()

  trainModel <- vector("list", length = length(modelTypes))
  testLik <- matrix(nrow = 1, ncol = length(modelTypes))
  nParamTrain <- matrix(nrow = 1, ncol = length(modelTypes))
  colnames(testLik) = modelTypes

    for (j in 1:length(modelTypes)) {

      if (indicator) cat(paste0("Model: ", modelTypes[j], ".   "))

      result <- do.call(fit_irt_model, c(list(data = responseData, model_type = modelTypes[j], verbose = FALSE), dots))
      trainModel[[j]] <- result

      # extract CV likelihood and number of parameters

      testLik[1,j] <- result$loglik
      nParamTrain[1,j] <- result$npar

    }


  # calculate model selection criteria for each replication and model type
  AICval <- cvAIC(loglikelihood = testLik, numParams = nParamTrain, method = "resub")
  AICCval <- cvAICc(loglikelihood = testLik, numParams = nParamTrain, n = nrow(responseData), method = "resub")
  BICval <- cvBIC(loglikelihood = testLik, numParams = nParamTrain, n = nrow(responseData), method = "resub")
  if (length(modelTypes) > 1){

    LLRtest <- cvLogLikRatio(loglikelihood = testLik, numParams = nParamTrain, models = modelTypes, method = "resub")

  } else {
    LLRtest <- NULL
  }

  endTime <- Sys.time()

  # create the results

  results <- list()
  results$call <- match.call()
  results$data <- responseData
  results$testLik <- testLik
  results$nModelParams <- nParamTrain
  results$AIC <- AICval
  results$AICc <- AICCval
  results$BIC <- BICval
  results$`-2 log-Likelihood Ratio Test` <- LLRtest
  results$warnings <- warningAttr(AICval, warningIndicator)
  results$time <- c(startTime, endTime)

  class(results) <- c("cvIRT", "cvIRTresub")

  return(results)

}
