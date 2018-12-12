cvAIC <- function(loglikelihood, numParams, method) {

  if (method == "loob"|method=="resub") {
    AIC = -2*loglikelihood[nrow(loglikelihood),] + 2*numParams[1,]
    AIC = matrix(AIC, ncol = ncol(loglikelihood))
    if (method == "loob") rownames(AIC) = "Leave One Out Bootstrap:" else rownames(AIC) = "Resubstitution: "
    colnames(AIC) = colnames(loglikelihood)
    attr(AIC, "warning") <- NULL
    return(AIC)
  } else if (method=="kfold bootstrap") {
    AIC = vector('list', length = length(loglikelihood))
    meanAIC = matrix(nrow = length(loglikelihood), ncol = ncol(loglikelihood[[1]]))
    for (i in 1:length(loglikelihood)) {
      AIC[[i]] = -2*loglikelihood[[i]] + 2*numParams[[i]]
      rownames(AIC[[i]]) = paste0("Fold ", i, " Boostrap Sample ", 1:nrow(AIC[[i]]), ":")
      meanAIC[i,] = colMeans(AIC[[i]])
    }
    colnames(meanAIC) = colnames(loglikelihood[[1]])
    rownames(meanAIC) = paste0("Fold ", 1:nrow(meanAIC), " Mean AIC:")
    final = rbind(meanAIC, "Mean:" = colMeans(meanAIC))
    return(final)
    } else {
  AIC = -2*loglikelihood + 2*numParams
  }

  if (method == "resubstitution") {
    rownames(AIC) = "Resubstitution AIC:"
    return(AIC)
  } else if (method == "holdout") {
  rownames(AIC) = paste0("Replication ", 1:nrow(AIC), ":")
  } else if (method == "kfold") {
    rownames(AIC) = paste0("Fold ", 1:nrow(AIC), ":")
  } else if (method == "bootstrap") {
    rownames(AIC) = paste0("Bootstrap Sample ", 1:nrow(AIC), ":")
  }

  warn <- NULL

  meanAIC = colMeans(AIC)
  if(is.na(sum(meanAIC))) {
    meanAIC = colMeans(AIC, na.rm = T)
    warn <- warning("There were missing log-likelihood values! Please check your data, model, and assumptions.\nThe problematic model type(s): \"", paste(names(which(is.na(colMeans(AIC)))), collapse = "\", \""), "\".\nThe missing log-likelihoods were removed from the AIC, AICc, BIC, and -2*log-likelihood test mean values.\nThese may not be reliable!",
            immediate. = F)
  }
  final <- rbind(AIC, "Mean:" = meanAIC)
  attr(final, "warning") <- warn

  return(final)
}

cvAICc <- function(loglikelihood, numParams, n, method) {
  if (method == "loob"|method=="resub") {
    AICc = -2*loglikelihood[nrow(loglikelihood),] + 2*numParams[1,] + (2*(numParams[1,]+2)*(numParams[1,]+1))/(n - numParams[1,]-2)
    AICc = matrix(AICc, ncol = ncol(loglikelihood))
    if (method=="loob") rownames(AICc) = "Leave One Out Bootstrap:" else rownames(AICc) = "Resubstitution:"
    colnames(AICc) = colnames(loglikelihood)
    return(AICc)
  } else if (method=="kfold bootstrap") {
    AICc = vector('list', length = length(loglikelihood))
    meanAICc = matrix(nrow = length(loglikelihood), ncol = ncol(loglikelihood[[1]]))
    for (i in 1:length(loglikelihood)) {
      AICc[[i]] = -2*loglikelihood[[i]] + 2*numParams[[i]] +
        (2*(numParams[[i]]+2)*(numParams[[i]]+1))/(n - numParams[[i]] -2)
      rownames(AICc[[i]]) = paste0("Fold ", i, " Boostrap Sample ", 1:nrow(AICc[[i]]), ":")
      meanAICc[i,] = colMeans(AICc[[i]])
    }
    colnames(meanAICc) = colnames(loglikelihood[[1]])
    rownames(meanAICc) = paste0("Fold ", 1:nrow(meanAICc), " Mean AICc:")
    final = rbind(meanAICc, "Mean:" = colMeans(meanAICc))
    return(final)
  } else {
    AICc = -2*loglikelihood + 2*numParams + (2*(numParams+2)*(numParams+1))/(n - numParams -2)
  }

  if (method == "resubstitution") {
    rownames(AICc) = "Resubstitution AICc:"
    return(AICc)
  } else if (method == "holdout") {
    rownames(AICc) = paste0("Replication ", 1:nrow(AICc), ":")
  } else if (method == "kfold") {
    rownames(AICc) = paste0("Fold ", 1:nrow(AICc), ":")
  } else if (method == "bootstrap") {
    rownames(AICc) = paste0("Bootstrap Sample ", 1:nrow(AICc), ":")
    }

  meanAICc = colMeans(AICc)
  if (is.na(sum(meanAICc))) {
    meanAICc = colMeans(AICc, na.rm = T)
  }
  final <- rbind(AICc, "Mean:" = meanAICc)

  return(final)
}

cvBIC <- function(loglikelihood, numParams, n, method) {
  if (method == "loob"|method=="resub") {
    BIC = -2*loglikelihood[nrow(loglikelihood),] + numParams[1,]*log(n)
    BIC = matrix(BIC, ncol = ncol(loglikelihood))
    if (method=="loob") rownames(BIC) = "Leave One Out Bootstrap:" else rownames(BIC) = "Resubstitution"
    colnames(BIC) = colnames(loglikelihood)
    return(BIC)
  } else if (method=="kfold bootstrap") {
    BIC = vector('list', length = length(loglikelihood))
    meanBIC = matrix(nrow = length(loglikelihood), ncol = ncol(loglikelihood[[1]]))
    for (i in 1:length(loglikelihood)) {
      BIC[[i]] = -2*loglikelihood[[i]] + numParams[[i]]*log(n)
      rownames(BIC[[i]]) = paste0("Fold ", i, " Boostrap Sample ", 1:nrow(BIC[[i]]), ":")
      meanBIC[i,] = colMeans(BIC[[i]])
    }
    colnames(meanBIC) = colnames(loglikelihood[[1]])
    rownames(meanBIC) = paste0("Fold ", 1:nrow(meanBIC), " Mean BIC:")
    final = rbind(meanBIC, "Mean:" = colMeans(meanBIC))
    return(final)
  } else {
    BIC = -2*loglikelihood + numParams*log(n)
  }

  if (method == "resubstitution") {
    rownames(BIC) = "Resubstitution BIC:"
    return(BIC)
  } else if (method == "holdout") {
    rownames(BIC) = paste0("Replication ", 1:nrow(BIC), ":")
  } else if (method == "kfold") {
    rownames(BIC) = paste0("Fold ", 1:nrow(BIC), ":")
  } else if (method == "bootstrap") {
    rownames(BIC) = paste0("Bootstrap Sample ", 1:nrow(BIC), ":")
  }

  meanBIC = colMeans(BIC)
  if (is.na(sum(meanBIC))) {
    meanBIC = colMeans(BIC, na.rm = T)
  }
  final <- rbind(BIC, "Mean:" = meanBIC)

  return(final)
}

cvLogLikRatio <- function(loglikelihood, numParams, models, method) {
  if (method %in% c("loob", "resub")) {
    testStat = degFree = matrix(nrow = 1, ncol = ncol(loglikelihood) - 1)
    modelsTested = vector(length = ncol(loglikelihood)-1)
    for (i in 2:ncol(loglikelihood)) {
    testStat[,i-1] = abs(2*(loglikelihood[nrow(loglikelihood),i-1] - loglikelihood[nrow(loglikelihood),i]))
    degFree[,i-1] = abs(numParams[1,i] - numParams[1,i-1])
    modelsTested[i-1] = paste0(models[i-1], " vs. ", models[i])
    }
    pval = stats::pchisq(testStat, df = degFree, lower.tail = F)
    results <- vector("list", length = length(modelsTested))
    for (i in 1:length(modelsTested)) {
      temp = cbind(testStat[,i], degFree[,i], pval[,i])
      colnames(temp) = c("Test Statistic", "degrees of freedom", "p-value")
      if (method=="loob") rownames(temp) = "Leave One Out Bootstrap:" else rownames(temp) = "Resubstitution:"
      results[[i]] <- temp
    }
    names(results) = modelsTested
    return(results)
  } else if (method == "kfold bootstrap") {

    results <- replicate(n = length(models),
                          matrix(nrow = length(loglikelihood),
                                 ncol = 3),
                          simplify = F)

    for (j in 1:length(loglikelihood)){

      testStat = degFree = matrix(nrow = 1, ncol = ncol(loglikelihood[[j]]) - 1)
      modelsTested = vector(length = ncol(loglikelihood[[j]])-1)

      for (i in 2:ncol(loglikelihood[[j]])) {

        testStat[,i-1] = abs(2*(loglikelihood[[j]][nrow(loglikelihood[[j]]),i-1] -
                                  loglikelihood[[j]][nrow(loglikelihood[[j]]),i]))
        degFree[,i-1] = abs(numParams[[j]][1,i] - numParams[[j]][1,i-1])
        modelsTested[i-1] = paste0(models[i-1], " vs. ", models[i])

      }

      pval = stats::pchisq(testStat, df = degFree, lower.tail = F)

      for (i in 1:length(modelsTested)) {
        temp = cbind(testStat[,i], degFree[,i], pval[,i])
        colnames(temp) = colnames(results[[i]]) = c("Test Statistic", "degrees of freedom", "p-value")
        rownames(results[[i]]) = paste0("Fold ", 1:nrow(results[[i]]), " Bootstrap:")
        results[[i]][j,] <- temp
      }
      names(results) = modelsTested
    }

    for (i in 1:length(results)) {
      meanTS = mean(results[[i]][,1])
      meanDF = mean(results[[i]][,2])
      meanP = stats::pchisq(meanTS, meanDF, lower.tail = F)
      results[[i]] = rbind(results[[i]], "Mean:" = c(meanTS, meanDF, meanP))
    }

    return(results)

  }
  testStat = degFree = matrix(nrow = nrow(loglikelihood), ncol = ncol(loglikelihood) - 1)
  modelsTested = vector(length = ncol(loglikelihood)-1)

  if (numParams[1,2] - numParams[1,1] < 0){

  for (i in 2:ncol(loglikelihood)){
    testStat[,i-1] = 2*(loglikelihood[,i-1] - loglikelihood[,i])
    degFree[,i-1] = numParams[,i] - numParams[,i-1]
    modelsTested[i-1] = paste0(models[i-1], " vs. ", models[i])

  }

    pval = stats::pchisq(testStat, df = numParams[,i-1] - numParams[,i], lower.tail = F)

    } else if (numParams[1,2] - numParams[1,1] > 0) {

      for (i in 2:ncol(loglikelihood)){

        testStat[,i-1] = 2*(loglikelihood[,i] - loglikelihood[,i-1])
        degFree[,i-1] = numParams[,i] - numParams[,i-1]
        modelsTested[i-1] = paste0(models[i-1], " vs. ", models[i])

      }

    pval = stats::pchisq(testStat, df = numParams[,i] - numParams[,i-1], lower.tail = F)

    }


  results <- vector("list", length = length(modelsTested))
  for (i in 1:length(modelsTested)) {
    temp = cbind(testStat[,i], degFree[,i], pval[,i])
    colnames(temp) = c("Test Statistic", "degrees of freedom", "p-value")

    if (method == "resubstitution") {
      rownames(temp) = "Resubstitution:"
    } else if (method == "holdout") {
      rownames(temp) = paste0("Replication ", 1:nrow(loglikelihood), ":")
    } else if (method == "kfold") {
      rownames(temp) = paste0("Fold ", 1:nrow(loglikelihood), ":")
    } else if (method == "bootstrap") {
      rownames(temp) = paste0("Bootstrap Sample ", 1:nrow(loglikelihood), ":")
    }

    temp = rbind(temp, "Mean:" = c(mean(testStat[,i], na.rm = T), mean(degFree[,i], na.rm = T), stats::pchisq(mean(testStat[,i]), df = mean(degFree[,i]), lower.tail = F)))
    results[[i]] = temp
  }
  names(results) <-  modelsTested
  return(results)
}

warningAttr <- function(x, warningIndicator) {
  warn <- c()
  for (i in warningIndicator) {
    warn <- c(warn, attr(x[[i]], 'warning'))
  }
  return(warn)
}

selectRow <- function(x, name) {
  # uses output from class "cvList"
  if (!is.list(x)) {
    index <- which(rownames(x) == name)
    return(x[index,])
  } else if (is.list(x)) {
    temp <- vector('list', length = length(x))
    for (i in 1:length(x)){

      index <- which(rownames(x[[i]]) == name)
      temp[[i]] <- x[[i]][index,]
    }
    names(temp) = paste("Mean", names(x))
    return(temp)

  }
}

selectNamedVector <- function(x, name) {
  # uses output from class "cvList"
  if (!is.list(x)) {
    index <- which(names(x) == name)
    return(x[index])
  } else if (is.list(x)) {
    temp <- vector('list', length = length(x))
    for (i in 1:length(x)){

      index <- which(names(x[[i]]) == name)
      temp[[i]] <- x[[i]][index]
    }
    names(temp) = paste(name, names(x))
    return(temp)

  }
}

selectListElement <- function(x) {
  # uses output from class "cvList"
  # and from the lapply function call for methods "print", "summary"
  elementNames <- names(x[[1]])

  temp <- vector('list', length = length(elementNames))
  names(temp) <- elementNames
  for (i in 1:length(elementNames)){
    for (j in 1:length(x)) {
      temp[[i]] <- rbind(temp[[i]], x[[j]][[i]])
    }
  }
  return(temp)
}

bestModel <- function(results, alpha = .05) {
  if ( !('cvIRT' %in% class(results)) ) {
    stop("The object is not of class 'cvIRT'!")
  }

  if ("cvIRTholdout" %in% class(results)) {
    AIC = results$AIC[2, which.min(results$AIC[2,])]
    names(AIC) = names(which.min(results$AIC[2,]))

    BIC = results$BIC[2, which.min(results$BIC[2,])]
    names(BIC) = names(which.min(results$BIC[2,]))

    AICc = results$AICc[2, which.min(results$AICc[2,])]
    names(AICc) = names(which.min(results$AICc[2,]))

    logLikResults <- lapply(results$`-2 log-Likelihood Ratio Test`, selectRow, "Mean:")

  } else if ("cvIRTkfold" %in% class(results) | "cvIRTbootstrap" %in% class(results)) {
    meanAIC = colMeans(do.call(rbind, lapply(results$AIC, selectRow, "Mean:")))
    AIC = meanAIC[which.min(meanAIC)]
    names(AIC) = names(which.min(meanAIC))

    meanBIC = colMeans(do.call(rbind, lapply(results$BIC, selectRow, "Mean:")))
    BIC = meanBIC[which.min(meanBIC)]
    names(BIC) = names(which.min(meanBIC))

    meanAICc = colMeans(do.call(rbind, lapply(results$AICc, selectRow, "Mean:")))
    AICc = meanAICc[which.min(meanAICc)]
    names(AICc) = names(which.min(meanAICc))

    logLikResults <- lapply(selectListElement(lapply(results$`-2 log-Likelihood Ratio Test`, selectRow, "Mean:")), colMeans)

  } else if ("cvIRTresub" %in% class(results)) {
    AIC = results$AIC[1, which.min(results$AIC[1,])]
    names(AIC) = names(which.min(results$AIC[1,]))

    BIC = results$BIC[1, which.min(results$BIC[1,])]
    names(BIC) = names(which.min(results$BIC[1,]))

    AICc = results$AICc[1, which.min(results$AICc[1,])]
    names(AICc) = names(which.min(results$AICc[1,]))

    logLikResults <- lapply(results$`-2 log-Likelihood Ratio Test`, selectRow, "Resubstitution:")

  }

  # Use the Holm-Bonferroni adjustment for multiple comparisons automatically
  logLikpval <- selectListElement(logLikResults)$`p-value`
  logLikResultsOrder <- logLikResults[order(logLikpval)]
  alphaHB <- alpha/length(logLikResults):1
  if (is.na(names(logLikResultsOrder[max(which(logLikpval[order(logLikpval),] < alphaHB))]))) {
    logLikBestResult <- logLikResults[1]
    logLikBest <- regmatches(names(logLikBestResult),
                             regexpr("[[:alnum:]]*(?= vs.)",
                                     names(logLikBestResult),
                                     perl = T))
  } else {
    logLikBestResult <- logLikResults[max(which(logLikpval < alpha))]
    logLikBest <- regmatches(names(logLikBestResult),
                           regexpr("(?<=vs. )[[:alnum:]]*",
                                   names(logLikBestResult),
                                   perl = T))
  }
  logLikBestResult[[1]] = c(logLikBestResult[[1]], "Best Model" = logLikBest)

  selectedModelInfo = c(AIC, BIC, AICc)
  selectedModelLRT = logLikBestResult

  final = list(info = selectedModelInfo, lrt = selectedModelLRT)
  class(final) = "cvIRT.bestModels"
  return(final)

}

print.cvIRT.bestModels <- function(x) {
  cat("The best model for each selection criteria:")
  cat(paste0("\nAIC: ", names(x$info[1]), ".\tValue: ", signif(x$info[1]), 3))
  cat(paste0("\nBIC: ", names(x$info[2]), ".\tValue: ", signif(x$info[2]), 3))
  cat(paste0("\nAICc: ", names(x$info[3]), ".\tValue: ", signif(x$info[3]), 3))
  cat(paste0("\nLRT: ", x$lrt[[1]][4], ".\tp-value: ", signif(as.numeric(x$lrt[[1]][3]),3), ".\tModels Compared: ", names(x$lrt)))
}

extract.cvIRT.bestModels <- function(x) {
  AIC <- names(x$info[1])
  BIC <- names(x$info[2])
  AICc <- names(x$info[3])
  LRT <- as.character(x$lrt[[1]][4])

  final <- c(AIC = AIC, BIC = BIC, AICc = AICc, LRT = LRT)

  return(final)
  }

### Removed ####
# loobLogLikEst <- function(fittedModels, loobMatrix, models, responses, bootstrapSamples) {
#   loobLogLik <- loobLogLikMean <- vector('list', length = length(models))
#   for (j in 1:length(models)) {
#     loobLogLik[[j]] <- matrix(nrow = nrow(loobMatrix), ncol = ncol(loobMatrix))
#     for (n in 1:nrow(loobMatrix)) {
#       for (b in 1:length(fittedModels)) {
#         cat(paste0("\rFitting observation ", n, " of ", nrow(loobMatrix), " with bootstrap sample ", b, " of ", length(fittedModels), " and the ", models[j], " model (", j, " of ", length(models), " models).      \t\t"))
#         if (loobMatrix[n,b]) {
#           if (models[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
#             tempLogLik <- tam.mml.loocv(matrix(responses[n,], ncol = ncol(responses)), irtmodel = models[j], xsi.fixed = fittedModels[[b]][[j]]$xsi.fixed.estimated, xsi.inits = fittedModels[[b]][[j]]$xsi.fixed.estimated, maxKiInput = rep(max(responses), times = ncol(responses)), verbose = F)$ic$loglike
#           } else if (models[j] %in% c("2PL", "GPCM", "2PL.groups")) {
#             tempLogLik <- tam.mml.2pl.loocv(matrix(responses[n,], ncol = ncol(responses)), irtmodel = models[j], xsi.fixed = fittedModels[[b]][[j]]$xsi.fixed.estimated, xsi.inits = fittedModels[[b]][[j]]$xsi.fixed.estimated, maxKiInput = rep(max(responses), times = ncol(responses)), verbose = F)$ic$loglike
#           }
#           loobLogLik[[j]][n,b] <- tempLogLik
#         }
#       }
#     }
#   }
#   loobLogLik <- lapply(loobLogLik, rowMeans, na.rm = T)
#   loobLogLikMean <- sapply(loobLogLik, cbind)
#   rownames(loobLogLikMean) <- paste0("Observation ", 1:nrow(loobLogLikMean), ":")
#   colnames(loobLogLikMean) <- models
#   loobLogLikMeanFinal <- rbind(loobLogLikMean, "Leave-One-Out Mean:" = colMeans(loobLogLikMean))
#   return(loobLogLikMeanFinal)
# }
# resubLogLikEst <- function(models, responses) {
#   logLik <- logLikMean <- vector('list', length = length(models))
#   for (j in 1:length(models)) {
#     if (models[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
#       tempModel <- TAM::tam.mml(responses, irtmodel = models[j], verbose = F)
#     } else if (models[j] %in% c("2PL", "GPCM", "2PL.groups")) {
#       tempModel <- TAM::tam.mml.2pl(responses, irtmodel = models[j], verbose = F)
#     }
#     for (n in 1:nrow(responses)) {
#       cat(paste0("\rFitting observation ", n, " of ", nrow(responses), " with the ", models[j], " model (", j, " of ", length(models), " models) for resubstitution.     \t\t\t "))
#       if (models[j] %in% c("1PL", "PCM", "PCM2", "RSM")) {
#         tempLik <- tam.mml.loocv(matrix(responses[n,], ncol = ncol(responses)), irtmodel = models[j], xsi.fixed = tempModel$xsi.fixed.estimated, xsi.inits = tempModel$xsi.fixed.estimated, maxKiInput = rep(max(responses), times = ncol(responses)), verbose = F)$ic$loglike
#       } else if (models[j] %in% c("2PL", "GPCM", "2PL.groups")) {
#         tempLik <- tam.mml.2pl.loocv(matrix(responses[n,], ncol = ncol(responses)), irtmodel = models[j], xsi.fixed = tempModel$xsi.fixed.estimated, xsi.inits = tempModel$xsi.fixed.estimated, maxKiInput = rep(max(responses), times = ncol(responses)), verbose = F)$ic$loglike
#       }
#       logLik[[j]][n] <- tempLik
#       names(logLik) <- models
#     }
#   }
#   resubLogLikMean <- sapply(logLik, cbind)
#   rownames(resubLogLikMean) <- paste0("Observation ", 1:nrow(resubLogLikMean), ":")
#   resubLogLikMeanFinal <- rbind(resubLogLikMean, "Resubstitution Mean:" = colMeans(resubLogLikMean))
#   return(resubLogLikMeanFinal)
# }
#
# loob632logLikEst <- function(loobEst, resubEst) {
#   results = t(as.matrix(.632*loobEst[nrow(loobEst),] - .368*resubEst[nrow(resubEst),]))
#   rownames(results) = ".632 Bootstrap log-Likelihood:"
#   return(results)
# }

