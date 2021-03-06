% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv_methods.R
\name{crossValidation}
\alias{crossValidation}
\title{*k*-Fold Cross-Validation}
\usage{
crossValidation(responseData, modelTypes, folds = 10, replications = 1,
  indicator = TRUE, ..., seed = NULL, type = "person")
}
\arguments{
\item{responseData}{The initial, full sample data as a matrix of item responses.}

\item{modelTypes}{A character vector specifying the model types to be compared. Uses the `TAM` package format, so must be one of the following: "1PL", "2PL", "PCM", "PCM2", "RSM", "GPCM", and "2PL.groups".}

\item{folds}{An integer value indicating the number of cross-validation folds to split the data into during the cross-validation process. If specified as `folds=nrow(responseData)` (or some other way equivalently), results in the leave-one-out cross-validation procedure (a.k.a., when *k*=*n*, *k*-fold CV == *n*-fold CV == LOOCV).}

\item{replications}{The number of replications of the *k*-fold CV procedure to perform, with data being split into the *k*-fold groups randomly each time.}

\item{indicator}{A logical value that controls the progress printing.}

\item{...}{Further arguments to be passed to the `tam` function.}

\item{seed}{Either a positive integer setting the random seed, or `NULL`.}

\item{type}{A character vector specifying whether the validation treats the "person" or the "item" as the unit of observation. Default is "person".}
}
\value{
An object of class "cvIRT" with the following values:
\item{seed}{The random seed that produced the results.}
\item{trainData}{A list of each replications' training data.}
\item{testData}{A list of each replications' testing data.}
\item{testLik}{A matrix of the loglikelihood values estimated on the testing data for each model within each holdout replication.}
\item{nModelParams}{A matrix of the number model parameters estimated on the training data within each holdout replication.}
\item{AIC}{A matrix of the AIC value for each holdout replication, as well as the mean value across each replication.}
\item{AICc}{A matrix of the AICc value for each holdout replication, as well as the mean value across each replication.}
\item{BIC}{A matrix of the BIC value for each holdout replication, as well as the mean value across each replication.}
\item{-2 log-Likelihood Ratio Test}{A list of the log-likelihood ratio test statistics, degrees of freedom, and p-values for each model comparison and each replication, as well as the test using the mean test statistics and degrees of freedom. If only one model is used, returns `NULL`.}
\item{warnings}{A character vector of any warnings incurred while the method runs.}
\item{time}{A vector of the start and end times of the function.}
}
\description{
*k*-Fold Cross-Validation
}
\examples{
#None.

}
\concept{cross-validation}
