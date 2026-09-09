# Integration tests that exercise end-to-end mirt model fitting.
# These are slower but cover the actual CV/holdout/bootstrap paths.

# Shared test data: 30 rows x 5 items, binary
set.seed(123)
test_data <- matrix(sample(0:1, 150, replace = TRUE), nrow = 30, ncol = 5)
colnames(test_data) <- paste0("I", 1:5)

# --- resubstitution with multiple models ---

test_that("resubstitution works with two models and produces LRT", {
  res <- resubstitution(test_data, modelTypes = c("1PL", "2PL"), indicator = FALSE)

  expect_s3_class(res, "cvIRT")
  expect_s3_class(res, "cvIRTresub")
  expect_equal(ncol(res$testLik), 2)
  expect_false(is.null(res$`-2 log-Likelihood Ratio Test`))
  expect_true(all(is.finite(res$testLik)))

  bm <- bestModel(res)
  expect_s3_class(bm, "cvIRT.bestModels")
})

# --- holdout with single model ---

test_that("holdout runs with 1PL model", {
  res <- holdout(test_data, modelTypes = "1PL", proportion = .3,
                 replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTholdout")
  expect_equal(ncol(res$testLik), 1)
  expect_true(is.finite(res$testLik[1, 1]))
  expect_null(res$`-2 log-Likelihood Ratio Test`)
})

# --- holdout with two models (exercises fixed_pars re-estimation on test data) ---

test_that("holdout with 1PL and 2PL fixes slopes on test data", {
  res <- holdout(test_data, modelTypes = c("1PL", "2PL"), proportion = .3,
                 replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTholdout")
  expect_equal(ncol(res$testLik), 2)
  expect_true(all(is.finite(res$testLik)))
  expect_false(is.null(res$`-2 log-Likelihood Ratio Test`))

  bm <- bestModel(res)
  expect_s3_class(bm, "cvIRT.bestModels")
})

# --- crossValidation with small k ---

test_that("crossValidation 5-fold works with 1PL", {
  res <- crossValidation(test_data, modelTypes = "1PL", folds = 5,
                         replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTkfold")
  expect_equal(length(res$foldAssignment), 1)
  expect_equal(length(res$`-2 log-Likelihood Ratio Test`), 0)
})

# --- crossValidation with two models ---

test_that("crossValidation 5-fold works with two models", {
  res <- crossValidation(test_data, modelTypes = c("1PL", "2PL"), folds = 5,
                         replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTkfold")
  expect_equal(ncol(res$testLik[[1]]), 2)
  expect_false(is.null(res$`-2 log-Likelihood Ratio Test`))

  bm <- bestModel(res)
  expect_s3_class(bm, "cvIRT.bestModels")
})

# --- crossValidation LOOCV (exercises single-row fixed_pars fitting) ---

test_that("crossValidation LOOCV exercises single-row fixed-parameter fits", {
  small_data <- test_data[1:10, ]
  res <- crossValidation(small_data, modelTypes = "1PL",
                         folds = nrow(small_data),
                         indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTkfold")
  expect_equal(nrow(res$testLik[[1]]), nrow(small_data))
})

# --- simpleBootstrap with small bootSize ---

test_that("simpleBootstrap runs with small bootSize", {
  res <- simpleBootstrap(test_data, modelTypes = "1PL", bootSize = 3,
                         replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTbootstrap")
  expect_equal(nrow(res$testLik[[1]]), 3)
  expect_true(all(is.finite(res$testLik[[1]])))
})

test_that("simpleBootstrap with two models produces LRT", {
  res <- simpleBootstrap(test_data, modelTypes = c("1PL", "2PL"), bootSize = 3,
                         replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTbootstrap")
  expect_false(is.null(res$`-2 log-Likelihood Ratio Test`))

  bm <- bestModel(res)
  expect_s3_class(bm, "cvIRT.bestModels")
})

test_that("simpleBootstrap generates a seed when NULL", {
  res <- simpleBootstrap(test_data, modelTypes = "1PL", bootSize = 3,
                         replications = 1, indicator = FALSE, seed = NULL)

  expect_true(is.numeric(res$seed))
  expect_true(res$seed > 0)
})

# --- simpleBootstrap leaveOneOut (.632 bootstrap) ---

test_that("simpleBootstrap leaveOneOut produces finite LOOB and .632 estimates", {
  loob_data <- test_data[1:15, ]
  res <- simpleBootstrap(loob_data, modelTypes = "1PL", bootSize = 20,
                         replications = 1, leaveOneOut = TRUE,
                         indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTloob")
  expect_true(is.finite(res$AIC[[1]]))
  expect_true(is.finite(res$AICc[[1]]))
  expect_true(is.finite(res$BIC[[1]]))

  expect_true(is.finite(res$resubstitution$`AIC resubstitution`[[1]]))
  expect_true(is.finite(res$`.632 Bootstrap Results`$`AIC .632 Bootstrap`[[1]]))
  expect_true(is.finite(res$`.632 Bootstrap Results`$`AICc .632 Bootstrap`[[1]]))
  expect_true(is.finite(res$`.632 Bootstrap Results`$`BIC .632 Bootstrap`[[1]]))
})

test_that("simpleBootstrap leaveOneOut with two models produces LOOB/resub/.632 LRT", {
  loob_data <- test_data[1:15, ]
  res <- simpleBootstrap(loob_data, modelTypes = c("1PL", "2PL"), bootSize = 20,
                         replications = 1, leaveOneOut = TRUE,
                         indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTloob")
  expect_true(all(is.finite(res$AIC[[1]])))
  expect_false(is.null(res$`-2 log-Likelihood Ratio Test`))
  expect_false(is.null(res$resubstitution$`-2 log-Likelihood Ratio Test resubstitution`))
  expect_false(is.null(res$`.632 Bootstrap Results`$`-2 log-Likelihood Ratio Test .632 Bootstrap`))
})

# --- kfoldBootstrap ---

test_that("kfoldBootstrap runs with small bootSize and folds", {
  res <- kfoldBootstrap(test_data, modelTypes = "1PL", bootSize = 3,
                        folds = 3, replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTbootstrap")
  expect_true(length(res$bootstrapSamples) > 0)
})

test_that("kfoldBootstrap with two models produces LRT", {
  res <- kfoldBootstrap(test_data, modelTypes = c("1PL", "2PL"), bootSize = 3,
                        folds = 3, replications = 1, indicator = FALSE, seed = 42)

  expect_s3_class(res, "cvIRTbootstrap")
  expect_false(is.null(res$`-2 log-Likelihood Ratio Test`))

  bm <- bestModel(res)
  expect_s3_class(bm, "cvIRT.bestModels")
})

test_that("kfoldBootstrap generates a seed when NULL", {
  res <- kfoldBootstrap(test_data, modelTypes = "1PL", bootSize = 3,
                        folds = 3, replications = 1, indicator = FALSE, seed = NULL)

  expect_true(is.numeric(res$seed))
  expect_true(res$seed > 0)
})

# --- data.frame input works ---

test_that("crossValidation accepts data.frame input", {
  df <- as.data.frame(test_data)
  res <- crossValidation(df, modelTypes = "1PL", folds = 5,
                         replications = 1, indicator = FALSE, seed = 42)
  expect_s3_class(res, "cvIRTkfold")
})

# --- NULL seed generates reproducible results ---

test_that("NULL seed is generated and stored", {
  res <- holdout(test_data, modelTypes = "1PL", proportion = .3,
                 replications = 1, indicator = FALSE, seed = NULL)
  expect_true(is.numeric(res$seed))
  expect_true(res$seed > 0)
})
