# Tests for bestModel, print.cvIRT.bestModels, extract.cvIRT.bestModels

bestModel <- cvIRT:::bestModel
extract.cvIRT.bestModels <- cvIRT:::extract.cvIRT.bestModels

# Helper to build a minimal holdout-style result object
make_holdout_result <- function(aic_vals = c(200, 190),
                                bic_vals = c(210, 205),
                                aicc_vals = c(201, 191),
                                lrt_pval = 0.01) {
  models <- c("1PL", "2PL")
  aic <- matrix(aic_vals, nrow = 1, ncol = 2, dimnames = list(NULL, models))
  aic <- rbind(aic, "Mean:" = aic_vals)
  bic <- matrix(bic_vals, nrow = 1, ncol = 2, dimnames = list(NULL, models))
  bic <- rbind(bic, "Mean:" = bic_vals)
  aicc <- matrix(aicc_vals, nrow = 1, ncol = 2, dimnames = list(NULL, models))
  aicc <- rbind(aicc, "Mean:" = aicc_vals)

  lrt_mat <- matrix(c(15.0, 5.0, lrt_pval), nrow = 1, ncol = 3)
  colnames(lrt_mat) <- c("Test Statistic", "degrees of freedom", "p-value")
  rownames(lrt_mat) <- "Mean:"
  lrt <- list("1PL vs. 2PL" = lrt_mat)

  results <- list(
    AIC = aic, BIC = bic, AICc = aicc,
    `-2 log-Likelihood Ratio Test` = lrt
  )
  class(results) <- c("cvIRT", "cvIRTholdout")
  results
}

make_resub_result <- function(aic_vals = c(200, 190),
                              bic_vals = c(210, 205),
                              aicc_vals = c(201, 191),
                              lrt_pval = 0.01) {
  models <- c("1PL", "2PL")
  aic <- matrix(aic_vals, nrow = 1, ncol = 2, dimnames = list("Resubstitution: ", models))
  bic <- matrix(bic_vals, nrow = 1, ncol = 2, dimnames = list("Resubstitution", models))
  aicc <- matrix(aicc_vals, nrow = 1, ncol = 2, dimnames = list("Resubstitution:", models))

  lrt_mat <- matrix(c(15.0, 5.0, lrt_pval), nrow = 1, ncol = 3)
  colnames(lrt_mat) <- c("Test Statistic", "degrees of freedom", "p-value")
  rownames(lrt_mat) <- "Resubstitution:"
  lrt <- list("1PL vs. 2PL" = lrt_mat)

  results <- list(
    AIC = aic, BIC = bic, AICc = aicc,
    `-2 log-Likelihood Ratio Test` = lrt
  )
  class(results) <- c("cvIRT", "cvIRTresub")
  results
}

# --- bestModel ---

test_that("bestModel rejects non-cvIRT objects", {
  expect_error(bestModel(list(a = 1)), "not of class")
})

test_that("bestModel holdout selects correct models by IC", {
  res <- make_holdout_result(aic_vals = c(200, 190), bic_vals = c(210, 205))
  bm <- bestModel(res)

  expect_equal(names(bm$info[1]), "2PL")
  expect_equal(names(bm$info[2]), "2PL")
  expect_equal(names(bm$info[3]), "2PL")
})

test_that("bestModel holdout selects simpler model when LRT not significant", {
  res <- make_holdout_result(lrt_pval = 0.50)
  bm <- bestModel(res)

  expect_equal(bm$lrt[[1]]["Best Model"], c("Best Model" = "1PL"))
})

test_that("bestModel holdout selects complex model when LRT is significant", {
  res <- make_holdout_result(lrt_pval = 0.001)
  bm <- bestModel(res)

  expect_equal(bm$lrt[[1]]["Best Model"], c("Best Model" = "2PL"))
})

test_that("bestModel resub works", {
  res <- make_resub_result(aic_vals = c(200, 190))
  bm <- bestModel(res)

  expect_equal(names(bm$info[1]), "2PL")
  expect_s3_class(bm, "cvIRT.bestModels")
})

test_that("bestModel resub selects simpler model when not significant", {
  res <- make_resub_result(lrt_pval = 0.80)
  bm <- bestModel(res)

  expect_equal(bm$lrt[[1]]["Best Model"], c("Best Model" = "1PL"))
})

# --- extract.cvIRT.bestModels ---

test_that("extract returns named character vector of best models", {
  res <- make_holdout_result()
  bm <- bestModel(res)
  ex <- extract.cvIRT.bestModels(bm)

  expect_type(ex, "character")
  expect_named(ex, c("AIC", "BIC", "AICc", "LRT"))
})

# --- print.cvIRT.bestModels ---

test_that("print.cvIRT.bestModels produces output without error", {
  res <- make_holdout_result()
  bm <- bestModel(res)

  output <- capture.output(print(bm))
  expect_true(any(grepl("AIC:", output)))
  expect_true(any(grepl("BIC:", output)))
  expect_true(any(grepl("LRT:", output)))
  # Verify signif digits arg is working (no stray "3" at end of line)
  expect_false(any(grepl("Value: [0-9.]+3$", output)))
})
