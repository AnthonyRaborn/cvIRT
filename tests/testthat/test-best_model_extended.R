# Tests for bestModel handlers: kfold, bootstrap, and loob class paths.

bestModel <- cvIRT:::bestModel
selectRow <- cvIRT:::selectRow
selectListElement <- cvIRT:::selectListElement

# --- kfold / bootstrap handler ---

make_kfold_result <- function(lrt_pval = 0.01) {
  models <- c("1PL", "2PL")

  aic1 <- matrix(c(200, 190), nrow = 1, dimnames = list("Fold 1:", models))
  aic1 <- rbind(aic1, "Mean:" = c(200, 190))
  aic2 <- matrix(c(205, 192), nrow = 1, dimnames = list("Fold 2:", models))
  aic2 <- rbind(aic2, "Mean:" = c(205, 192))

  bic1 <- matrix(c(210, 205), nrow = 1, dimnames = list("Fold 1:", models))
  bic1 <- rbind(bic1, "Mean:" = c(210, 205))
  bic2 <- matrix(c(215, 208), nrow = 1, dimnames = list("Fold 2:", models))
  bic2 <- rbind(bic2, "Mean:" = c(215, 208))

  aicc1 <- matrix(c(201, 191), nrow = 1, dimnames = list("Fold 1:", models))
  aicc1 <- rbind(aicc1, "Mean:" = c(201, 191))
  aicc2 <- matrix(c(206, 193), nrow = 1, dimnames = list("Fold 2:", models))
  aicc2 <- rbind(aicc2, "Mean:" = c(206, 193))

  lrt_mat <- matrix(c(15.0, 5.0, lrt_pval), nrow = 1, ncol = 3)
  colnames(lrt_mat) <- c("Test Statistic", "degrees of freedom", "p-value")
  rownames(lrt_mat) <- "Mean:"
  lrt1 <- list("1PL vs. 2PL" = lrt_mat)
  lrt2 <- list("1PL vs. 2PL" = lrt_mat)

  results <- list(
    AIC = list(aic1, aic2),
    BIC = list(bic1, bic2),
    AICc = list(aicc1, aicc2),
    `-2 log-Likelihood Ratio Test` = list(lrt1, lrt2)
  )
  class(results) <- c("cvIRT", "cvIRTkfold")
  results
}

test_that("bestModel kfold selects correct models by IC", {
  res <- make_kfold_result()
  bm <- bestModel(res)

  expect_s3_class(bm, "cvIRT.bestModels")
  expect_equal(names(bm$info[1]), "2PL")
  expect_equal(names(bm$info[2]), "2PL")
  expect_equal(names(bm$info[3]), "2PL")
})

test_that("bestModel kfold selects simpler model when LRT not significant", {
  res <- make_kfold_result(lrt_pval = 0.80)
  bm <- bestModel(res)

  expect_equal(bm$lrt[[1]]["Best Model"], c("Best Model" = "1PL"))
})

test_that("bestModel kfold selects complex model when LRT is significant", {
  res <- make_kfold_result(lrt_pval = 0.001)
  bm <- bestModel(res)

  expect_equal(bm$lrt[[1]]["Best Model"], c("Best Model" = "2PL"))
})

# bootstrap class uses the same handler as kfold
test_that("bestModel bootstrap handler works", {
  res <- make_kfold_result()
  class(res) <- c("cvIRT", "cvIRTbootstrap")
  bm <- bestModel(res)

  expect_s3_class(bm, "cvIRT.bestModels")
  expect_equal(names(bm$info[1]), "2PL")
})

# --- loob handler ---

make_loob_result <- function(lrt_pval = 0.01) {
  models <- c("1PL", "2PL")

  aic1 <- matrix(c(200, 190), nrow = 1, dimnames = list("Leave One Out Bootstrap:", models))
  aic2 <- matrix(c(205, 192), nrow = 1, dimnames = list("Leave One Out Bootstrap:", models))

  bic1 <- matrix(c(210, 205), nrow = 1, dimnames = list("Leave One Out Bootstrap:", models))
  bic2 <- matrix(c(215, 208), nrow = 1, dimnames = list("Leave One Out Bootstrap:", models))

  aicc1 <- matrix(c(201, 191), nrow = 1, dimnames = list("Leave One Out Bootstrap:", models))
  aicc2 <- matrix(c(206, 193), nrow = 1, dimnames = list("Leave One Out Bootstrap:", models))

  lrt_mat <- matrix(c(15.0, 5.0, lrt_pval), nrow = 1, ncol = 3)
  colnames(lrt_mat) <- c("Test Statistic", "degrees of freedom", "p-value")
  rownames(lrt_mat) <- "Leave One Out Bootstrap:"
  lrt1 <- list("1PL vs. 2PL" = lrt_mat)
  lrt2 <- list("1PL vs. 2PL" = lrt_mat)

  results <- list(
    AIC = list(aic1, aic2),
    BIC = list(bic1, bic2),
    AICc = list(aicc1, aicc2),
    `-2 log-Likelihood Ratio Test` = list(lrt1, lrt2)
  )
  class(results) <- c("cvIRT", "cvIRTloob")
  results
}

test_that("bestModel loob selects correct models by IC", {
  res <- make_loob_result()
  bm <- bestModel(res)

  expect_s3_class(bm, "cvIRT.bestModels")
  expect_equal(names(bm$info[1]), "2PL")
  expect_equal(names(bm$info[2]), "2PL")
  expect_equal(names(bm$info[3]), "2PL")
})

test_that("bestModel loob selects simpler model when LRT not significant", {
  res <- make_loob_result(lrt_pval = 0.80)
  bm <- bestModel(res)

  expect_equal(bm$lrt[[1]]["Best Model"], c("Best Model" = "1PL"))
})
