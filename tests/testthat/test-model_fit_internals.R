# Tests for internal computation functions in model_fit_internals.R
# These functions are not exported, so access via :::

cvAIC <- cvIRT:::cvAIC
cvAICc <- cvIRT:::cvAICc
cvBIC <- cvIRT:::cvBIC
cvLogLikRatio <- cvIRT:::cvLogLikRatio
warningAttr <- cvIRT:::warningAttr
selectRow <- cvIRT:::selectRow
selectNamedVector <- cvIRT:::selectNamedVector
selectListElement <- cvIRT:::selectListElement

# --- Test fixtures ---

make_ll_matrix <- function(nreps = 3, nmodels = 2) {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = nreps, ncol = nmodels)
  colnames(ll) <- c("1PL", "2PL")
  ll
}

make_np_matrix <- function(nreps = 3, nmodels = 2) {
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = nreps, ncol = nmodels)
  colnames(np) <- c("1PL", "2PL")
  np
}

# --- cvAIC ---

test_that("cvAIC holdout method computes -2LL + 2k with correct structure", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  result <- cvAIC(ll, np, method = "holdout")

  expected_aic <- -2 * ll + 2 * np
  expect_equal(nrow(result), nrow(ll) + 1)
  expect_equal(unname(result[1:3, ]), unname(expected_aic))
  expect_equal(rownames(result)[4], "Mean:")
  expect_equal(unname(result[4, ]), unname(colMeans(expected_aic)))
})

test_that("cvAIC kfold method computes correctly", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  result <- cvAIC(ll, np, method = "kfold")

  expected_aic <- -2 * ll + 2 * np
  expect_equal(unname(result[1:3, ]), unname(expected_aic))
  expect_true(grepl("^Fold ", rownames(result)[1]))
})

test_that("cvAIC bootstrap method computes correctly", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  result <- cvAIC(ll, np, method = "bootstrap")

  expect_true(grepl("^Bootstrap Sample ", rownames(result)[1]))
  expect_equal(rownames(result)[4], "Mean:")
})

test_that("cvAIC resub method uses last row of loglikelihood", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  result <- cvAIC(ll, np, method = "resub")

  expected <- -2 * ll[nrow(ll), ] + 2 * np[1, ]
  expect_equal(as.numeric(result), unname(expected))
  expect_equal(nrow(result), 1)
})

test_that("cvAIC handles NA log-likelihoods with warning", {
  ll <- make_ll_matrix()
  ll[2, 1] <- NA
  np <- make_np_matrix()

  expect_warning(
    result <- cvAIC(ll, np, method = "holdout"),
    "missing log-likelihood"
  )
  expect_false(is.null(attr(result, "warning")))
  expect_false(any(is.na(result[nrow(result), ])))
})

# --- cvAICc ---

test_that("cvAICc holdout adds correction term", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  n <- 100
  result <- cvAICc(ll, np, n = n, method = "holdout")

  correction <- (2 * (np + 2) * (np + 1)) / (n - np - 2)
  expected <- -2 * ll + 2 * np + correction
  expect_equal(unname(result[1:3, ]), unname(expected))
})

test_that("cvAICc resub uses last row", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  n <- 100
  result <- cvAICc(ll, np, n = n, method = "resub")

  k <- np[1, ]
  expected <- -2 * ll[nrow(ll), ] + 2 * k + (2 * (k + 2) * (k + 1)) / (n - k - 2)
  expect_equal(as.numeric(result), unname(expected))
})

# --- cvBIC ---

test_that("cvBIC holdout computes -2LL + k*log(n)", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  n <- 100
  result <- cvBIC(ll, np, n = n, method = "holdout")

  expected <- -2 * ll + np * log(n)
  expect_equal(unname(result[1:3, ]), unname(expected))
  expect_equal(unname(result[4, ]), unname(colMeans(expected)))
})

test_that("cvBIC resub uses last row", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  n <- 100
  result <- cvBIC(ll, np, n = n, method = "resub")

  expected <- -2 * ll[nrow(ll), ] + np[1, ] * log(n)
  expect_equal(as.numeric(result), unname(expected))
})

# --- cvLogLikRatio ---

test_that("cvLogLikRatio holdout produces correct test statistics", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  models <- c("1PL", "2PL")

  result <- cvLogLikRatio(ll, np, models, method = "holdout")

  expect_type(result, "list")
  expect_length(result, 1)
  expect_equal(names(result), "1PL vs. 2PL")
  expect_equal(colnames(result[[1]]), c("Test Statistic", "degrees of freedom", "p-value"))
  # Mean row at the end
  expect_equal(rownames(result[[1]])[nrow(result[[1]])], "Mean:")
  # All p-values between 0 and 1
  expect_true(all(result[[1]][, "p-value"] >= 0 & result[[1]][, "p-value"] <= 1))
})

test_that("cvLogLikRatio resub produces single-row result", {
  ll <- make_ll_matrix()
  np <- make_np_matrix()
  models <- c("1PL", "2PL")

  result <- cvLogLikRatio(ll, np, models, method = "resub")

  expect_type(result, "list")
  expect_length(result, 1)
  expect_equal(nrow(result[[1]]), 1)
})

test_that("cvLogLikRatio handles three models", {
  ll <- matrix(c(-100, -95, -90, -80, -75, -70, -70, -65, -60),
               nrow = 3, ncol = 3)
  colnames(ll) <- c("1PL", "PCM", "GPCM")
  np <- matrix(c(5, 5, 5, 10, 10, 10, 15, 15, 15),
               nrow = 3, ncol = 3)
  colnames(np) <- c("1PL", "PCM", "GPCM")

  result <- cvLogLikRatio(ll, np, c("1PL", "PCM", "GPCM"), method = "holdout")

  expect_length(result, 2)
  expect_equal(names(result), c("1PL vs. PCM", "PCM vs. GPCM"))
})

# --- warningAttr ---

test_that("warningAttr extracts warning attributes from indexed elements", {
  x <- list(
    structure(matrix(1), warning = "warn1"),
    matrix(2),
    structure(matrix(3), warning = "warn3")
  )
  attr(x[[1]], "warning") <- "warn1"
  attr(x[[3]], "warning") <- "warn3"

  result <- warningAttr(x, c(1, 3))
  expect_equal(result, c("warn1", "warn3"))
})

test_that("warningAttr returns empty on empty indicator", {
  x <- list(matrix(1), matrix(2))
  result <- warningAttr(x, c())
  expect_null(result)
})

# --- selectRow ---

test_that("selectRow extracts named row from matrix", {
  m <- matrix(1:6, nrow = 3, ncol = 2)
  rownames(m) <- c("a", "b", "Mean:")
  colnames(m) <- c("X", "Y")

  result <- selectRow(m, "Mean:")
  expect_equal(as.numeric(result), c(3, 6))
})

test_that("selectRow extracts named row from list of matrices", {
  m1 <- matrix(1:4, nrow = 2, ncol = 2)
  rownames(m1) <- c("a", "Mean:")
  m2 <- matrix(5:8, nrow = 2, ncol = 2)
  rownames(m2) <- c("a", "Mean:")

  result <- selectRow(list(comp1 = m1, comp2 = m2), "Mean:")
  expect_type(result, "list")
  expect_length(result, 2)
})

# --- selectNamedVector ---

test_that("selectNamedVector extracts named element from vector", {
  v <- c(a = 1, b = 2, c = 3)
  result <- selectNamedVector(v, "b")
  expect_equal(result, c(b = 2))
})

test_that("selectNamedVector works on list of vectors", {
  v1 <- c(a = 1, b = 2)
  v2 <- c(a = 3, b = 4)
  result <- selectNamedVector(list(x = v1, y = v2), "b")
  expect_type(result, "list")
  expect_equal(result[[1]], c(b = 2))
  expect_equal(result[[2]], c(b = 4))
})

# --- selectListElement ---

test_that("selectListElement transposes list-of-named-vectors into named-list-of-vectors", {
  x <- list(
    list("Test Statistic" = 10, "p-value" = 0.01),
    list("Test Statistic" = 12, "p-value" = 0.005)
  )
  result <- selectListElement(x)

  expect_equal(names(result), c("Test Statistic", "p-value"))
  expect_equal(nrow(result$`Test Statistic`), 2)
})

# --- kfold bootstrap methods ---

test_that("cvAIC kfold bootstrap method works with list inputs", {
  ll <- list(
    matrix(c(-100, -105, -90, -95), nrow = 2, ncol = 2,
           dimnames = list(NULL, c("1PL", "2PL"))),
    matrix(c(-98, -103, -88, -93), nrow = 2, ncol = 2,
           dimnames = list(NULL, c("1PL", "2PL")))
  )
  np <- list(
    matrix(c(5, 5, 10, 10), nrow = 2, ncol = 2),
    matrix(c(5, 5, 10, 10), nrow = 2, ncol = 2)
  )

  result <- cvAIC(ll, np, method = "kfold bootstrap")

  expect_equal(nrow(result), 3)
  expect_equal(rownames(result)[3], "Mean:")
  expect_true(all(grepl("Mean AIC", rownames(result)[1:2])))
})
