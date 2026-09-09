# Extended tests for model_fit_internals.R covering uncovered branches.

cvAIC <- cvIRT:::cvAIC
cvAICc <- cvIRT:::cvAICc
cvBIC <- cvIRT:::cvBIC
cvLogLikRatio <- cvIRT:::cvLogLikRatio
loob632logLikEst <- cvIRT:::loob632logLikEst

# --- Fixtures ---

make_kfold_bootstrap_ll <- function() {
  list(
    matrix(c(-100, -105, -90, -95), nrow = 2, ncol = 2,
           dimnames = list(NULL, c("1PL", "2PL"))),
    matrix(c(-98, -103, -88, -93), nrow = 2, ncol = 2,
           dimnames = list(NULL, c("1PL", "2PL")))
  )
}

make_kfold_bootstrap_np <- function() {
  list(
    matrix(c(5, 5, 10, 10), nrow = 2, ncol = 2),
    matrix(c(5, 5, 10, 10), nrow = 2, ncol = 2)
  )
}

# --- cvAICc: AICc guard (n <= numParams + 2) ---

test_that("cvAICc returns NA with warning when n <= numParams + 2 (holdout)", {
  ll <- matrix(c(-100, -90), nrow = 1, ncol = 2, dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 10), nrow = 1, ncol = 2)
  expect_warning(
    result <- cvAICc(ll, np, n = 7, method = "holdout"),
    "AICc correction is undefined"
  )
  expect_true(all(is.na(result[1, ])))
})

test_that("cvAICc returns NA with warning when n <= numParams + 2 (resub)", {
  ll <- matrix(c(-100, -90), nrow = 1, ncol = 2, dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 10), nrow = 1, ncol = 2)
  expect_warning(
    result <- cvAICc(ll, np, n = 7, method = "resub"),
    "AICc correction is undefined"
  )
  expect_true(all(is.na(result)))
})

test_that("cvAICc returns NA with warning when n <= numParams + 2 (loob)", {
  ll <- matrix(c(-100, -90), nrow = 1, ncol = 2, dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 10), nrow = 1, ncol = 2)
  expect_warning(
    result <- cvAICc(ll, np, n = 1, method = "loob"),
    "AICc correction is undefined"
  )
  expect_true(all(is.na(result)))
  expect_equal(rownames(result), "Leave One Out Bootstrap:")
})

test_that("cvAICc returns NA with warning when n <= numParams + 2 (kfold bootstrap)", {
  ll <- make_kfold_bootstrap_ll()
  np <- make_kfold_bootstrap_np()
  expect_warning(
    result <- cvAICc(ll, np, n = 7, method = "kfold bootstrap"),
    "AICc correction is undefined"
  )
  expect_true(all(is.na(result)))
})

# --- cvAICc: kfold bootstrap normal path ---

test_that("cvAICc kfold bootstrap computes with correction", {
  ll <- make_kfold_bootstrap_ll()
  np <- make_kfold_bootstrap_np()
  result <- cvAICc(ll, np, n = 100, method = "kfold bootstrap")

  expect_equal(nrow(result), 3)
  expect_equal(rownames(result)[3], "Mean:")
  expect_true(all(grepl("Mean AICc", rownames(result)[1:2])))
})

# --- cvAICc: kfold method ---

test_that("cvAICc kfold method computes correctly", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  result <- cvAICc(ll, np, n = 100, method = "kfold")

  expect_equal(nrow(result), 4)
  expect_true(grepl("^Fold ", rownames(result)[1]))
  expect_equal(rownames(result)[4], "Mean:")
})

test_that("cvAICc bootstrap method computes correctly", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  result <- cvAICc(ll, np, n = 100, method = "bootstrap")

  expect_equal(nrow(result), 4)
  expect_true(grepl("^Bootstrap Sample ", rownames(result)[1]))
})

test_that("cvAICc handles NA in mean calculation", {
  ll <- matrix(c(-100, NA, -90, -95), nrow = 2, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 10, 10), nrow = 2, ncol = 2)
  result <- cvAICc(ll, np, n = 100, method = "holdout")

  expect_false(all(is.na(result[nrow(result), ])))
})

# --- cvBIC: kfold bootstrap ---

test_that("cvBIC kfold bootstrap method works", {
  ll <- make_kfold_bootstrap_ll()
  np <- make_kfold_bootstrap_np()
  result <- cvBIC(ll, np, n = 100, method = "kfold bootstrap")

  expect_equal(nrow(result), 3)
  expect_equal(rownames(result)[3], "Mean:")
})

# --- cvBIC: kfold ---

test_that("cvBIC kfold method works", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  result <- cvBIC(ll, np, n = 100, method = "kfold")

  expect_equal(nrow(result), 4)
  expect_true(grepl("^Fold ", rownames(result)[1]))
})

# --- cvBIC: bootstrap ---

test_that("cvBIC bootstrap method works", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  result <- cvBIC(ll, np, n = 100, method = "bootstrap")

  expect_equal(nrow(result), 4)
  expect_true(grepl("^Bootstrap Sample ", rownames(result)[1]))
})

test_that("cvBIC handles NA in mean calculation", {
  ll <- matrix(c(-100, NA, -90, -95), nrow = 2, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 10, 10), nrow = 2, ncol = 2)
  result <- cvBIC(ll, np, n = 100, method = "holdout")

  expect_false(all(is.na(result[nrow(result), ])))
})

# --- cvAIC: loob ---

test_that("cvAIC loob method uses last row", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  result <- cvAIC(ll, np, method = "loob")

  expected <- -2 * ll[3, ] + 2 * np[1, ]
  expect_equal(as.numeric(result), unname(expected))
  expect_equal(rownames(result), "Leave One Out Bootstrap:")
})

# --- cvBIC: loob ---

test_that("cvBIC loob method works", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  result <- cvBIC(ll, np, n = 100, method = "loob")

  expected <- -2 * ll[3, ] + np[1, ] * log(100)
  expect_equal(as.numeric(result), unname(expected))
  expect_equal(rownames(result), "Leave One Out Bootstrap:")
})

# --- cvBIC: resub rownames ---

test_that("cvBIC resub has correct rowname", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  result <- cvBIC(ll, np, n = 100, method = "resub")

  expect_equal(rownames(result), "Resubstitution")
})

# --- cvLogLikRatio: kfold bootstrap ---

test_that("cvLogLikRatio kfold bootstrap method works", {
  ll <- make_kfold_bootstrap_ll()
  np <- make_kfold_bootstrap_np()
  models <- c("1PL", "2PL")

  result <- cvLogLikRatio(ll, np, models, method = "kfold bootstrap")

  expect_type(result, "list")
  expect_true("1PL vs. 2PL" %in% names(result))
  expect_equal(rownames(result[["1PL vs. 2PL"]])[nrow(result[["1PL vs. 2PL"]])], "Mean:")
  expect_equal(colnames(result[["1PL vs. 2PL"]]), c("Test Statistic", "degrees of freedom", "p-value"))
})

# --- cvLogLikRatio: loob ---

test_that("cvLogLikRatio loob method works", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  models <- c("1PL", "2PL")

  result <- cvLogLikRatio(ll, np, models, method = "loob")

  expect_type(result, "list")
  expect_length(result, 1)
  expect_equal(nrow(result[[1]]), 1)
  expect_equal(rownames(result[[1]]), "Leave One Out Bootstrap:")
})

# --- cvLogLikRatio: kfold ---

test_that("cvLogLikRatio kfold method works", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  models <- c("1PL", "2PL")

  result <- cvLogLikRatio(ll, np, models, method = "kfold")

  expect_true(grepl("^Fold ", rownames(result[[1]])[1]))
  expect_equal(rownames(result[[1]])[nrow(result[[1]])], "Mean:")
})

# --- cvLogLikRatio: bootstrap ---

test_that("cvLogLikRatio bootstrap method works", {
  ll <- matrix(c(-100, -105, -98, -90, -95, -88), nrow = 3, ncol = 2,
               dimnames = list(NULL, c("1PL", "2PL")))
  np <- matrix(c(5, 5, 5, 10, 10, 10), nrow = 3, ncol = 2)
  models <- c("1PL", "2PL")

  result <- cvLogLikRatio(ll, np, models, method = "bootstrap")

  expect_true(grepl("^Bootstrap Sample ", rownames(result[[1]])[1]))
})

# --- loob632logLikEst ---

test_that("loob632logLikEst computes .632 weighted combination", {
  loob <- matrix(c(-100, -95, -90, -85), nrow = 2, ncol = 2,
                 dimnames = list(c("Obs 1:", "Leave-One-Out Mean:"), c("1PL", "2PL")))
  resub <- matrix(c(-80, -75, -70, -65), nrow = 2, ncol = 2,
                  dimnames = list(c("Obs 1:", "Resubstitution Mean:"), c("1PL", "2PL")))

  result <- loob632logLikEst(loob, resub)

  expected <- .632 * loob[2, ] + .368 * resub[2, ]
  expect_equal(as.numeric(result), unname(expected))
  expect_equal(rownames(result), ".632 Bootstrap log-Likelihood:")
})
