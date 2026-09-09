# Tests for input validation in exported functions.
# These test the guard clauses before any TAM model fitting runs.

test_that("crossValidation rejects non-matrix input", {
  expect_error(crossValidation("not a matrix", modelTypes = "1PL"),
               "matrix or data.frame")
})

test_that("crossValidation warns on non-numeric seed", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      crossValidation(m, modelTypes = "1PL", seed = "bad"),
      error = function(e) NULL
    ),
    "numeric value"
  )
})

test_that("crossValidation warns on invalid folds", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      crossValidation(m, modelTypes = "1PL", folds = -1, seed = 42),
      error = function(e) NULL
    ),
    "folds"
  )
})

test_that("crossValidation warns on invalid replications", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      crossValidation(m, modelTypes = "1PL", replications = -1, seed = 42),
      error = function(e) NULL
    ),
    "replications"
  )
})

test_that("crossValidation forces replications=1 for LOOCV", {
  m <- matrix(sample(0:1, 30, replace = TRUE), nrow = 5)
  expect_warning(
    tryCatch(
      crossValidation(m, modelTypes = "1PL", folds = 5, replications = 3, seed = 42),
      error = function(e) NULL
    ),
    "replications.*forced"
  )
})

test_that("holdout rejects non-matrix input", {
  expect_error(holdout("not a matrix", modelTypes = "1PL", proportion = .3, replications = 1),
               "matrix or data.frame")
})

test_that("holdout warns on invalid proportion", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      holdout(m, modelTypes = "1PL", proportion = 1.5, replications = 1, seed = 42),
      error = function(e) NULL
    ),
    "proportion"
  )
  expect_warning(
    tryCatch(
      holdout(m, modelTypes = "1PL", proportion = 0, replications = 1, seed = 42),
      error = function(e) NULL
    ),
    "proportion"
  )
})

test_that("holdout warns on non-numeric seed", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      holdout(m, modelTypes = "1PL", proportion = .3, replications = 1, seed = "bad"),
      error = function(e) NULL
    ),
    "numeric value"
  )
})

test_that("simpleBootstrap rejects non-matrix input", {
  expect_error(simpleBootstrap("not a matrix", modelTypes = "1PL"),
               "matrix or data.frame")
})

test_that("simpleBootstrap warns on invalid bootSize", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      simpleBootstrap(m, modelTypes = "1PL", bootSize = -5, seed = 42),
      error = function(e) NULL
    ),
    "bootSize"
  )
})

test_that("kfoldBootstrap rejects non-matrix input", {
  expect_error(kfoldBootstrap("not a matrix", modelTypes = "1PL"),
               "matrix or data.frame")
})

test_that("kfoldBootstrap warns on invalid bootSize", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", bootSize = -5, seed = 42),
      error = function(e) NULL
    ),
    "bootSize"
  )
})

test_that("resubstitution rejects non-matrix input", {
  expect_error(resubstitution("not a matrix", modelTypes = "1PL"),
               "matrix or data.frame")
})

test_that("resubstitution warns on non-logical indicator", {
  m <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2)
  expect_warning(
    tryCatch(
      resubstitution(m, modelTypes = "1PL", indicator = "yes"),
      error = function(e) NULL
    ),
    "indicator"
  )
})
