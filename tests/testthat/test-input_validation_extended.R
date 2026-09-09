# Extended input validation tests covering new guard clauses.

# --- modelTypes validation (all 5 entry points) ---

test_that("crossValidation rejects invalid modelTypes", {
  m <- matrix(sample(0:1, 30, replace = TRUE), nrow = 5)
  expect_error(crossValidation(m, modelTypes = c("1PL", "FAKE"), seed = 42),
               "Unrecognized modelTypes.*FAKE")
})

test_that("holdout rejects invalid modelTypes", {
  m <- matrix(sample(0:1, 30, replace = TRUE), nrow = 5)
  expect_error(holdout(m, modelTypes = "INVALID", proportion = .3, replications = 1, seed = 42),
               "Unrecognized modelTypes.*INVALID")
})

test_that("resubstitution rejects invalid modelTypes", {
  m <- matrix(sample(0:1, 30, replace = TRUE), nrow = 5)
  expect_error(resubstitution(m, modelTypes = "BAD"),
               "Unrecognized modelTypes.*BAD")
})

test_that("simpleBootstrap rejects invalid modelTypes", {
  m <- matrix(sample(0:1, 30, replace = TRUE), nrow = 5)
  expect_error(simpleBootstrap(m, modelTypes = "NOPE", seed = 42),
               "Unrecognized modelTypes.*NOPE")
})

test_that("kfoldBootstrap rejects invalid modelTypes", {
  m <- matrix(sample(0:1, 30, replace = TRUE), nrow = 5)
  expect_error(kfoldBootstrap(m, modelTypes = "ZZZ", seed = 42),
               "Unrecognized modelTypes.*ZZZ")
})

# --- folds integer coercion ---

test_that("crossValidation warns and coerces non-integer folds", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      crossValidation(m, modelTypes = "1PL", folds = 2.7, seed = 42),
      error = function(e) NULL
    ),
    "not an integer.*Coercing"
  )
})

test_that("kfoldBootstrap warns and coerces non-integer folds", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", folds = 3.5, seed = 42),
      error = function(e) NULL
    ),
    "not an integer.*Coercing"
  )
})

# --- replications = 0 now rejected ---

test_that("crossValidation warns on replications = 0", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      crossValidation(m, modelTypes = "1PL", replications = 0, seed = 42),
      error = function(e) NULL
    ),
    "replications"
  )
})

test_that("holdout warns on replications = 0", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      holdout(m, modelTypes = "1PL", proportion = .3, replications = 0, seed = 42),
      error = function(e) NULL
    ),
    "replications"
  )
})

test_that("simpleBootstrap warns on replications = 0", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      simpleBootstrap(m, modelTypes = "1PL", replications = 0, seed = 42),
      error = function(e) NULL
    ),
    "replications"
  )
})

test_that("kfoldBootstrap warns on replications = 0", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", replications = 0, seed = 42),
      error = function(e) NULL
    ),
    "replications"
  )
})

# --- indicator validation in functions not yet tested ---

test_that("crossValidation warns on non-logical indicator", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      crossValidation(m, modelTypes = "1PL", indicator = "yes", seed = 42),
      error = function(e) NULL
    ),
    "indicator"
  )
})

test_that("holdout warns on non-logical indicator", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      holdout(m, modelTypes = "1PL", proportion = .3, replications = 1,
              indicator = 42, seed = 42),
      error = function(e) NULL
    ),
    "indicator"
  )
})

test_that("simpleBootstrap warns on non-logical indicator", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      simpleBootstrap(m, modelTypes = "1PL", indicator = "no", seed = 42),
      error = function(e) NULL
    ),
    "indicator"
  )
})

test_that("kfoldBootstrap warns on non-logical indicator", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", indicator = 0, seed = 42),
      error = function(e) NULL
    ),
    "indicator"
  )
})

# --- seed validation ---

test_that("simpleBootstrap warns on non-numeric seed", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      simpleBootstrap(m, modelTypes = "1PL", seed = "bad"),
      error = function(e) NULL
    ),
    "numeric value"
  )
})

test_that("kfoldBootstrap warns on non-numeric seed", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", seed = "bad"),
      error = function(e) NULL
    ),
    "numeric value"
  )
})

# --- kfoldBootstrap folds validation ---

test_that("kfoldBootstrap warns on invalid folds", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", folds = -1, seed = 42),
      error = function(e) NULL
    ),
    "folds"
  )
})

# --- holdout replications validation ---

test_that("holdout warns on non-numeric replications", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      holdout(m, modelTypes = "1PL", proportion = .3, replications = "a", seed = 42),
      error = function(e) NULL
    ),
    "replications"
  )
})

# --- kfoldBootstrap bootSize validation ---

test_that("kfoldBootstrap warns on invalid bootSize", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", bootSize = 0, seed = 42),
      error = function(e) NULL
    ),
    "bootSize"
  )
})

# --- kfoldBootstrap replications validation ---

test_that("kfoldBootstrap warns on non-numeric replications", {
  m <- matrix(sample(0:1, 50, replace = TRUE), nrow = 10)
  expect_warning(
    tryCatch(
      kfoldBootstrap(m, modelTypes = "1PL", replications = "a", seed = 42),
      error = function(e) NULL
    ),
    "replications"
  )
})
