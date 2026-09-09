# Tests for the mirt adapter layer in R/mirt_adapter.R
# These functions are not exported, so access via :::

fit_irt_model <- cvIRT:::fit_irt_model
extract_fixed_pars <- cvIRT:::extract_fixed_pars
model_type_to_mirt <- cvIRT:::model_type_to_mirt

# --- Test fixtures ---

set.seed(123)
dich_data <- matrix(sample(0:1, 300, replace = TRUE), ncol = 10)
poly_data <- matrix(sample(0:2, 300, replace = TRUE), ncol = 10)

# --- fit_irt_model() across all model types ---

test_that("fit_irt_model works for each dichotomous model type", {
  for (mt in c("Rasch", "1PL", "2PL", "3PL")) {
    result <- fit_irt_model(dich_data, mt, verbose = FALSE)
    expect_true(is.finite(result$loglik), info = mt)
    expect_true(result$npar > 0, info = mt)
    expect_s4_class(result$fit, "SingleGroupClass")
  }
})

test_that("fit_irt_model works for each polytomous model type", {
  for (mt in c("PCM", "RSM", "GPCM")) {
    result <- fit_irt_model(poly_data, mt, verbose = FALSE)
    expect_true(is.finite(result$loglik), info = mt)
    expect_true(result$npar > 0, info = mt)
  }
})

# --- fixed_pars re-estimation ---

test_that("fit_irt_model with fixed_pars returns a finite loglik", {
  train <- fit_irt_model(dich_data, "2PL", verbose = FALSE)
  test <- fit_irt_model(dich_data, "2PL", fixed_pars = train$fixed_pars, verbose = FALSE)
  expect_true(is.finite(test$loglik))
})

test_that("fit_irt_model works on single-row data with fixed_pars (the LOOCV case)", {
  for (mt in c("Rasch", "1PL", "2PL", "3PL")) {
    train <- fit_irt_model(dich_data, mt, verbose = FALSE)
    single <- fit_irt_model(matrix(dich_data[1, ], nrow = 1), mt,
                            fixed_pars = train$fixed_pars, verbose = FALSE)
    expect_true(is.finite(single$loglik), info = mt)
  }
  for (mt in c("PCM", "RSM", "GPCM")) {
    train <- fit_irt_model(poly_data, mt, verbose = FALSE)
    single <- fit_irt_model(matrix(poly_data[1, ], nrow = 1), mt,
                            fixed_pars = train$fixed_pars, verbose = FALSE)
    expect_true(is.finite(single$loglik), info = mt)
  }
})

# --- extract_fixed_pars() ---

test_that("extract_fixed_pars pins previously-free parameters to a single value", {
  fit <- fit_irt_model(dich_data, "2PL", verbose = FALSE)
  pars <- fit$fixed_pars

  expect_s3_class(pars, "data.frame")
  free <- pars$est
  expect_true(any(free))
  expect_equal(pars$lbound[free], pars$value[free])
  expect_equal(pars$ubound[free], pars$value[free])
})

test_that("extract_fixed_pars carries a K attribute matching category counts", {
  fit <- fit_irt_model(poly_data, "GPCM", verbose = FALSE)
  K <- attr(fit$fixed_pars, "K")
  expect_length(K, ncol(poly_data))
  expect_true(all(K == 3))
})

# --- model_type_to_mirt() ---

test_that("model_type_to_mirt maps all valid model names correctly", {
  expect_equal(model_type_to_mirt("Rasch", dich_data), "Rasch")
  expect_equal(model_type_to_mirt("1PL", dich_data), "2PL")
  expect_equal(model_type_to_mirt("2PL", dich_data), "2PL")
  expect_equal(model_type_to_mirt("3PL", dich_data), "3PL")
  expect_equal(model_type_to_mirt("PCM", poly_data), "gpcm")
  expect_equal(model_type_to_mirt("PCM2", poly_data), "gpcm")
  expect_equal(model_type_to_mirt("RSM", poly_data), "rsm")
  expect_equal(model_type_to_mirt("GPCM", poly_data), "gpcm")
})

test_that("model_type_to_mirt does not map 1PL to Rasch", {
  expect_false(identical(model_type_to_mirt("1PL", dich_data),
                         model_type_to_mirt("Rasch", dich_data)))
})

test_that("model_type_to_mirt errors on an unknown model type", {
  expect_error(model_type_to_mirt("BOGUS", dich_data), "Unknown model_type")
})

# --- Rasch vs 1PL distinctness ---

test_that("Rasch and 1PL produce different npar on the same data", {
  rasch <- fit_irt_model(dich_data, "Rasch", verbose = FALSE)
  onepl <- fit_irt_model(dich_data, "1PL", verbose = FALSE)

  expect_false(rasch$npar == onepl$npar)
  expect_true(rasch$npar < onepl$npar)
})
