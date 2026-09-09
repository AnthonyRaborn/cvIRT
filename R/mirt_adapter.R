#' Valid model type strings accepted by cvIRT's cross-validation methods.
#' @keywords internal
cvirt_valid_models <- c("Rasch", "1PL", "2PL", "3PL", "PCM", "PCM2", "RSM", "GPCM")

#' Map a cvIRT model type string to a mirt itemtype
#'
#' "1PL" intentionally maps to "2PL" here, not "Rasch" -- the 1PL uses a
#' shared, freely estimated discrimination (see build_constraints()), while
#' Rasch fixes discrimination at 1. Conflating the two collapses them into the
#' same model, which is what this mapping (and fit_irt_model()) is written to
#' avoid.
#'
#' @keywords internal
model_type_to_mirt <- function(model_type, data) {
  switch(model_type,
    "Rasch" = "Rasch",
    "1PL"   = "2PL",
    "2PL"   = "2PL",
    "3PL"   = "3PL",
    "PCM"   = "gpcm",
    "PCM2"  = "gpcm",
    "RSM"   = "rsm",
    "GPCM"  = "gpcm",
    stop("Unknown model_type: ", model_type)
  )
}

#' Build the mirt `constrain` argument for models needing an equality
#' constraint across items.
#'
#' Only the 1PL needs this: its discrimination is estimated (not fixed at 1
#' like Rasch), but shared across all items. The `pars = "values"` call below
#' only returns the parameter table structure -- it does not fit a model.
#'
#' @keywords internal
build_constraints <- function(model_type, data) {
  if (model_type == "1PL") {
    pars <- mirt::mirt(as.data.frame(data), 1, itemtype = "2PL", pars = "values")
    a1_indices <- which(pars$name == "a1" & pars$est)
    if (length(a1_indices) > 1) {
      return(list(a1_indices))
    }
  }
  NULL
}

#' Build a starting parameter table for models that fix specific parameters
#' rather than tying them via `constrain`.
#'
#' - PCM/PCM2: gpcm with discrimination (a1) fixed at 1 for every item.
#' - Rasch: mirt's "Rasch" itemtype fixes item discrimination at 1 but, by
#'   default, freely estimates the latent variance for identification. That
#'   makes it a reparameterization of the 1PL above (same loglik, same
#'   parameter count), erasing the distinction between the two models. Fixing
#'   the latent variance at 1 (the standard Rasch identification convention)
#'   keeps Rasch and 1PL as genuinely different models with different npar.
#'
#' @keywords internal
build_start_pars <- function(model_type, data, itemtype) {
  if (model_type %in% c("PCM", "PCM2")) {
    pars <- mirt::mirt(as.data.frame(data), 1, itemtype = itemtype, pars = "values")
    pars$value[pars$name == "a1"] <- 1
    pars$est[pars$name == "a1"] <- FALSE
    return(pars)
  }
  if (model_type == "Rasch") {
    pars <- mirt::mirt(as.data.frame(data), 1, itemtype = itemtype, pars = "values")
    pars$value[pars$name == "COV_11"] <- 1
    pars$est[pars$name == "COV_11"] <- FALSE
    return(pars)
  }
  NULL
}

#' Extract a fixed-parameter table from a fitted mirt object.
#'
#' Used to re-estimate a model on held-out data (holdout, k-fold, LOOCV)
#' with item parameters frozen at their training-data values. Rather than
#' setting `est <- FALSE` on every row, previously-free parameters are pinned
#' by collapsing their bounds to a single point. Several itemtypes (`rsm`,
#' and the 1PL's shared-a1 constraint) register internal equality
#' constraints keyed off the `est` flag; marking those rows non-free breaks
#' mirt's own constraint validation ("Equality constraints can only be
#' applied to freely estimated parameters"). Bound-pinning keeps the rows
#' formally free so those internal constraints stay valid, while preventing
#' the optimizer from moving the values away from the training fit.
#'
#' The `K` attribute (categories per item) is required by fit_irt_model()
#' when refitting on a single row of data, where mirt cannot otherwise infer
#' how many response categories an item has.
#'
#' @keywords internal
extract_fixed_pars <- function(fit) {
  pars <- mirt::mod2values(fit)
  free <- pars$est
  pars$lbound[free] <- pars$value[free]
  pars$ubound[free] <- pars$value[free]
  attr(pars, "K") <- mirt::extract.mirt(fit, "K")
  pars
}

#' Fit a single IRT model via mirt, returning a standardized result list.
#'
#' This is the single entry point every cvIRT CV method calls instead of
#' `TAM::` directly, isolating the mirt dependency to this file.
#'
#' @param data Response matrix (persons x items).
#' @param model_type One of `cvirt_valid_models`.
#' @param fixed_pars `NULL` for free estimation, or a fixed-parameter table
#'   from `extract_fixed_pars()` to fit on new data with item parameters held
#'   at previously estimated values (holdout/k-fold/LOOCV test fits).
#' @param verbose Passed through to `mirt::mirt(verbose = ...)`.
#' @param ... Additional arguments forwarded to `mirt::mirt()`.
#'
#' @return A list with `fit` (the mirt object), `loglik`, `npar`, and
#'   `fixed_pars` (ready to pass to a subsequent `fit_irt_model()` call).
#'
#' @keywords internal
fit_irt_model <- function(data, model_type, fixed_pars = NULL, verbose = FALSE, ...) {
  itemtype <- model_type_to_mirt(model_type, data)
  model <- 1
  dots <- list(...)

  if (!is.null(fixed_pars)) {
    technical <- dots$technical
    technical$customK <- attr(fixed_pars, "K")
    dots$technical <- technical
    fit <- do.call(mirt::mirt, c(list(
      data = as.data.frame(data), model = model,
      itemtype = itemtype, pars = fixed_pars,
      verbose = verbose
    ), dots))
  } else {
    start_pars <- build_start_pars(model_type, data, itemtype)
    if (!is.null(start_pars)) {
      fit <- do.call(mirt::mirt, c(list(
        data = as.data.frame(data), model = model,
        itemtype = itemtype, pars = start_pars,
        verbose = verbose
      ), dots))
    } else {
      constrain_list <- build_constraints(model_type, data)
      fit <- do.call(mirt::mirt, c(list(
        data = as.data.frame(data), model = model,
        itemtype = itemtype, constrain = constrain_list,
        verbose = verbose
      ), dots))
    }
  }

  loglik <- as.numeric(mirt::extract.mirt(fit, "logLik"))
  npar   <- mirt::extract.mirt(fit, "nest")

  list(
    fit = fit,
    loglik = loglik,
    npar = npar,
    fixed_pars = extract_fixed_pars(fit)
  )
}
