# cvIRT

Cross-validation and bootstrap methods for selecting between competing IRT (Item Response Theory) and Rasch models. Given a matrix of item responses, `cvIRT` fits candidate models using the [TAM](https://CRAN.R-project.org/package=TAM) package and compares them with AIC, AICc, BIC, and the log-likelihood ratio test with Holm–Bonferroni correction.

## Supported models

Any combination of the following TAM model types can be compared:

- **1PL** (Rasch)
- **2PL**
- **PCM** (Partial Credit Model)
- **PCM2**
- **RSM** (Rating Scale Model)
- **GPCM** (Generalized Partial Credit Model)
- **2PL.groups**

## Installation

```r
# install.packages("devtools")
devtools::install_github("AnthonyRaborn/cvIRT")
```

## Quick start

Every method takes a response matrix (persons × items) and a character vector of model types to compare. All examples below use a small simulated binary dataset:

```r
library(cvIRT)

set.seed(123)
resp <- matrix(sample(0:1, 500, replace = TRUE), nrow = 50, ncol = 10)
colnames(resp) <- paste0("Item", 1:10)
```

### Resubstitution

Fits each model on the full dataset (no data splitting). Fast but optimistically biased — useful as a baseline.

```r
resub_result <- resubstitution(resp, modelTypes = c("1PL", "2PL"), indicator = FALSE)
resub_result
```

### Holdout validation

Randomly splits data into training and test sets. The `proportion` argument controls the fraction held out.

```r
ho_result <- holdout(
  resp,
  modelTypes  = c("1PL", "2PL"),
  proportion  = 0.3,
  replications = 5,
  indicator   = FALSE,
  seed        = 42
)
ho_result
```

### *k*-fold cross-validation

Splits data into *k* folds; each fold serves as the test set once. Setting `folds = nrow(resp)` performs leave-one-out cross-validation (LOOCV).

```r
cv_result <- crossValidation(
  resp,
  modelTypes  = c("1PL", "2PL"),
  folds       = 5,
  replications = 3,
  indicator   = FALSE,
  seed        = 42
)
cv_result
```

### Simple bootstrap

Draws bootstrap samples from the full dataset and fits each model on every sample.

```r
boot_result <- simpleBootstrap(
  resp,
  modelTypes = c("1PL", "2PL"),
  bootSize   = 20,
  indicator  = FALSE,
  seed       = 42
)
boot_result
```

### *k*-fold bootstrap

Combines *k*-fold splitting with bootstrap resampling within each fold.

```r
kfboot_result <- kfoldBootstrap(
  resp,
  modelTypes = c("1PL", "2PL"),
  bootSize   = 10,
  folds      = 5,
  indicator  = FALSE,
  seed       = 42
)
kfboot_result
```

### Model selection

`bestModel()` extracts the best model under each criterion from any result object. It uses the Holm–Bonferroni correction for the log-likelihood ratio tests.

```r
best <- bestModel(cv_result)
best

# Extract just the model names
extract.cvIRT.bestModels(best)
```

## Included data

The package ships with `TIMSSParameters`, a data frame of polytomous item parameters from the TIMSS 2015 4th-grade math assessment. It contains difficulty, step, and discrimination parameters for 64 items under the RSM, PCM, and GPCM models. This dataset is useful for understanding IRT parameterizations but is not response data — the package functions require a person-by-item response matrix.

```r
data(TIMSSParameters)
str(TIMSSParameters)
```

## License

GPL-2
