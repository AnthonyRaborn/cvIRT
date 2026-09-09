# cvIRT (development version)

## Breaking changes

* The IRT estimation backend has moved from `TAM` to `mirt`
  (`Imports: mirt, stats`). `...` arguments passed to `resubstitution()`,
  `holdout()`, `crossValidation()`, `simpleBootstrap()`, and
  `kfoldBootstrap()` are now forwarded to `mirt::mirt()`, not
  `TAM::tam.mml()`/`TAM::tam.mml.2pl()`. TAM-specific control arguments
  (e.g. `control = list(snodes = 2000)`) have no effect under mirt and
  must be replaced with mirt equivalents (e.g.
  `technical = list(NCYCLES = 2000)`).
* `modelTypes` now accepts `"Rasch"` and `"3PL"` (new) and no longer
  accepts `"2PL.groups"` (removed -- it required TAM's group-specific
  discrimination machinery with no straightforward mirt equivalent
  short of `mirt::multipleGroup()`, which changes the calling
  convention).
* `"Rasch"` and `"1PL"` are now distinct models. Previously `"1PL"` was
  estimated identically to Rasch (all discriminations fixed at 1).
  `"1PL"` now estimates a single shared discrimination across items
  (mirt's `"2PL"` itemtype with an equality constraint on `a1`), while
  `"Rasch"` fixes discrimination at 1 and estimates item difficulties
  only. The two models have different numbers of parameters and will
  generally produce different loglik/AIC/BIC/AICc values on the same
  data.
* Loglikelihood, AIC, AICc, and BIC values will differ slightly from
  previous versions even for models with the same name (e.g. `"2PL"`,
  `"PCM"`), because TAM and mirt use different estimation algorithms
  (TAM: marginal MML with EM; mirt: EM with adaptive quadrature by
  default). Model *selection* should generally agree, but exact numeric
  equality with TAM-based results should not be expected.

## Other changes

* Internal fitted-model objects stored in `trainModel`/`testModel`
  elements are no longer TAM fit objects (class `"tam.mml"`); they are
  now lists with `fit` (a mirt object), `loglik`, `npar`, and
  `fixed_pars`. This only matters if code accesses those internal
  objects directly -- the exported API (`bestModel()`,
  `extract.cvIRT.bestModels()`, and the IC/LRT values in the returned
  `cvIRT` objects) is unaffected.

# cvIRT 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
* First uploaded package to GitHub.
* NOTE: STILL IN ALPHA - functions are more or less set, but methods (show, print, summary) are still being built.
