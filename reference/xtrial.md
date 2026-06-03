# Main interface for cross-trial comparison estimators

Main interface for cross-trial comparison estimators

## Usage

``` r
xtrial(
  data,
  approach = c("traditional", "observational", "direct", "indirect"),
  estimator = c("ipw", "gcomp", "onestep", "tmle"),
  nboot = 100,
  nuisance_type = c("simple", "flexible")
)
```

## Arguments

- data:

  A data frame containing the trial indicator, treatment, outcome, and
  baseline covariates.

- approach:

  Identification strategy. One of "traditional", "observational",
  "direct", or "indirect".

- estimator:

  Estimation method. One of "ipw", "gcomp", "onestep", or "tmle".

- nboot:

  Number of bootstrap samples for bootstrap-based estimators.

- nuisance_type:

  Nuisance estimation type. Currently mainly used for one-step/EIF and
  TMLE estimators; one of "simple" or "flexible".

## Value

An `xtrial` object.
