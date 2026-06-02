# EIF-based direct estimator

## Usage

``` r
eif_direct_est(
  dataset,
  a_val = "a",
  b_val = "b",
  clamp = c(1e-06, 1 - 1e-06),
  nuisance_type = c("simple", "flexible")
)
```

## Arguments

- dataset:

  A data frame containing the outcome, treatment, trial indicator, and
  baseline covariates. Must include variables `y`, `treatment`, `s`,
  `L1`, and `L2`.

- a_val:

  Character value indicating the active treatment in trial 1. Default is
  `"a"`.

- b_val:

  Character value indicating the active treatment in trial 2. Default is
  `"b"`.

- clamp:

  Numeric vector of length 2 giving the lower and upper bounds used to
  truncate estimated probabilities. Default is `c(1e-6, 1 - 1e-6)`.

- nuisance_type:

  Character string specifying the nuisance estimation strategy. Options
  are `"simple"` and `"flexible"`.

## Value

A list with the following components:

- eif:

  The estimated EIF contribution for each observation.

- estimate:

  The point estimate of the direct contrast.

- se_if:

  Influence-function-based standard error.

- CI_if:

  Wald-type 95\\ componentsA data frame containing the four EIF
  components. nuisance_typeThe nuisance estimation strategy used.

Estimates the direct cross-trial contrast using an efficient influence
function (EIF)-based estimator. The estimator combines outcome
regressions, within-trial treatment mechanisms, and the trial membership
mechanism. This function estimates the direct cross-trial contrast
comparing treatment `a_val` in trial 1 with treatment `b_val` in trial
2. It fits the following nuisance functions:

- treatment mechanism in trial 1,

- treatment mechanism in trial 2,

- outcome regression in trial 1,

- outcome regression in trial 2,

- trial membership mechanism.

The standard error is computed as the empirical standard deviation of
the estimated EIF divided by the square root of the sample size.
