# TMLE estimator of the direct effect

## Usage

``` r
tmle_direct(
  data,
  nuisance_type = c("simple", "flexible"),
  a_val = "a",
  b_val = "b",
  covars = c("L1", "L2"),
  outcome = "y",
  known_pi = NULL,
  trim = 0.01,
  clamp = c(1e-06, 1 - 1e-06)
)
```

## Arguments

- data:

  A data frame containing the observed data. Must include variables `y`,
  `treatment`, `s`, `L1`, and `L2`.

- nuisance_type:

  Nuisance estimation strategy. Either `"simple"` or `"flexible"`.

- a_val:

  Active treatment in trial 1.

- b_val:

  Active treatment in trial 2.

- covars:

  Baseline covariates used in nuisance estimation.

- outcome:

  Outcome variable name.

- known_pi:

  Optional known treatment assignment probabilities.

- trim:

  Truncation level for estimated propensity scores.

- clamp:

  Bounds used to stabilize predicted probabilities.

## Value

A list containing:

- estimate:

  TMLE estimate of the direct effect.

- se:

  Influence-function-based standard error.

- CI:

  95\\ eifEstimated efficient influence function values.
  nuisance_typeNuisance estimation strategy used.

Estimates the direct cross-trial contrast using targeted maximum
likelihood estimation (TMLE). The estimator combines:

- outcome regressions within each trial,

- treatment assignment mechanisms within each trial,

- a trial membership mechanism,

- a targeting step based on the efficient influence function.

Outcome regressions are first estimated separately within each trial.
These initial estimates are then updated through a logistic fluctuation
step targeting the efficient influence function corresponding to the
direct effect parameter.Standard errors are obtained from the empirical
variance of the estimated efficient influence function.
