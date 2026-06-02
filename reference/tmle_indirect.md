# TMLE estimator of the indirect effect

## Usage

``` r
tmle_indirect(
  data,
  nuisance_type = c("simple", "flexible"),
  a_val = "a",
  b_val = "b",
  c1_val = "sc1",
  c2_val = "sc2",
  covars = c("L1", "L2"),
  outcome = "y",
  known_pi = NULL,
  trim = 0.01,
  clamp = c(1e-06, 1 - 1e-06)
)
```

## Arguments

- data:

  A data frame containing the observed data.

- nuisance_type:

  Nuisance estimation strategy.

- a_val:

  Active treatment in trial 1.

- b_val:

  Active treatment in trial 2.

- c1_val:

  Standard-of-care treatment in trial 1.

- c2_val:

  Standard-of-care treatment in trial 2.

- covars:

  Baseline covariates used in nuisance estimation.

- outcome:

  Outcome variable name.

- known_pi:

  Optional known treatment probabilities.

- trim:

  Truncation level for estimated propensity scores.

- clamp:

  Bounds used to stabilize predicted probabilities.

## Value

A list containing:

- estimate:

  TMLE estimate of the indirect effect.

- se:

  Influence-function-based standard error.

- CI:

  95\\ eifEstimated efficient influence function values for the indirect
  effect. componentsEstimates and inference for the \\\theta\\ and
  \\\phi\\ components. nuisance_typeNuisance estimation strategy used.

Estimates the indirect effect using the decomposition \$\$ \delta =
\theta + \phi, \$\$where:

- \\\theta\\ is estimated by
  [`tmle_theta()`](https://CONGJIANG.github.io/tbcompare/reference/tmle_theta.md),

- \\\phi\\ is estimated by
  [`tmle_direct()`](https://CONGJIANG.github.io/tbcompare/reference/tmle_direct.md)
  evaluated at the two standard-of-care treatments.

The indirect effect estimator is constructed as\$\$ \hat\delta =
\hat\theta\_{\mathrm{TMLE}} + \hat\phi\_{\mathrm{TMLE}}. \$\$Because
both components are estimated from the same sample, inference is based
on the combined efficient influence function\$\$ D\_{\delta} =
D\_{\theta} + D\_{\phi}, \$\$which automatically accounts for the
covariance between \\\hat\theta\\ and \\\hat\phi\\.
[`tmle_theta()`](https://CONGJIANG.github.io/tbcompare/reference/tmle_theta.md),
[`tmle_direct()`](https://CONGJIANG.github.io/tbcompare/reference/tmle_direct.md)
