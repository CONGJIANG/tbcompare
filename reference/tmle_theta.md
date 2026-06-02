# TMLE estimator of the pathway-specific contrast (\\\theta\\)

## Usage

``` r
tmle_theta(
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

  TMLE estimate of \\\theta\\.

- se:

  Influence-function-based standard error.

- CI:

  95\\ eifEstimated efficient influence function values.
  nuisance_typeNuisance estimation strategy used.

Estimates the pathway-specific contrast \$\$ \theta = \\E\[Y^a -
Y^{c_1}\]\\\_{S=1 \rightarrow 2} - \\E\[Y^b - Y^{c_2}\]\\\_{S=2
\rightarrow 1}. \$\$using targeted maximum likelihood estimation
(TMLE).Four treatment-specific outcome regressions are estimated and
targeted: \\\mu\_{1,a}\\, \\\mu\_{1,c_1}\\, \\\mu\_{2,b}\\, and
\\\mu\_{2,c_2}\\.The resulting TMLE estimator solves the empirical
efficient influence function estimating equation for the parameter
\\\theta\\.
