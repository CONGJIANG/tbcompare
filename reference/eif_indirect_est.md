# EIF-based estimator of the indirect effect

## Usage

``` r
eif_indirect_est(
  dataset,
  a = "a",
  b = "b",
  c1 = "sc1",
  c2 = "sc2",
  clamp = c(1e-06, 1 - 1e-06),
  nuisance_type = c("simple", "flexible")
)
```

## Arguments

- dataset:

  A data frame containing the observed data. Must include variables `y`,
  `treatment`, `s`, `L1`, and `L2`.

- a:

  Character value identifying the active treatment in trial 1.

- b:

  Character value identifying the active treatment in trial 2.

- c1:

  Character value identifying the standard-of-care treatment in trial 1.

- c2:

  Character value identifying the standard-of-care treatment in trial 2.

- clamp:

  Numeric vector of length two specifying lower and upper truncation
  bounds for estimated probabilities.

- nuisance_type:

  Character string indicating nuisance estimation strategy. Must be one
  of `"simple"` or `"flexible"`.

## Value

A list containing:

- indirect:

  Estimated indirect effect.

- se_if:

  Influence-function-based standard error.

- CI_if:

  Wald-type 95\\ eifEstimated efficient influence function values.
  thetaEstimated pathway-specific component \\\theta\\. phiEstimated
  direct-effect component \\\phi\\. nuisance_typeNuisance estimation
  strategy used.

Estimates the indirect effect by combining the EIF-based estimator of
the pathway-specific contrast \\\theta\\ with the EIF-based estimator of
the direct contrast \\\phi\\. Specifically, \$\$
\psi\_{\mathrm{indirect}} = \theta + \phi. \$\$The indirect effect is
estimated as the sum of:

- the EIF-based estimator of \\\theta\\, obtained from
  [`eif_theta_est()`](https://CONGJIANG.github.io/xtrial/reference/eif_theta_est.md),
  and

- the EIF-based estimator of the direct contrast \\\phi\\, obtained from
  [`eif_direct_est()`](https://CONGJIANG.github.io/xtrial/reference/eif_direct_est.md).

The efficient influence function for the indirect effect is constructed
as the sum of the estimated influence functions for \\\theta\\ and
\\\phi\\. Standard errors and confidence intervals are then obtained
using the empirical variance of this combined EIF.
[`eif_theta_est()`](https://CONGJIANG.github.io/xtrial/reference/eif_theta_est.md),
[`eif_direct_est()`](https://CONGJIANG.github.io/xtrial/reference/eif_direct_est.md)
