# EIF-based estimator of the indirect pathway contrast (\\\theta\\)

Estimates the intermediate contrast \$\$ \theta = \\E\[Y^a -
Y^{c_1}\]\\\_{S=1 \rightarrow 2} - \\E\[Y^b - Y^{c_2}\]\\\_{S=2
\rightarrow 1}, \$\$ using an efficient influence function (EIF)-based
estimator.

## Usage

``` r
eif_theta_est(
  dataset,
  a = "a",
  c1 = "sc1",
  b = "b",
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

- c1:

  Character value identifying the standard-of-care treatment in trial 1.

- b:

  Character value identifying the active treatment in trial 2.

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

- theta:

  Estimated value of \\\theta\\.

- eif:

  Estimated efficient influence function values.

- se:

  Influence-function-based standard error.

- nuisance_type:

  Nuisance estimation strategy used.

## Details

The estimator combines:

- treatment propensity models within each trial,

- outcome regression models within each trial,

- a trial membership model \\P(S=1\mid L)\\.

Outcome regressions are fit separately within each trial and evaluated
under each treatment level of interest. Treatment propensity models are
estimated separately within each trial, while the trial membership model
is estimated using the pooled data.

Standard errors are computed from the empirical variance of the
estimated EIF divided by the sample size.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- eif_theta_est(
  dataset = dat,
  nuisance_type = "simple"
)

fit$theta
fit$se
} # }
```
