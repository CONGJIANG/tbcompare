# Cross-trial g-computation prediction helper

Fits an outcome regression model within one trial and predicts the
counterfactual mean outcomes under two treatment values in the other
trial.

## Usage

``` r
gcomp_cross(
  data,
  s_val,
  treat_active,
  treat_control,
  covars = c("L1", "L2"),
  outcome = "y",
  quiet = TRUE
)
```

## Arguments

- data:

  A pooled trial data frame.

- s_val:

  Trial indicator value used to fit the outcome model.

- treat_active:

  Treatment label for the active regimen.

- treat_control:

  Treatment label for the comparator regimen.

- covars:

  Character vector of covariate names included in the outcome model.

- outcome:

  Name of the binary outcome variable.

- quiet:

  Logical; if `FALSE`, prints the fitted outcome model coefficients.

## Value

A list with predicted mean outcomes under the active and comparator
treatment values in the opposite trial population.

## Details

This is an internal helper used by the g-computation estimators.
