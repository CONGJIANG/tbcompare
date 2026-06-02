# Unified IPW estimators for cross-trial comparisons

Computes inverse probability weighted estimators for four cross-trial
comparison strategies:

## Usage

``` r
estimate_ipw_unified(
  dataset,
  a_val = "a",
  b_val = "b",
  c1_val = "sc1",
  c2_val = "sc2",
  clamp = c(1e-06, 1 - 1e-06),
  compute_boot = FALSE,
  B = 200
)
```

## Arguments

- dataset:

  A data frame containing the pooled trial data.

- a_val:

  Treatment label for regimen A.

- b_val:

  Treatment label for regimen B.

- c1_val:

  Standard-of-care comparator in trial 1.

- c2_val:

  Standard-of-care comparator in trial 2.

- clamp:

  Bounds used to truncate estimated probabilities.

- compute_boot:

  Logical; whether bootstrap inference should be computed.

- B:

  Number of bootstrap replicates.

## Value

A list containing estimates, standard errors, confidence intervals, and
influence-function quantities for all implemented comparison strategies.

## Details

- Traditional indirect comparison

- Observational comparison

- Direct pooled IPD comparison

- Indirect decomposition comparison

## See also

[`tb_compare()`](https://CONGJIANG.github.io/tbcompare/reference/tb_compare.md)
