# Generate cross-trial simulation data

This is the core data-generating mechanism for the cross-trial
comparison setting. The function generates two trials, baseline
covariates, country, treatment assignment, and a binary outcome.

## Usage

``` r
dgp_xtrial_data(
  n = 1000,
  beta0 = -0.5,
  betaL1 = 0.4,
  betaL2 = -0.3,
  betaA_a = 0.8,
  betaA_b = 0.5,
  betaA_sc1 = 0.2,
  delta_c = 0.15,
  betaC2 = 0.3,
  betaC3 = 0.5,
  betaC4 = 0.2
)
```

## Arguments

- n:

  Integer. Sample size.

- beta0:

  Numeric. Outcome-model intercept.

- betaL1:

  Numeric. Coefficient for baseline covariate `L1`.

- betaL2:

  Numeric. Coefficient for baseline covariate `L2`.

- betaA_a:

  Numeric. Treatment effect for treatment `a`.

- betaA_b:

  Numeric. Treatment effect for treatment `b`.

- betaA_sc1:

  Numeric. Treatment effect for standard care version `sc1`.

- delta_c:

  Numeric. Difference between `sc1` and `sc2`; internally,
  `betaA_sc2 = betaA_sc1 - delta_c`.

- betaC2:

  Numeric. Reserved country parameter for compatibility with earlier
  simulation scripts.

- betaC3:

  Numeric. Reserved country parameter for compatibility with earlier
  simulation scripts.

- betaC4:

  Numeric. Reserved country parameter for compatibility with earlier
  simulation scripts.

## Value

A data frame with columns `s`, `L1`, `L2`, `country`, `treatment`, and
`y`.

## Examples

``` r
dat <- dgp_xtrial_data(n = 100, delta_c = 0.10)
head(dat)
#>   s          L1 L2 country treatment y
#> 1 1  0.86208648  1       2       sc2 0
#> 2 1 -0.24323674  0       1         a 1
#> 3 2 -1.20608719  1       1         b 1
#> 4 1  0.01917759  1       3       sc1 1
#> 5 1  0.02956075  0       1         a 1
#> 6 2 -0.45017246  1       4         b 1
```
