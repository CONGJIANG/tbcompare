# Simulate cross-trial comparison data

User-facing wrapper for generating cross-trial simulation data. This
function uses the scenario choices `"delta0"`, `"delta1"`, `"delta2"`,
and `"delta3"` to set the standard-care contrast.

## Usage

``` r
simulate_xtrial_data(n = 1000, dgp_choice = "delta1", ...)
```

## Arguments

- n:

  Integer. Sample size.

- dgp_choice:

  Character. One of `"delta0"`, `"delta1"`, `"delta2"`, or `"delta3"`.

- ...:

  Additional arguments passed to
  [`dgp_xtrial_data()`](https://CONGJIANG.github.io/xtrial/reference/dgp_xtrial_data.md).

## Value

A data frame with simulated cross-trial comparison data.

## Examples

``` r
dat <- simulate_xtrial_data(n = 1000, dgp_choice = "delta1")
table(dat$s, dat$treatment)
#>    
#>       a   b sc1 sc2
#>   1 255   0 103 122
#>   2   0 263 137 120
```
