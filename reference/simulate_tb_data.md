# Simulate TB comparison data

User-facing wrapper for generating TB cross-trial simulation data. This
function uses the scenario choices `"delta0"`, `"delta1"`, `"delta2"`,
and `"delta3"` to set the standard-care contrast.

## Usage

``` r
simulate_tb_data(n = 1000, dgp_choice = "delta1", seed = NULL, ...)
```

## Arguments

- n:

  Integer. Sample size.

- dgp_choice:

  Character. One of `"delta0"`, `"delta1"`, `"delta2"`, or `"delta3"`.

- seed:

  Optional integer random seed.

- ...:

  Additional arguments passed to
  [`dgp_tb_trial_data()`](https://CONGJIANG.github.io/tbcompare/reference/dgp_tb_trial_data.md).

## Value

A data frame with simulated TB comparison data.

## Examples

``` r
dat <- simulate_tb_data(n = 1000, dgp_choice = "delta1", seed = 1)
table(dat$s, dat$treatment)
#>    
#>       a   b sc1 sc2
#>   1 261   0 124 117
#>   2   0 262 109 127
```
