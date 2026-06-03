# Get cross-trial DGP parameters

Get cross-trial DGP parameters

## Usage

``` r
get_xtrial_dgp_parameters(dgp_choice = "delta1")
```

## Arguments

- dgp_choice:

  Character. One of `"delta0"`, `"delta1"`, `"delta2"`, or `"delta3"`.

## Value

A list containing scenario-specific DGP parameters.

## Examples

``` r
get_xtrial_dgp_parameters("delta1")
#> $delta_c
#> [1] 0.1
#> 
```
