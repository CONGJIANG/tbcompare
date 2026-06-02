# Get TB DGP parameters

Get TB DGP parameters

## Usage

``` r
get_tb_dgp_parameters(dgp_choice = "delta1")
```

## Arguments

- dgp_choice:

  Character. One of `"delta0"`, `"delta1"`, `"delta2"`, or `"delta3"`.

## Value

A list containing scenario-specific DGP parameters.

## Examples

``` r
get_tb_dgp_parameters("delta1")
#> $delta_c
#> [1] 0.1
#> 
```
