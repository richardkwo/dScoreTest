# Weight for a fitted CATE model

Weight for a fitted CATE model

## Usage

``` r
weight_fun(fit, X, ...)

# S3 method for class 'CATE'
weight_fun(fit, X, ...)
```

## Arguments

- fit:

  A fitted model object.

- X:

  Covariate matrix.

- ...:

  Additional arguments passed to methods.

## Value

Numeric vector of weights.

## Methods (by class)

- `weight_fun(CATE)`: Weight for a `"CATE"` object: constant 1 for each
  row of `X`, since the outcome model has a squared-error loss.

## See also

`weight_fun.CATE`
