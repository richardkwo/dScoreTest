# Score for a fitted CATE model

Score for a fitted CATE model

## Usage

``` r
score_fun(fit, y, X, ...)

# S3 method for class 'CATE'
score_fun(fit, y, X, ...)
```

## Arguments

- fit:

  A fitted model object.

- y:

  Numeric response vector.

- X:

  Covariate matrix.

- ...:

  Additional arguments passed to methods.

## Value

Numeric vector of scores.

## Methods (by class)

- `score_fun(CATE)`: Score for a `"CATE"` object: the residual
  \\\hat{\mathbb{E}}\[Y \| T, Z\] - Y\\, i.e. predicted value minus `y`.

## See also

`score_fun.CATE`
