# Predict from a fitted CATE model

Computes the outcome mean \\\mathbb{E}\[Y \| T, Z\] = \mu_0(Z) +
T\\\tau(Z)\\ for a `"CATE"` object returned by
[`fit_CATE`](https://unbiased.co.in/dScoreTest/reference/fit_CATE.md).

## Usage

``` r
# S3 method for class 'CATE'
predict(object, X, ...)
```

## Arguments

- object:

  A `"CATE"` object from
  [`fit_CATE`](https://unbiased.co.in/dScoreTest/reference/fit_CATE.md).

- X:

  Matrix \\X = \[T, Z\]\\ of dimension m x (p+1), where the first column
  is the binary treatment T and the remaining are the covariates Z.

- ...:

  Unused, for S3 consistency.

## Value

Numeric vector of length `nrow(X)` giving \\\mu_0(Z) + T\\\tau(Z)\\.

## See also

[`fit_CATE`](https://unbiased.co.in/dScoreTest/reference/fit_CATE.md)
