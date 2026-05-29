# Weighted-least-squares hunting

Returns the hunted signal h as a function

## Usage

``` r
hunt_wls(
  wls_alt_method,
  resids,
  X,
  X.cols = 1:ncol(X),
  trim.outlier = TRUE,
  arg.wls_alt_method = NULL,
  predict_fun_alt = stats::predict
)
```

## Arguments

- wls_alt_method:

  Function with signature `wls_alt_method(y, X, w, ...)` that returns a
  fitted *alternative model* \\\hat{g} \in \mathcal{G}\\ by minimizing
  \\\sum_i w_i (y_i - g(x_i))^2\\. The returned object must support
  `predict_fun_alt(g, X)` for evaluation.

- resids:

  Residuals (i.e., negative scores) of length n from the null model.

- X:

  Covariates of dim n x p.

- X.cols:

  Subset of covariates to hunt. (Default: `1:ncol(X)`)

- trim.outlier:

  If `TRUE`, outliers in \\\hat{h}(X)\\ will be trimmed from the hunted
  \\\hat{h}\\ using Tukey's IQR rule.

- arg.wls_alt_method:

  Named list of additional arguments passed to `wls_alt_method` (default
  to `NULL`).

- predict_fun_alt:

  Function with signature `predict_fun_alt(fit, X)` returning a numeric
  vector of predictions from the alternative-model fit. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html).

## Value

A function h with signature `h(X)`.
