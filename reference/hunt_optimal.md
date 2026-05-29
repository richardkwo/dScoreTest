# Optimal hunting

Returns the hunted signal h as a function

## Usage

``` r
hunt_optimal(
  wls_alt_method,
  wls_method,
  score_fun,
  weight_fun,
  fit,
  y,
  X,
  X.cols = 1:ncol(X),
  binary.y = FALSE,
  trim.outlier = TRUE,
  arg.wls_alt_method = NULL,
  arg.wls_method = NULL,
  predict_fun = stats::predict,
  predict_fun_alt = stats::predict
)
```

## Arguments

- wls_alt_method:

  Function with signature `wls_alt_method(y, X, w, ...)` that returns a
  fitted *alternative model* \\\hat{g} \in \mathcal{G}\\ by minimizing
  \\\sum_i w_i (y_i - g(x_i))^2\\. The returned object must support
  `predict_fun_alt(g, X)` for evaluation.

- wls_method:

  Function with signature `wls_method(y, X, w, ...)` that fits the *null
  model* \\\hat{f} \in \mathcal{F}\\ with weighted least squares, i.e.,
  minimizing \\\sum_i w_i (f(x_i) - y_i)^2\\. The returned object must
  support `predict_fun(f, X)` for evaluation.

- score_fun:

  Function with signature `score_fun(fit, y, X)` returning a vector of
  scores \\l'(\hat{f}(x_i), y_i)\\.

- weight_fun:

  Function with signature `weight_fun(fit, X)` that computes the weight
  \\\mathbb{E}\[l''(\hat{f}(x_i), y_i) \| x_i\]\\ for each row \\x_i\\
  of X.

- fit:

  Fitted null model. Must support `predict_fun(fit, X)`.

- y:

  Response vector of length n.

- X:

  Covariates of dim n x p.

- X.cols:

  Subset of covariates to hunt. Default `1:ncol(X)`.

- binary.y:

  Set to `TRUE` only if y is binary (default: `FALSE`). This only
  affects how the variance function is estimated. When `TRUE`,
  `predict_fun(fit, X)` must return the conditional probability \\P(y =
  1 \| x)\\.

- trim.outlier:

  If `TRUE`, outliers in \\\hat{h}(X)\\ will be trimmed from the hunted
  \\\hat{h}\\ using Tukey's IQR rule.

- arg.wls_alt_method:

  Named list of additional arguments passed to `wls_alt_method` (default
  to `NULL`).

- arg.wls_method:

  Named list of additional arguments passed to `wls_method` (default to
  `NULL`).

- predict_fun:

  Function with signature `predict_fun(fit, X)` returning a numeric
  vector of predictions from null-model fits. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html).

- predict_fun_alt:

  Function with signature `predict_fun_alt(fit, X)` returning a numeric
  vector of predictions from the alternative-model fit. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html).

## Value

A function h with signature `h(X)`.

## Details

If y is binary, then set `binary.y=TRUE`. Meanwhile,
`predict_fun(fit, X, type="response")` must output the predicted
probabilities.
