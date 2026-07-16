# Weighted-least-squares hunting

Hunt by fitting residuals on X, trained by solving a weighted least
squares.

## Usage

``` r
hunt_wls(
  wls_hunt_method,
  resids,
  X,
  X.cols = 1:ncol(X),
  trim.outlier = TRUE,
  arg.wls_hunt_method = NULL,
  predict_fun_hunt = stats::predict
)
```

## Arguments

- wls_hunt_method:

  Function with signature `wls_hunt_method(y, X, w, ...)` that returns a
  fitted *alternative model* \\\hat{g} \in \mathcal{G}\\ by minimizing
  \\\sum_i w_i (y_i - g(x_i))^2\\. The returned object must support
  `predict_fun_hunt(g, X)` for evaluation.

- resids:

  Residuals (i.e., negative scores) of length n from the null model.

- X:

  Covariates of dim n x p.

- X.cols:

  Subset of covariates to hunt. (Default: `1:ncol(X)`)

- trim.outlier:

  If `TRUE`, outliers in \\\hat{h}(X)\\ will be trimmed from the hunted
  \\\hat{h}\\ using Tukey's IQR rule.

- arg.wls_hunt_method:

  Named list of additional arguments passed to `wls_hunt_method`
  (default to `NULL`).

- predict_fun_hunt:

  Function with signature `predict_fun_hunt(fit, X)` returning a numeric
  vector of predictions from the alternative-model fit. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html).

## Value

An object of class `"hunt"`, a list with elements:

- `hunt.fit`:

  The fitted hunt model produced by `wls_hunt_method`.

- `trim.bounds`:

  The Tukey IQR trimming bounds, or `c(-Inf, Inf)` when
  `trim.outlier = FALSE`.

- `predict_fun_hunt`:

  The prediction function for `hunt.fit`, as supplied.

- `X.cols`:

  The columns of `X` used for the hunt, as supplied.

- `h`:

  A function with signature `h(X)` giving the hunted signal.
