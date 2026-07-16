# Vanilla hunting

Hunt by fitting residuals on X.

## Usage

``` r
hunt_vanilla(
  fit_hunt_method,
  resids,
  X,
  X.cols = 1:ncol(X),
  trim.outlier = TRUE,
  arg.fit_hunt_method = NULL,
  predict_fun_hunt = stats::predict
)
```

## Arguments

- fit_hunt_method:

  Function with signature `fit_hunt_method(y, X, ...)` that returns a
  fitted *alternative model* \\\hat{g} \in \mathcal{G}\\. For a fitted
  `g`, it must support `predict_fun_hunt(g, X)` for evaluation.

- resids:

  Residuals (i.e., negative scores) of length n from the null model.

- X:

  Covariates of dim n x p.

- X.cols:

  Subset of covariates to hunt. (Default: `1:ncol(X)`)

- trim.outlier:

  If `TRUE`, outliers in \\\hat{h}(X)\\ will be trimmed from the hunted
  \\\hat{h}\\ using Tukey's IQR rule.

- arg.fit_hunt_method:

  Named list of additional arguments passed to `fit_hunt_method`
  (default to `NULL`).

- predict_fun_hunt:

  Function with signature `predict_fun_hunt(fit, X)` returning a numeric
  vector of predictions from the alternative-model fit. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html).

## Value

An object of class `"hunt"`, a list with elements:

- `hunt.fit`:

  The fitted hunt model produced by `fit_hunt_method`.

- `trim.bounds`:

  The Tukey IQR trimming bounds, or `c(-Inf, Inf)` when
  `trim.outlier = FALSE`.

- `predict_fun_hunt`:

  The prediction function for `hunt.fit`, as supplied.

- `X.cols`:

  The columns of `X` used for the hunt, as supplied.

- `h`:

  A function with signature `h(X)` giving the hunted signal.
