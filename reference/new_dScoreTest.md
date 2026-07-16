# Constructor for the debiased score test

Internal worker that builds a `dScoreTest` object from a fixed three-way
(or two-way) sample split. Called by
[`dScoreTest`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md).
Splits, hunting algorithm, and predict semantics are taken as fully
resolved arguments — no defaults are inferred from the data.

## Usage

``` r
new_dScoreTest(
  y,
  X,
  idx.hunt,
  idx.debias,
  idx.test,
  score_fun,
  weight_fun,
  fit_method,
  wls_method,
  hunt.style = "optimal",
  hunt.method = "customized",
  debias.method = "standard",
  debias_fun = debias_standard,
  fit_hunt_method = NULL,
  wls_hunt_method = NULL,
  X.cols.hunt = 1:ncol(X),
  binary.y = FALSE,
  trim.outlier.hunt = TRUE,
  predict_fun = stats::predict,
  predict_fun_hunt = stats::predict,
  arg.fit_method = NULL,
  arg.wls_method = NULL,
  arg.fit_hunt_method = NULL,
  arg.wls_hunt_method = NULL
)
```

## Arguments

- y:

  Numeric response vector of length n.

- X:

  Numeric covariate matrix of dimension n x p.

- idx.hunt:

  Integer indices into `1:n` for the hunting subsample.

- idx.debias:

  Integer indices into `1:n` for the debiasing subsample (refitting the
  null model and projecting the hunted direction).

- idx.test:

  Integer indices into `1:n` for the test subsample on which the test
  statistic is evaluated. May coincide with `idx.debias` (two-way split)
  or be disjoint (three-way split).

- score_fun:

  Function with signature `score_fun(fit, y, X)` returning a vector of
  scores \\l'(\hat{f}(x_i), y_i)\\.

- weight_fun:

  Function with signature `weight_fun(fit, X)` returning the weight
  \\\mathbb{E}\[l''(\hat{f}(x_i), y_i) \| x_i\]\\ for each row of `X`.

- fit_method:

  Function with signature `fit_method(y, X, ...)` returning a fitted
  null model. The returned object must support `predict_fun(fit, X)`.

- wls_method:

  Function with signature `wls_method(y, X, w, ...)` that fits the null
  model by weighted least squares. The returned object must support
  `predict_fun(fit, X)`.

- hunt.style:

  One of `"optimal"` (default), `"wls"`, or `"vanilla"`. Selects which
  `hunt_*` routine is used.

- hunt.method:

  String for the hunting method.

- debias.method:

  String for the debiasing method, recorded in the returned object's
  `Call`.

- debias_fun:

  Function performing the debiasing, with the same signature as
  [`debias_standard`](https://unbiased.co.in/dScoreTest/reference/debias_standard.md)
  (the default). Must return a list with an element `h`, the debiased
  hunted function.

- fit_hunt_method:

  Required when `hunt.style = "vanilla"`; ignored otherwise. Function
  with signature `fit_hunt_method(y, X, ...)` returning a fitted
  alternative model that supports `predict_fun_hunt(g, X)`.

- wls_hunt_method:

  Required when `hunt.style %in% c("optimal", "wls")`; ignored
  otherwise. Function with signature `wls_hunt_method(y, X, w, ...)`
  returning a fitted alternative model that supports
  `predict_fun_hunt(g, X)`.

- X.cols.hunt:

  Integer or name vector selecting which columns of `X` drive the hunt.
  Default: all columns.

- binary.y:

  Logical. When `TRUE`, the optimal hunter computes \\\mathrm{Var}(l' \|
  x)\\ from `predict_fun(fit, X, type='response')` assuming a Bernoulli
  response. Only consulted by `hunt.style = "optimal"`.

- trim.outlier.hunt:

  Logical. Passed to the chosen `hunt_*` routine as `trim.outlier`. If
  `TRUE` (default), extreme values in the hunted function will be
  removed using Tukey's IQR rule.

- predict_fun:

  Function with signature `predict_fun(fit, X)` returning predictions
  from a fitted null model. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html).

- predict_fun_hunt:

  Function with signature `predict_fun_hunt(fit, X)` returning
  predictions from a model fitted under alternative. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html).

- arg.fit_method, arg.wls_method, arg.fit_hunt_method,
  arg.wls_hunt_method:

  Named lists of additional arguments forwarded to the corresponding
  fitter via `do.call`. Default `NULL`.

## Value

A list of class `"dScoreTest"` with elements:

- `t.stat`:

  Debiased test statistic
  \\\sqrt{n\_{\mathrm{test}}}\\\bar{L}/\hat{\sigma}\_L\\.

- `p.val`:

  One-sided p-value (right tail of the standard normal).

- `resids`:

  Score residuals on the test subsample.

- `h`:

  Orthogonalised hunted direction on the test subsample.

- `h.raw`:

  Hunted direction before the outer debias projection.

- `hunted_fun`:

  The debiased hunted function \\\hat{h} - \hat{m}\_{\hat{h}}\\, a
  function that can be applied to X.

- `Data`:

  List with `X`, `y`, and the three index vectors.

- `Call`:

  Named list of methods, `hunt.style`, `hunt.method`, `debias.method`,
  both predict functions, and the four `arg.*` lists.

## See also

[`dScoreTest`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md),
[`hunt_optimal`](https://unbiased.co.in/dScoreTest/reference/hunt_optimal.md),
[`hunt_wls`](https://unbiased.co.in/dScoreTest/reference/hunt_wls.md),
[`hunt_vanilla`](https://unbiased.co.in/dScoreTest/reference/hunt_vanilla.md)
