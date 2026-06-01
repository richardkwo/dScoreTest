# Debiased score test: goodness-of-fit test and model comparison

Test whether a semiparametric (e.g., GAM) or parametric (e.g., glm)
regression model is well-specified. The test is a debiased
(Neyman-orthogonalized) score test computed via sample splitting: on a
held-out hunt sample, the null model is fit and a flexible ML algorithm
is used to hunt for a direction in which the null model's score seems
positive; on an independent test sample, that direction's score is
evaluated to assess the significance. The test employs orthogonalization
to eliminate plug-in bias from estimating the null model, so the
resulting test statistic is asymptotically standard normal under the
null without requiring a parametric form for the alternative.

## Usage

``` r
dScoreTest(
  y,
  X,
  score_fun,
  weight_fun,
  fit_method,
  wls_method,
  hunt.style = "optimal",
  hunt.method = "grf",
  hunt_fun = NULL,
  trim.outlier.hunt = TRUE,
  X.cols.hunt = 1:ncol(X),
  splits = c(0.5, 0.5),
  arg.fit_method = NULL,
  arg.wls_method = NULL,
  arg.hunt_fun = NULL,
  predict_fun = stats::predict,
  predict_fun_alt = NULL,
  verbose = FALSE
)
```

## Arguments

- y:

  Numeric response vector of length n.

- X:

  Numeric covariate matrix of dimension n x p.

- score_fun:

  Function with signature `score_fun(fit, y, X)` returning a vector of
  scores \\l'(\hat{f}(x_i), y_i)\\, which can be viewed as negative
  residuals.

- weight_fun:

  Function with signature `weight_fun(fit, X)` that computes the weight
  \\\mathbb{E}\[l''(\hat{f}(x_i), y_i) \| x_i\]\\ for each row \\x_i\\
  of X.

- fit_method:

  Function with signature `fit_method(y, X, ...)` that returns a fitted
  null model \\\hat{f} \in \mathcal{F}\\ by minimizing the loss \\\sum_i
  l(f(x_i), y_i)\\. For a fitted `f`, it must support
  `predict_fun(f, X)` for evaluation.

- wls_method:

  Function with signature `wls_method(y, X, w, ...)` that fits the null
  model \\\hat{f} \in \mathcal{F}\\ with weighted least squares, i.e.,
  minimizing \\\sum_i w_i (f(x_i) - y_i)^2\\. For a fitted `f`, it must
  support `predict_fun(f, X)` for evaluation.

- hunt.style:

  Hunting algorithm with the following options.

  - `'optimal'`: optimal hunting (default). See
    [`hunt_optimal`](https://unbiased.co.in/dScoreTest/reference/hunt_optimal.md).

  - `'wls'`: a simpler hunting using weighted least squares, which can
    be less powerful. See
    [`hunt_wls`](https://unbiased.co.in/dScoreTest/reference/hunt_wls.md).

  - `'vanilla'`: a basic hunting; not recommended unless unable to fit
    an alternative model with weighted least squares. See
    [`hunt_vanilla`](https://unbiased.co.in/dScoreTest/reference/hunt_vanilla.md).

- hunt.method:

  Built-in method for hunting. Currently available:

  - `'grf'`: regression forest from package `grf`.

  When this is set to any other value, arguments `hunt_fun` and
  `predict_fun_alt` must be set properly to supply a customized hunting
  method.

- hunt_fun:

  Default `NULL`. When `hunt.method` is not set to a built-in method,
  this is a customized function for hunting. When `hunt.style` is
  `'optimal'` or `'wls'`, this function must have signature
  `hunt_fun(y, X, w, ...)` that returns a fitted *alternative model*
  \\\hat{g} \in \mathcal{G}\\ via weighted least squares, i.e., by
  minimizing \\\sum_i w_i (y_i - g(x_i))^2\\; otherwise, for `'vanilla'`
  hunting, this function must have signature `hunt_fun(y, X, ...)` that
  returns an *alternative model* fitted in any fashion. The returned
  object `g` must support `predict_fun_alt(g, X)` for evaluation.

- trim.outlier.hunt:

  If `TRUE` (default), extreme values produced by the hunted function
  will be trimmed using Tukey's IQR rule.

- X.cols.hunt:

  Integer vector selecting which columns of `X` drive the hunt. Default
  `1:ncol(X)`. This is modified only in special settings, e.g., when
  there is an offset in the null model.

- splits:

  Numeric vector of length 2 or 3 giving the relative sizes of the
  sample splits; rescaled internally to sum to one. Default is
  `c(0.5, 0.5)`, which splits data into two halves for hunt and test
  respectively. Though typically unnecessary in practice, one can also
  specify a 3-way split for hunt, debiasing and test respectively.

- arg.fit_method:

  Named list of additional arguments passed to `fit_method` (default to
  `NULL`).

- arg.wls_method:

  Named list of additional arguments passed to `wls_method` (default to
  `NULL`).

- arg.hunt_fun:

  Extra arguments (default `NULL`) passed to the customized `hunt.fun`.

- predict_fun:

  Function with signature `predict_fun(fit, X)` returning a numeric
  vector of predictions from a fitted null model, which is produced by
  `fit_method()` and `wls_method()`. Note that if `fit` is \\\hat{f}\\,
  this function should return \\\hat{f}(X)\\. Default
  [`stats::predict`](https://rdrr.io/r/stats/predict.html). When y is
  binary, it must also support signature
  `predict_fun(fit, X, type='response')` for returning probabilities.

- predict_fun_alt:

  Default `NULL`. When `hunt.method` is not set to a built-in method,
  this is a function with signature `predict_fun_alt(fit, X)` returning
  a numeric vector of predictions from a fitted alternative model
  produced by `hunt_fun()`.

- verbose:

  Default `FALSE`; information is printed if set to `TRUE`.

## Value

An object of class `"dScoreTest"`: a list whose key elements are the
debiased test statistic `t.stat` and the one-sided p-value `p.val`
(right tail of the standard normal), along with the test-set score
residuals, the hunted direction, and the call. It has
[`print`](https://unbiased.co.in/dScoreTest/reference/print.dScoreTest.md),
[`summary`](https://unbiased.co.in/dScoreTest/reference/summary.dScoreTest.md)
and
[`plot`](https://unbiased.co.in/dScoreTest/reference/plot.dScoreTest.md)
methods.

## Details

For most scenarios, use one of these methods instead:

- Use
  [`gof_test`](https://unbiased.co.in/dScoreTest/reference/gof_test.md)
  to test whether a fitted model is well-specified against a
  nonparametric alternative. S3 methods are provided for `glm`
  ([`gof_test.glm`](https://unbiased.co.in/dScoreTest/reference/gof_test.glm.md)),
  `lm`
  ([`gof_test.lm`](https://unbiased.co.in/dScoreTest/reference/gof_test.lm.md))
  and [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)
  ([`gof_test.gam`](https://unbiased.co.in/dScoreTest/reference/gof_test.gam.md)).

- Use
  [`compare_models`](https://unbiased.co.in/dScoreTest/reference/compare_models.md)
  to test a null model `fit.0` against an alternative supermodel `fit.1`
  in the same model class. Similar to
  [`anova`](https://rdrr.io/r/stats/anova.html), method can be used to
  conduct a significance test of one or more predictors. In contrast
  with
  [`gof_test`](https://unbiased.co.in/dScoreTest/reference/gof_test.md),
  this method targets the alternative `fit.1`. S3 methods are provided
  for `glm`
  ([`compare_models.glm`](https://unbiased.co.in/dScoreTest/reference/compare_models.glm.md)),
  `lm`
  ([`compare_models.lm`](https://unbiased.co.in/dScoreTest/reference/compare_models.lm.md))
  and [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)
  ([`compare_models.gam`](https://unbiased.co.in/dScoreTest/reference/compare_models.gam.md)).

Use `dScoreTest` directly for full control over the score, weight, refit
and hunt routines: this is the underlying engine that the S3 methods
wrap.

## See also

Useful links:

- <https://unbiased.co.in/dScoreTest>

- <https://github.com/richardkwo/dScoreTest>

- Report bugs at <https://github.com/richardkwo/dScoreTest/issues>

[`plot.dScoreTest`](https://unbiased.co.in/dScoreTest/reference/plot.dScoreTest.md),
[`summary.dScoreTest`](https://unbiased.co.in/dScoreTest/reference/summary.dScoreTest.md),
[`hunt_optimal`](https://unbiased.co.in/dScoreTest/reference/hunt_optimal.md),
[`hunt_wls`](https://unbiased.co.in/dScoreTest/reference/hunt_wls.md),
[`hunt_vanilla`](https://unbiased.co.in/dScoreTest/reference/hunt_vanilla.md),
[`new_dScoreTest`](https://unbiased.co.in/dScoreTest/reference/new_dScoreTest.md)

## Author

**Maintainer**: F. Richard Guo <ricguo@umich.edu>
([ORCID](https://orcid.org/0000-0002-2081-7398)) \[copyright holder\]

Authors:

- Aditya Dhawan <ad950@cam.ac.uk>
