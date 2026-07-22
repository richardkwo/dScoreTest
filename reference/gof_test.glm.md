# Goodness-of-fit test for GLM

Debiased score test for goodness of fit of GLM.

## Usage

``` r
# S3 method for class 'glm'
gof_test(
  object,
  hunt.style = "optimal",
  hunt.method = "grf",
  hunt_fun = NULL,
  trim.outlier.hunt = TRUE,
  X.cols.exclude = NULL,
  splits = c(0.5, 0.5),
  arg.hunt_fun = NULL,
  predict_fun_hunt = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  Fitted glm object.

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

  When this is set to any other value, arguments `hunt_fun`,
  `arg.hunt_fun` and `predict_fun_hunt` are used to specify a customized
  hunting method.

- hunt_fun:

  Default `NULL`. When `hunt.method` is not set to a built-in method,
  this is a customized function for hunting. When `hunt.style` is
  `'optimal'` or `'wls'`, this function must have signature
  `hunt_fun(y, X, w, ...)` that returns a fitted *alternative model*
  \\\hat{g} \in \mathcal{G}\\ via weighted least squares, i.e., by
  minimizing \\\sum_i w_i (y_i - g(x_i))^2\\; otherwise, for `'vanilla'`
  hunting, this function must have signature `hunt_fun(y, X, ...)` that
  returns an *alternative model* fitted in any fashion. The returned
  object `g` must support `predict_fun_hunt(g, X)` for evaluation.

- trim.outlier.hunt:

  If `TRUE` (default), extreme values produced by the hunted function
  will be trimmed using Tukey's IQR rule.

- X.cols.exclude:

  Columns in `stats::model.matrix(object)` to be excluded when hunting
  for alternative signal. Default `NULL`.

- splits:

  Numeric vector of length 2 or 3 giving the relative sizes of the
  sample splits; rescaled internally to sum to one. Default is
  `c(0.5, 0.5)`, which splits data into two halves for hunt and test
  respectively. Though typically unnecessary in practice, one can also
  specify a 3-way split for hunt, debiasing and test respectively.

- arg.hunt_fun:

  Extra arguments (default `NULL`) passed to the customized `hunt.fun`.

- predict_fun_hunt:

  When a customized `hunt.fun` is used, this is a function with
  signature `predict_fun_hunt(fit, X)` returning a numeric vector of
  predictions from a fitted alternative model produced by `hunt_fun()`.

- verbose:

  Default `FALSE`; information is printed if set to `TRUE`.

- ...:

  Unused; present for S3 generic/method consistency.

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

## Examples

``` r
 set.seed(42)
 n <- 500
 X <- matrix(rnorm(n * 3), nrow = n)
 # log(E[y]) ~ X well-specified
 y0 <- 5 * exp(X[,1] + X[,3]) + rnorm(n) * 3
 fit.0 <- glm(y0 ~ X, family = gaussian(link = "log"), start=rep(1,4))
 gof_test(fit.0)
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), X1, X2, X3.
#> (hunt.style = optimal, hunt.method = grf, debias.method = standard)
#> n = 500, two-way split: hunt = 250, debias & test = 250
#> 
#> T = -0.2426, p-value = 0.595837
 # log(E[y]) ~ X misspecified
 y1 <- y0 + exp(6 * cos(X[,1]/6)^2) / sqrt(n)
 fit.1 <- glm(y1 ~ X, family = gaussian(link = "log"), start=rep(1,4))
 gof_test(fit.1)
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), X1, X2, X3.
#> (hunt.style = optimal, hunt.method = grf, debias.method = standard)
#> n = 500, two-way split: hunt = 250, debias & test = 250
#> 
#> T = 11.2969, p-value = 6.79367e-30
```
