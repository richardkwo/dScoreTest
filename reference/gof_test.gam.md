# Goodness-of-fit test for a GAM

Debiased score test for goodness of fit of an
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) fit.

## Usage

``` r
# S3 method for class 'gam'
gof_test(
  object,
  hunt.style = "optimal",
  hunt.method = "grf",
  hunt_fun = NULL,
  trim.outlier.hunt = TRUE,
  X.cols.exclude = NULL,
  splits = c(0.5, 0.5),
  arg.hunt_fun = NULL,
  predict_fun_alt = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  Fitted [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) object.

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
  `arg.hunt_fun` and `predict_fun_alt` are used to specify a customized
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
  object `g` must support `predict_fun_alt(g, X)` for evaluation.

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

- predict_fun_alt:

  When a customized `hunt.fun` is used, this is a function with
  signature `predict_fun_alt(fit, X)` returning a numeric vector of
  predictions from a fitted alternative model produced by `hunt_fun()`.

- verbose:

  Default `FALSE`; information is printed if set to `TRUE`.

- ...:

  Unused; present for S3 generic/method consistency.

## Details

Only the numeric predictors appearing in `stats::model.frame(object)`
are exposed to the hunt; `X.cols.exclude` indexes into these predictor
variables (not basis columns). Factor-by smooths and other non-numeric
predictors are not currently supported. Formulas using
[`offset()`](https://rdrr.io/r/stats/offset.html) terms, a `weights`
argument, or a multi-column response (e.g. `cbind(succ, fail) ~ ...`)
are also not supported.

## Examples

``` r
 set.seed(42)
 dat <- mgcv::gamSim(eg=1, n=400, dist="normal", scale=2, verbose = FALSE)
 dat.0 <- dat[,1:5]
 
 # well-specified
 fit.0 <- mgcv::gam(y~s(x0)+s(x1)+s(x2)+s(x3),data=dat.0)
 test.0 <- gof_test(fit.0)
 # f3=0, also well-specified
 fit.1 <- mgcv::gam(y~s(x0)+s(x1)+s(x2),data=dat.0)
 test.1 <- gof_test(fit.1)
 plot(test.1)


 # misspecified
 dat.1 <- dat.0
 dat.1$y <- dat.1$y * dat$f0 
 fit.2 <- mgcv::gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat.1)
 test.2 <- gof_test(fit.2)
 plot(test.2)


```
