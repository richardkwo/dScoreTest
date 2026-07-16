# Debiased score test for goodness of fit

Debiased score test for goodness of fit

## Usage

``` r
gof_test(object, ...)
```

## Arguments

- object:

  A fitted model object. Methods are provided for `glm`, `lm` and
  [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) fits.

- ...:

  Additional arguments passed to the dispatched method.

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

## See also

[`gof_test.glm`](https://unbiased.co.in/dScoreTest/reference/gof_test.glm.md),
[`gof_test.lm`](https://unbiased.co.in/dScoreTest/reference/gof_test.lm.md),
[`gof_test.gam`](https://unbiased.co.in/dScoreTest/reference/gof_test.gam.md),
[`compare_models`](https://unbiased.co.in/dScoreTest/reference/compare_models.md),
[`dScoreTest`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md)

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
#> T = 9.0599, p-value = 6.53019e-20
```
