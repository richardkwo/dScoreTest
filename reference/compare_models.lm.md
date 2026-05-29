# Compare two fitted linear models

Debiased score test of the null `fit.0` against the alternative `fit.1`.
Both are internally refit as Gaussian-family GLMs and the call is
dispatched to
[`compare_models.glm`](https://richardkwo.github.io/dScoreTest/reference/compare_models.glm.md).

## Usage

``` r
# S3 method for class 'lm'
compare_models(fit.0, fit.1, ...)
```

## Arguments

- fit.0:

  The null model as a fitted `lm` object.

- fit.1:

  The alternative model as a fitted `lm` object. Must be a supermodel of
  `fit.0`; see
  [`compare_models.glm`](https://richardkwo.github.io/dScoreTest/reference/compare_models.glm.md).

- ...:

  Additional arguments passed to
  [`compare_models.glm`](https://richardkwo.github.io/dScoreTest/reference/compare_models.glm.md).

## Examples

``` r
 set.seed(42)
 n <- 500
 dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
 dat$x3 <- dat$x3 + (dat$x1 + dat$x2) / 3
 dat$y  <- 1 + dat$x1 + 2 * dat$x3 + rnorm(n)
 fit.0 <- lm(y ~ x1 + x3,        data = dat)
 fit.1 <- lm(y ~ x1 + x2 + x3,   data = dat)
 
 # test fit.0 against fit.1: should not be rejected
 compare_models(fit.0, fit.1)
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), x1, x2, x3.
#> (hunt.style = optimal, hunt.method = glm)
#> n = 500, two-way split: hunt = 250, debias & test = 250
#> 
#> T = -0.3495, p-value = 0.636634
 anova(fit.0, fit.1)
#> Analysis of Variance Table
#> 
#> Model 1: y ~ x1 + x3
#> Model 2: y ~ x1 + x2 + x3
#>   Res.Df    RSS Df Sum of Sq      F Pr(>F)
#> 1    497 502.73                           
#> 2    496 502.72  1 0.0097927 0.0097 0.9217

 # misspecified model: should be rejected
 fit.00 <- lm(y ~ x1 + x2, data = dat)
 compare_models(fit.00, fit.1)
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), x1, x2, x3.
#> (hunt.style = optimal, hunt.method = glm)
#> n = 500, two-way split: hunt = 250, debias & test = 250
#> 
#> T = 10.5216, p-value = 3.43451e-26
 plot(compare_models(fit.00, fit.1))


 compare_models(fit.00, fit.1, hunt.style="wls")
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), x1, x2, x3.
#> (hunt.style = wls, hunt.method = glm)
#> n = 500, two-way split: hunt = 250, debias & test = 250
#> 
#> T = 10.4095, p-value = 1.12168e-25
 anova(fit.00, fit.1)
#> Analysis of Variance Table
#> 
#> Model 1: y ~ x1 + x2
#> Model 2: y ~ x1 + x2 + x3
#>   Res.Df     RSS Df Sum of Sq      F    Pr(>F)    
#> 1    497 2295.63                                  
#> 2    496  502.72  1    1792.9 1768.9 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
