# Getting started with dScoreTest

``` r

library(dScoreTest)
```

## What the test does

`dScoreTest` implements a **debiased (Neyman-orthogonalized) score
test**. Given a fitted null model, it asks: is there a direction in
which the model’s score is systematically non-zero? If so, the model is
misspecified.

The test is computed by **sample splitting**:

1.  On a held-out *hunt* sample, a flexible auxiliary fit searches for a
    promising direction `h(x)` of misspecification.
2.  On an independent *test* sample, the score along that direction is
    evaluated and standardized.

The orthogonalization step absorbs the plug-in bias from estimating `h`
on a finite sample, so the test statistic is asymptotically standard
normal under the null — without assuming a parametric form for the
alternative. The p-value is one-sided (power lives in the right tail).

## Two entry points

- [`gof_test()`](https://richardkwo.github.io/dScoreTest/reference/gof_test.md)
  — is a fitted model well-specified, against a nonparametric
  alternative?
- [`compare_models()`](https://richardkwo.github.io/dScoreTest/reference/compare_models.md)
  — does a nested alternative capture signal the null model misses? (An
  [`anova()`](https://rdrr.io/r/stats/anova.html)-style comparison.)

Both dispatch on the fitted object, with methods for `lm`, `glm`, and
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html). Both return a
`dScoreTest` object that supports
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## A first example

We use a small linear model so this vignette builds quickly. Take a
truth that is linear in `x1` and `x2`:

``` r

set.seed(3)
n <- 300
dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
dat$y <- 1 + dat$x1 + dat$x2 + rnorm(n)

fit <- lm(y ~ x1 + x2, data = dat)
gof_test(fit)
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), x1, x2.
#> (hunt.style = optimal, hunt.method = grf)
#> n = 300, two-way split: hunt = 150, debias & test = 150
#> 
#> T = -1.3137, p-value = 0.905519
```

The p-value is large: the linear model is correctly specified, and the
test does not reject.

Now introduce a quadratic effect the linear model cannot capture:

``` r

dat$y2 <- 1 + dat$x1 + dat$x1^2 + dat$x2 + rnorm(n)
fit.mis <- lm(y2 ~ x1 + x2, data = dat)
gof_test(fit.mis)
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), x1, x2.
#> (hunt.style = optimal, hunt.method = grf)
#> n = 300, two-way split: hunt = 150, debias & test = 150
#> 
#> T = 6.7450, p-value = 7.64993e-12
```

The p-value is small: the test detects that `E[y2 | x1, x2]` is not
linear in `x1`.

## Comparing nested models

[`compare_models()`](https://richardkwo.github.io/dScoreTest/reference/compare_models.md)
tests a null model against a richer alternative that contains it. We add
the quadratic term as an explicit column:

``` r

dat$x1sq <- dat$x1^2
fit.0 <- lm(y2 ~ x1 + x2,        data = dat)   # null (linear)
fit.1 <- lm(y2 ~ x1 + x1sq + x2, data = dat)   # alternative (superset)
compare_models(fit.0, fit.1)
#> Debiased score test: 
#> y ~ X, with X consists of (Intercept), x1, x1sq, x2.
#> (hunt.style = optimal, hunt.method = glm)
#> n = 300, two-way split: hunt = 150, debias & test = 150
#> 
#> T = 5.3831, p-value = 3.65976e-08
```

The alternative’s quadratic term captures the signal, so the comparison
rejects.

([`compare_models()`](https://richardkwo.github.io/dScoreTest/reference/compare_models.md)
refits the models from their model frames, so supply extra terms as
plain columns — e.g. `x1sq` above — rather than as in-formula
transformations like `poly(x1, 2)` or `I(x1^2)`.)

## Choosing the hunt

The search for a direction of misspecification (the “hunt”) has three
styles, passed via `hunt.style`:

- `"optimal"` (default) — the asymptotically optimal direction.
- `"wls"` — a simpler weighted-least-squares hunt; can be less powerful.
- `"vanilla"` — a basic hunt; a fallback.

## Where to go next

The default hunt for
[`gof_test()`](https://richardkwo.github.io/dScoreTest/reference/gof_test.md)
uses a regression forest (`grf`), which makes it a genuinely
nonparametric goodness-of-fit test. See the *Goodness of fit and model
comparison for GAMs* article for worked examples with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html), including
diagnostics via
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).
