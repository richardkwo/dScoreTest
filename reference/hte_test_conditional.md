# Test for detecting heterogeneity in the conditional treatment effect (CATE)

For binary treatment \\T\\, covariates \\Z\\ and outcome \\Y\\, it tests
the hypothesis \$\$H_0^c(S): \tau(Z) \text{ only depends on } Z \text{
through } Z\_{S^c},\$\$ where \\S\\ is a subset of covariates and
\\\tau(Z)=\mathbb{E}\[Y \| T=1, Z\] - \mathbb{E}\[Y \| T=0, Z\]\\ is the
CATE. When \\S\\ is the full set, this amounts to testing \\\tau(Z)\\ is
a constant; when \\S\\ is a proper subset, this amounts to testing
\\Z_S\\ has no further effect modification while holding \\Z\_{S^c}\\
fixed.

## Usage

``` r
hte_test_conditional(
  y,
  Tr,
  Z,
  S = 1:ncol(Z),
  hunt.style = "optimal",
  folds.crossfit = 5,
  trim.outlier.hunt = TRUE,
  splits = c(0.5, 0.5),
  arg.hunt_grf = list(honesty = FALSE, tune.parameters = "all"),
  verbose = FALSE,
  randomized = FALSE
)
```

## Arguments

- y:

  Numeric response vector of length n.

- Tr:

  Binary (0/1) treatment vector of length n.

- Z:

  Numeric covariate matrix of dimension n x p.

- S:

  A subset of `1:ncol(Z)` giving the covariates whose effect
  modification is tested (they have none under the null, given the
  rest). Default `1:ncol(Z)` tests for *any* heterogeneity (constant
  CATE under the null); `S = NULL` leaves the CATE unrestricted.

- hunt.style:

  One of `"optimal"` (default), `"wls"` or `"vanilla"`, selecting the
  hunting algorithm in
  [`dScoreTest`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md).
  The hunted alternative is a `grf` outcome model fitted separately for
  the treated and control arms (T-learner).

- folds.crossfit:

  Number of cross-fitting folds passed to
  [`fit_CATE`](https://unbiased.co.in/dScoreTest/reference/fit_CATE.md)
  when estimating the CATE. Default 5.

- trim.outlier.hunt, splits, verbose:

  Passed through to
  [`dScoreTest`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md);
  see there for details.

- arg.hunt_grf:

  Arguments passed to
  [`grf::regression_forest()`](https://rdrr.io/pkg/grf/man/regression_forest.html)
  for hunting.

- randomized:

  If `FALSE` (default), the propensity \\e(Z) = \mathbb{E}\[T \mid Z\]\\
  used by
  [`fit_CATE`](https://unbiased.co.in/dScoreTest/reference/fit_CATE.md)
  and by the debiasing step is estimated with
  [`grf::probability_forest`](https://rdrr.io/pkg/grf/man/probability_forest.html).
  If `TRUE`, `T` is assumed randomized (independent of `Z`), so \\e(Z)\\
  is taken to be the constant `mean(T)`, fitted upfront without
  cross-fitting.

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

[`fit_CATE`](https://unbiased.co.in/dScoreTest/reference/fit_CATE.md),
[`dScoreTest`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md)

## Examples

``` r
 set.seed(1)
 n <- 600
 Z  <- matrix(rnorm(n * 3), n, 3)
 Tr <- rbinom(n, 1, plogis(Z[, 1]))
 y  <- Z[, 2] + Z[, 3] + Tr * (1 + Z[, 1]) + rnorm(n)  # CATE varies with Z1
 # \donttest{
 # allow modification by Z1 (S = {2,3}): well-specified, should not reject
 hte_test_conditional(y, Tr, Z, S = c(2, 3))
#> Debiased score test: 
#> y ~ X, with X consists of T, Z1, Z2, Z3.
#> (hunt.style = optimal, hunt.method = grf.hte, debias.method = hte.conditional)
#> n = 600, two-way split: hunt = 300, debias & test = 300
#> 
#> T = 0.5666, p-value = 0.285481
 # forbid all modification (constant CATE): misspecified, should reject
 hte_test_conditional(y, Tr, Z, S = 1:3)
#> Debiased score test: 
#> y ~ X, with X consists of T, Z1, Z2, Z3.
#> (hunt.style = optimal, hunt.method = grf.hte, debias.method = hte.conditional)
#> n = 600, two-way split: hunt = 300, debias & test = 300
#> 
#> T = 3.8425, p-value = 6.08863e-05
 # }
```
