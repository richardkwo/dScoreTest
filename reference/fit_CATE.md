# Fit the conditional treatment effect (CATE) function

Returns a fitted CATE \\\tau(Z)\\ where covariates \\Z_S\\ has no effect
modification given the remaining covariates. It employs R-loss and cross
fitting to estimate \$\$\mathbb{E}\[Y\|T,Z\] = \mu_0(Z) + T \cdot
\tau(Z), \quad \mu_0(Z) := \mathbb{E}\[Y \mid T=0, Z\].\$\$

## Usage

``` r
fit_CATE(
  y,
  X,
  w = rep(1, nrow(X)),
  S = 1:(ncol(X) - 1),
  folds.crossfit = 5,
  randomized = FALSE
)
```

## Arguments

- y:

  Numeric response vector of length n.

- X:

  Matrix \\X = \[T, Z\]\\ of dimension n x (p+1), where the first column
  is the binary treatment T and the remaining are the covariates Z.

- w:

  Non-negative numeric weight vector of length n. Defaults to
  `rep(1, nrow(X))`. Applied when fitting the CATE (with weights
  \\\tilde{T}^2 w\\) and the control mean (with weights \\w\\). When it
  is not constant, this amounts to fitting CATE (or rather
  \\\mathbb{E}\[Y \| T, Z\]\\) using weighted least squares.

- S:

  A subset of `1:(ncol(X)-1)` such that `X[,-1][,S]` gives the
  covariates S which have no effect modification conditionally. When
  `S=1:(ncol(X)-1)`, CATE must be a constant; when `S=NULL`, CATE is
  unrestricted.

- folds.crossfit:

  An integer for the number of folds in cross fitting using the R-loss
  to estimate CATE. When it is 1, no cross fitting is used.

- randomized:

  If `FALSE` (default), the propensity \\e(Z) = \mathbb{E}\[T \mid Z\]\\
  is estimated by a cross-fitted
  [`grf::probability_forest`](https://rdrr.io/pkg/grf/man/probability_forest.html).
  If `TRUE`, \\T\\ is assumed to be randomized (independent of `Z`), so
  \\e(Z)\\ is taken to be the constant `mean(T)`, fitted upfront on the
  full sample without cross-fitting.

## Value

An object of class `"CATE"`:

- `control_mean_fun`:

  \\\mu_0(Z) = \mathbb{E}\[Y \| T=0, Z\]\\

- `CATE_fun`:

  \\\tau(Z) = \mathbb{E}\[Y \| T=1, Z\] - \mathbb{E}\[Y \| T=0, Z\]\\,
  which only depends on \\Z\\ through \\Z\_{S^c}\\.

- `S`:

  `S` as specified.

- `p`:

  Number of covariates, which equals `ncol(Z)`.

## Examples

``` r
## A randomized trial in which the treatment effect depends on Z1 only.
set.seed(2)
n <- 500
Z  <- matrix(rnorm(n * 2), n, 2, dimnames = list(NULL, c("Z1", "Z2")))
Tr <- rbinom(n, 1, 0.5)                       # randomized treatment
tau <- Z[, "Z1"]                              # true CATE
y  <- Z[, "Z1"] + Z[, "Z2"] + Tr * tau + rnorm(n)

# \donttest{
## S = NULL leaves the CATE unrestricted, so it may depend on Z1 and Z2.
fit <- fit_CATE(y, cbind(Tr, Z), S = NULL, randomized = TRUE,
                folds.crossfit = 2)
cor(fit$CATE_fun(Z), tau)      # close to 1: the true CATE is recovered
#> [1] 0.9059014
sd(fit$CATE_fun(Z))
#> [1] 1.105731

## S = 1 bars Z1 from modifying the effect. Because the true CATE depends
## on Z1 alone, the fitted CATE then collapses to nearly a constant.
fit0 <- fit_CATE(y, cbind(Tr, Z), S = 1, randomized = TRUE,
                 folds.crossfit = 2)
sd(fit0$CATE_fun(Z))           # much smaller than above
#> [1] 0.1190038

## predict() gives the fitted outcome mean mu0(Z) + T * tau(Z).
head(predict(fit, cbind(Tr, Z)))
#> [1] -2.7260819  0.6204713  1.6418427 -1.4231552  0.5249579  1.4561778
# }
```
