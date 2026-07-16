# Fit the conditional treatment effect (CATE) function

Returns a fitted CATE \\\tau(Z)\\ where covariates \\Z_S\\ has no effect
modification given the remaining covariates. It employs R-loss and cross
fitting to estimate \$\$\mathbb{E}\[Y\|T,Z\] = \mu_0(Z) + T \cdot
\tau(Z), \quad \mu_0(Z) := \mathbb{E}\[Y \mid T=0, Z\].\$\$

## Usage

``` r
fit_CATE(y, X, w = rep(1, nrow(X)), S = 1:(ncol(X) - 1), folds.crossfit = 5)
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
