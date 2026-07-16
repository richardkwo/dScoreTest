# Standard debiasing

Transform a hunted function \\\hat{h}\\ into a debiased function
\\\hat{h} - \hat{m}\_{\hat{h}}\\, where \\\hat{m}\_{\hat{h}}\\ is the
projection of \\\hat{h}\\ onto the null model.

## Usage

``` r
debias_standard(
  h.hat,
  X.debias,
  fit.debias,
  predict_fun,
  weight_fun,
  wls_method,
  arg.wls_method
)
```

## Arguments

- h.hat:

  A list as returned by one of
  [`hunt_optimal()`](https://unbiased.co.in/dScoreTest/reference/hunt_optimal.md),
  [`hunt_wls()`](https://unbiased.co.in/dScoreTest/reference/hunt_wls.md),
  [`hunt_vanilla()`](https://unbiased.co.in/dScoreTest/reference/hunt_vanilla.md).

- X.debias:

  Part of X for debiasing.

- fit.debias:

  Null model fitted on the debiasing sample of X and y.

- predict_fun, weight_fun, wls_method, arg.wls_method:

  They must be compatible with `fit.debias`; see
  [`dScoreTest()`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md)
  for details.

## Value

A list with elements:

- `m.h.fit`:

  The null model fitted (over all columns of `X`) to project and debias
  \\\hat{h}\\.

- `h`:

  The debiased hunt function \\\hat{h} - \hat{m}\_{\hat{h}}\\ with
  signature `h(X)`.

## Details

The projection \\\hat{m}\_{\hat{h}}\\ is obtained by fitting the null
model (via `wls_method`, weighted by `weight_fun(fit.debias, X.debias)`)
with the hunted values `h.hat$h(X.debias)` as response. This projection
uses **all** columns of `X`, even when the hunt itself is driven by only
a subset of covariates (`h.hat$X.cols`).
