# Customized debiasing for [`hte_test_conditional()`](https://unbiased.co.in/dScoreTest/reference/hte_test_conditional.md)

For a hunt \\\hat{h}(t,z)\\, it can be rewritten as \$\$\hat{h}(t,z) =
\hat{h}\_0(z) + t \cdot \hat{h}\_{\Delta}(z),\$\$ where
\\\hat{h}\_{\Delta}(z) := \hat{h}(1,z) - \hat{h}(0,z)\\. Let the
projection of \\\hat{h}(t,z)\\ onto the null space be
\$\$m\_{\hat{h}}(t,z) = m_0(z) + t \cdot m\_{\Delta}(z\_{S^c}).\$\$ This
function returns the debiased hunt function \\\hat{h} -
\hat{m}\_{\hat{h}}\\.

## Usage

``` r
debias_hte_conditional(
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

  Object of class `hunt` produced by
  [`hunt_optimal()`](https://unbiased.co.in/dScoreTest/reference/hunt_optimal.md),
  [`hunt_wls()`](https://unbiased.co.in/dScoreTest/reference/hunt_wls.md)
  or
  [`hunt_vanilla()`](https://unbiased.co.in/dScoreTest/reference/hunt_vanilla.md),
  where `h.hat$h(X)` is \\\hat{h}(t, Z)\\ for \\X=\[t, Z\]\\.

- X.debias, fit.debias, predict_fun, weight_fun, wls_method,
  arg.wls_method:

  See
  [`dScoreTest()`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md)
  for details. Argument `arg.wls_method` must contain field `S` that
  defines the model space.

## Value

A list with elements:

- `m.h.fit`:

  The null model fitted (over all columns of `X`) to project and debias
  \\\hat{h}\\.

- `h`:

  The debiased hunt function \\\hat{h} - \hat{m}\_{\hat{h}}\\ with
  signature `h(X)`.

## Details

Function \\\hat{m}\_{\Delta}\\ is fitted by minimizing \$\$\sum_i (T_i -
\hat{e}(Z_i))^2 \\ (\hat{h}\_{\Delta}(Z_i) -
m\_{\Delta}(Z\_{i,S^c}))^2,\$\$ where \\e(z):= \mathbb{E}(T \mid Z=z)\\.
The debiased hunt function is given by \$\$(\hat{h} -
\hat{m}\_{\hat{h}})(t,z) = (t - \hat{e}(z))\\(\hat{h}\_{\Delta}(z) -
\hat{m}\_{\Delta}(z\_{S^c})).\$\$
