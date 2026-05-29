# Compare two models using the debiased score test

Compare two models using the debiased score test

## Usage

``` r
compare_models(fit.0, ...)
```

## Arguments

- fit.0:

  A fitted null-model object. Methods are provided for `glm`, `lm` and
  [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) fits.

- ...:

  Additional arguments passed to the dispatched method, notably the
  alternative-model fit `fit.1`.

## See also

[`compare_models.glm`](https://unbiased.co.in/dScoreTest/reference/compare_models.glm.md),
[`compare_models.lm`](https://unbiased.co.in/dScoreTest/reference/compare_models.lm.md),
[`compare_models.gam`](https://unbiased.co.in/dScoreTest/reference/compare_models.gam.md),
[`gof_test`](https://unbiased.co.in/dScoreTest/reference/gof_test.md),
[`dScoreTest`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md)
