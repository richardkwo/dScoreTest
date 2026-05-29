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

[`compare_models.glm`](https://richardkwo.github.io/dScoreTest/reference/compare_models.glm.md),
[`compare_models.lm`](https://richardkwo.github.io/dScoreTest/reference/compare_models.lm.md),
[`compare_models.gam`](https://richardkwo.github.io/dScoreTest/reference/compare_models.gam.md),
[`gof_test`](https://richardkwo.github.io/dScoreTest/reference/gof_test.md),
[`dScoreTest`](https://richardkwo.github.io/dScoreTest/reference/dScoreTest.md)
