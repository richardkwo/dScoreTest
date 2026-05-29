# Package index

## Goodness-of-fit tests

Test whether a fitted model is well-specified against a nonparametric
alternative.

- [`gof_test()`](https://unbiased.co.in/dScoreTest/reference/gof_test.md)
  : Debiased score test for goodness of fit
- [`gof_test(`*`<glm>`*`)`](https://unbiased.co.in/dScoreTest/reference/gof_test.glm.md)
  : Goodness-of-fit test for GLM
- [`gof_test(`*`<lm>`*`)`](https://unbiased.co.in/dScoreTest/reference/gof_test.lm.md)
  : Goodness-of-fit test for a linear model
- [`gof_test(`*`<gam>`*`)`](https://unbiased.co.in/dScoreTest/reference/gof_test.gam.md)
  : Goodness-of-fit test for a GAM

## Model comparison

Test a null model against a nested alternative in the same model class.

- [`compare_models()`](https://unbiased.co.in/dScoreTest/reference/compare_models.md)
  : Compare two models using the debiased score test
- [`compare_models(`*`<glm>`*`)`](https://unbiased.co.in/dScoreTest/reference/compare_models.glm.md)
  : Compare two fitted GLM models
- [`compare_models(`*`<lm>`*`)`](https://unbiased.co.in/dScoreTest/reference/compare_models.lm.md)
  : Compare two fitted linear models
- [`compare_models(`*`<gam>`*`)`](https://unbiased.co.in/dScoreTest/reference/compare_models.gam.md)
  : Compare two fitted GAM models

## Engine and methods

The underlying debiased score test and its S3 methods.

- [`dScoreTest()`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md)
  : Debiased score test: goodness-of-fit test and model comparison
- [`new_dScoreTest()`](https://unbiased.co.in/dScoreTest/reference/new_dScoreTest.md)
  : Constructor for the debiased score test
- [`print(`*`<dScoreTest>`*`)`](https://unbiased.co.in/dScoreTest/reference/print.dScoreTest.md)
  : Print the score test
- [`plot(`*`<dScoreTest>`*`)`](https://unbiased.co.in/dScoreTest/reference/plot.dScoreTest.md)
  : Plot the score test
- [`summary(`*`<dScoreTest>`*`)`](https://unbiased.co.in/dScoreTest/reference/summary.dScoreTest.md)
  [`print(`*`<summary.dScoreTest>`*`)`](https://unbiased.co.in/dScoreTest/reference/summary.dScoreTest.md)
  : Summary of the score test

## Hunt algorithms

Routines that hunt for a direction of misspecification.

- [`hunt_optimal()`](https://unbiased.co.in/dScoreTest/reference/hunt_optimal.md)
  : Optimal hunting
- [`hunt_wls()`](https://unbiased.co.in/dScoreTest/reference/hunt_wls.md)
  : Weighted-least-squares hunting
- [`hunt_vanilla()`](https://unbiased.co.in/dScoreTest/reference/hunt_vanilla.md)
  : Vanilla hunting

## Package

- [`dScoreTest()`](https://unbiased.co.in/dScoreTest/reference/dScoreTest.md)
  : Debiased score test: goodness-of-fit test and model comparison
