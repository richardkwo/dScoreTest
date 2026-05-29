# Goodness-of-fit test for a linear model

Debiased score test for goodness of fit of an `lm`. Internally refits
the model as a Gaussian-family GLM and dispatches to
[`gof_test.glm`](https://richardkwo.github.io/dScoreTest/reference/gof_test.glm.md).

## Usage

``` r
# S3 method for class 'lm'
gof_test(object, ...)
```

## Arguments

- object:

  Fitted `lm` object.

- ...:

  Additional arguments passed to
  [`gof_test.glm`](https://richardkwo.github.io/dScoreTest/reference/gof_test.glm.md).

## Examples

``` r
 set.seed(42)
 n <- 500
 X <- matrix(rnorm(n * 3), nrow = n)
 X[,3] <- X[,3] + X[,1] + X[,2] / 2
 y0 <- 1 + X %*% c(1,1,2) + rnorm(n)  # well-specified
 fit.0 <- lm(y0 ~ X)
 test.0 <- gof_test(fit.0)
 plot(test.0)


 y1 <- y0 + cos(X[,1])  # mis-specified
 fit.1 <- lm(y1 ~ X)
 test.1 <- gof_test(fit.1)
 plot(test.1)


 
```
