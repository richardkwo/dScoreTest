# Summary of the score test

Reports the headline statistic and p-value, the sample-split sizes, a
raw-vs-debiased comparison of the test statistic (so the effect of the
outer-projection debiasing step is visible), and
[`summary`](https://rdrr.io/r/base/summary.html) digests of the
diagnostic vectors `L`, `L.raw`, `h`, `h.raw` and `resids`.

## Usage

``` r
# S3 method for class 'dScoreTest'
summary(object, ...)

# S3 method for class 'summary.dScoreTest'
print(x, ...)
```

## Arguments

- object:

  A `dScoreTest` object.

- ...:

  Unused, for S3 consistency.

- x:

  A `summary.dScoreTest` object.

## Value

A list of class `"summary.dScoreTest"`.
