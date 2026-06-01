# Plot the score test

Diagonotic plots for the test:

1.  Histogram of \\\\L_i\\\\, where \\L_i = resid_i \times h_i\\.

2.  \\\\L_i\\\\ against the index \\i\\, where \\i\\ refers to the
    \\i\\-th observation in the full dataset. Only those \\i\\'s in the
    test split are drawn. The mean is drawn as a horizontal line.
    Extremes values under the null can result in bad normal
    approximation. In this case, consider setting
    `trim.outlier.hunt=TRUE`.

3.  Residuals (negative scores) versus the hunted signal. A horizontal
    segment is drawn between each pair of raw hunted signal and the
    debiased hunted signal. If debiased gets higher, colored in red;
    otherwise colored in green. A regression line (blue) with a large
    positive slope indicates the model is misspecified.

4.  Normalized \\\\L_i\\\\ drawn in order.

## Usage

``` r
# S3 method for class 'dScoreTest'
plot(x, ...)
```

## Arguments

- x:

  A `dScoreTest` object.

- ...:

  Further graphical parameters passed to underlying plotting functions.

## Value

No return value; called for its side effect of producing the diagnostic
plots described above.
