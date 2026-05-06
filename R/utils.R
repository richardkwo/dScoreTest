tukey_iqr_bounds <- function(y) {
    q1          <- stats::quantile(y, 0.25)
    q3          <- stats::quantile(y, 0.75)
    iqr         <- q3 - q1
    lower.bound <- q1 - 1.5 * iqr
    upper.bound <- q3 + 1.5 * iqr
    return(c(lower.bound, upper.bound))
}
