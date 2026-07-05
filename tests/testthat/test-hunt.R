# tukey_iqr_bounds -----------------------------------------------------------

test_that("tukey_iqr_bounds returns Q1-1.5IQR and Q3+1.5IQR", {
    x <- 1:9                       # Q1 = 3, Q3 = 7, IQR = 4
    b <- dScoreTest:::tukey_iqr_bounds(x)
    expect_equal(unname(b), c(3 - 1.5 * 4, 7 + 1.5 * 4))
})

# hunt_vanilla with a stub fitter --------------------------------------------

# fitter returns an "object" (just X[,1]) and predict_fun_hunt re-extracts it
.stub_fit_hunt <- function(y, X) {
    out <- list(coef = 1, name = colnames(X)[1])
    class(out) <- "stubfit"
    out
}
.stub_predict_hunt <- function(fit, X) X[, 1]

test_that("hunt_vanilla returns a list with hunt.fit, trim.bounds, and h()", {
    set.seed(1)
    n <- 40
    X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
    resids <- rnorm(n)
    out <- hunt_vanilla(.stub_fit_hunt, resids, X,
                        trim.outlier = FALSE,
                        predict_fun_hunt = .stub_predict_hunt)
    expect_named(out, c("hunt.fit", "trim.bounds", "predict_fun_hunt", "X.cols", "h"))
    expect_s3_class(out$hunt.fit, "stubfit")
    expect_identical(out$predict_fun_hunt, .stub_predict_hunt)
    expect_equal(out$trim.bounds, c(-Inf, Inf))
    expect_true(is.function(out$h))
    # untrimmed: should reproduce X[,1] exactly
    expect_equal(out$h(X), X[, 1], ignore_attr = TRUE)
})

test_that("hunt_vanilla validates input shapes", {
    X <- matrix(rnorm(20), 10, 2)
    expect_error(
        hunt_vanilla(.stub_fit_hunt, resids = rnorm(5), X = X,
                     predict_fun_hunt = .stub_predict_hunt),
        "num of obs")
    expect_error(
        hunt_vanilla(.stub_fit_hunt, resids = rnorm(10), X = 1:10,
                     predict_fun_hunt = .stub_predict_hunt),
        "matrix")
})

test_that("hunt_vanilla trims outliers via Tukey IQR rule", {
    set.seed(1)
    n <- 60
    X <- cbind(x1 = c(rnorm(n - 1), 100), x2 = rnorm(n))  # x1 has an outlier
    resids <- rnorm(n)
    out <- hunt_vanilla(.stub_fit_hunt, resids, X,
                        trim.outlier = TRUE,
                        predict_fun_hunt = .stub_predict_hunt)
    bounds <- dScoreTest:::tukey_iqr_bounds(X[, 1])
    expect_equal(out$trim.bounds, bounds)
    pred   <- out$h(X)
    expect_true(all(pred >= bounds[1] & pred <= bounds[2]))
})

# hunt_wls residual flooring -------------------------------------------------

test_that("hunt_wls floors near-zero residuals without flipping sign", {
    # Patch wls_hunt_method to capture the y/w it receives so we can inspect
    # what happens to near-zero residuals.
    captured <- new.env()
    capture_fitter <- function(y, X, w) {
        captured$y <- y
        captured$w <- w
        .stub_fit_hunt(y, X)
    }
    n <- 20
    X <- cbind(x1 = rnorm(n), x2 = rnorm(n))
    resids <- c(rep(1e-12, 2), rnorm(n - 2))  # two near-zero residuals
    hunt_wls(capture_fitter, resids, X,
             trim.outlier = FALSE,
             predict_fun_hunt = .stub_predict_hunt)
    # 1/resids was passed as y; tiny residuals become 1 / (sign * 1e-9)
    # i.e. ±1e9, not Inf.
    expect_true(all(is.finite(captured$y)))
    expect_true(all(captured$w >= 0))
})
