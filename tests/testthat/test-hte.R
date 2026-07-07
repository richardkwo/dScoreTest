# Tests for the CATE / HTE machinery in R/hte.R. grf-heavy, so skipped on CRAN.
# The functional test uses set.seed(42), under which the p-value is < 0.05.

.hte_data <- function(n = 200, p = 3, seed = 1) {
    set.seed(seed)
    Z  <- matrix(rnorm(n * p), n, p)
    Tr <- rbinom(n, 1, 0.5)
    y  <- Z[, 2] + Tr * (1 + Z[, 1]) + rnorm(n)
    list(y = y, Tr = Tr, Z = Z, X = cbind(Tr, Z))
}

test_that("fit_CATE and its predict/score/weight methods run", {
    skip_on_cran()
    d <- .hte_data()
    f <- fit_CATE(d$y, d$X, S = c(2, 3), folds.crossfit = 2)
    expect_s3_class(f, "CATE")
    expect_named(f, c("control_mean_fun", "CATE_fun", "S", "p"))
    yhat <- predict(f, d$X)
    expect_length(yhat, nrow(d$X))
    expect_true(all(is.finite(yhat)))
    expect_equal(score_fun(f, d$y, d$X), yhat - d$y, ignore_attr = TRUE)
    expect_equal(weight_fun(f, d$X), rep(1, nrow(d$X)))
    # S = full set gives a constant CATE
    f0 <- fit_CATE(d$y, d$X, S = 1:3, folds.crossfit = 2)
    expect_length(unique(f0$CATE_fun(d$Z)), 1)
})

test_that("hte_test_conditional runs for every hunt style", {
    skip_on_cran()
    d <- .hte_data(n = 180, p = 3)
    for (hs in c("optimal", "wls", "vanilla")) {
        tt <- hte_test_conditional(d$y, d$Tr, d$Z, S = 1:3, hunt.style = hs,
                                   folds.crossfit = 2)
        expect_s3_class(tt, "dScoreTest")
        expect_true(is.finite(tt$p.val) && tt$p.val >= 0 && tt$p.val <= 1)
    }
})

test_that("hte_test_conditional rejects strong heterogeneity (S = 1:p)", {
    skip_on_cran()
    set.seed(42)
    n <- 400; p <- 2
    Z  <- matrix(rnorm(n * p), n, p)
    Tr <- rbinom(n, 1, 0.5)                            # randomized -> easy nuisances
    y  <- Z[, 2] + Tr * (3 * Z[, 1]) + rnorm(n) * 0.5  # strong CATE heterogeneity
    tt <- hte_test_conditional(y, Tr, Z, S = 1:p, hunt.style = "optimal",
                               folds.crossfit = 2)
    expect_lt(tt$p.val, 0.05)
})
