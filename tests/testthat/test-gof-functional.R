# Functional tests for gof_test methods. Each test sets set.seed(42) fresh
# (the user has verified the p-values fall on the expected side of the
# thresholds under this seed). All are skipped on CRAN because they exercise
# grf-heavy code paths.

# gof_test.glm ---------------------------------------------------------------

.gof_glm <- function(misspec, hunt.style) {
    set.seed(42)
    n <- 500
    X <- matrix(rnorm(n * 3), nrow = n)
    y0 <- 5 * exp(X[, 1] + X[, 3]) + rnorm(n) * 3
    y  <- if (misspec) y0 + exp(6 * cos(X[, 1] / 6)^2) / sqrt(n) else y0
    fit <- glm(y ~ X, family = gaussian(link = "log"), start = rep(1, 4))
    gof_test(fit, hunt.style = hunt.style)
}

test_that("gof_test.glm: well-specified, optimal -> p > 0.05", {
    skip_on_cran()
    expect_gt(.gof_glm(misspec = FALSE, "optimal")$p.val, 0.05)
})
test_that("gof_test.glm: well-specified, wls -> p > 0.05", {
    skip_on_cran()
    expect_gt(.gof_glm(misspec = FALSE, "wls")$p.val, 0.05)
})
test_that("gof_test.glm: mis-specified, optimal -> p < 0.05", {
    skip_on_cran()
    expect_lt(.gof_glm(misspec = TRUE, "optimal")$p.val, 0.05)
})
test_that("gof_test.glm: mis-specified, wls -> p < 0.05", {
    skip_on_cran()
    expect_lt(.gof_glm(misspec = TRUE, "wls")$p.val, 0.05)
})

# gof_test.lm ----------------------------------------------------------------

.gof_lm <- function(misspec, hunt.style) {
    set.seed(42)
    n <- 500
    X <- matrix(rnorm(n * 3), nrow = n)
    X[, 3] <- X[, 3] + X[, 1] + X[, 2] / 2
    y0 <- 1 + X %*% c(1, 1, 2) + rnorm(n)
    y  <- if (misspec) y0 + cos(X[, 1]) else y0
    fit <- lm(y ~ X)
    gof_test(fit, hunt.style = hunt.style)
}

test_that("gof_test.lm: well-specified, optimal -> p > 0.05", {
    skip_on_cran()
    expect_gt(.gof_lm(misspec = FALSE, "optimal")$p.val, 0.05)
})
test_that("gof_test.lm: well-specified, wls -> p > 0.05", {
    skip_on_cran()
    expect_gt(.gof_lm(misspec = FALSE, "wls")$p.val, 0.05)
})
test_that("gof_test.lm: mis-specified, optimal -> p < 0.05", {
    skip_on_cran()
    expect_lt(.gof_lm(misspec = TRUE, "optimal")$p.val, 0.05)
})
test_that("gof_test.lm: mis-specified, wls -> p < 0.05", {
    skip_on_cran()
    expect_lt(.gof_lm(misspec = TRUE, "wls")$p.val, 0.05)
})

# gof_test.gam ---------------------------------------------------------------

.gof_gam <- function(misspec, hunt.style) {
    set.seed(42)
    dat <- mgcv::gamSim(eg = 1, n = 400, dist = "normal",
                        scale = 2, verbose = FALSE)
    dat.0 <- dat[, 1:5]
    if (misspec) {
        dat.0$y <- dat.0$y * dat$f0
    }
    fit <- mgcv::gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat.0)
    gof_test(fit, hunt.style = hunt.style)
}

test_that("gof_test.gam: well-specified, optimal -> p > 0.05", {
    skip_on_cran()
    expect_gt(.gof_gam(misspec = FALSE, "optimal")$p.val, 0.05)
})
test_that("gof_test.gam: well-specified, wls -> p > 0.05", {
    skip_on_cran()
    expect_gt(.gof_gam(misspec = FALSE, "wls")$p.val, 0.05)
})
test_that("gof_test.gam: mis-specified, optimal -> p < 0.05", {
    skip_on_cran()
    expect_lt(.gof_gam(misspec = TRUE, "optimal")$p.val, 0.05)
})
test_that("gof_test.gam: mis-specified, wls -> p < 0.05", {
    skip_on_cran()
    expect_lt(.gof_gam(misspec = TRUE, "wls")$p.val, 0.05)
})
