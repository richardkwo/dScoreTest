# Functional tests for compare_models methods. Each test sets set.seed(42)
# fresh (the user has verified the p-values fall on the expected side of the
# thresholds under this seed). All are skipped on CRAN because they exercise
# grf-heavy code paths.

# compare_models.glm ---------------------------------------------------------

.cmp_glm <- function(misspec, hunt.style) {
    set.seed(42)
    n <- 500
    dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
    dat$x3 <- dat$x3 + (dat$x1 + dat$x2) / 3
    dat$y  <- 5 * exp(dat$x1 + dat$x3) + rnorm(n) * 3
    fit.1 <- glm(y ~ x1 + x2 + x3, family = gaussian(link = "log"),
                 data = dat, start = rep(1, 4))
    fit.0 <- if (misspec) {
        glm(y ~ x2, family = gaussian(link = "log"),
            data = dat, start = rep(1, 2))
    } else {
        glm(y ~ x1 + x3, family = gaussian(link = "log"),
            data = dat, start = rep(1, 3))
    }
    compare_models(fit.0, fit.1, hunt.style = hunt.style)
}

test_that("compare_models.glm: well-specified, optimal -> p > 0.05", {
    skip_on_cran()
    expect_gt(.cmp_glm(misspec = FALSE, "optimal")$p.val, 0.05)
})
test_that("compare_models.glm: well-specified, wls -> p > 0.05", {
    skip_on_cran()
    expect_gt(.cmp_glm(misspec = FALSE, "wls")$p.val, 0.05)
})
test_that("compare_models.glm: mis-specified, optimal -> p < 0.05", {
    skip_on_cran()
    expect_lt(.cmp_glm(misspec = TRUE, "optimal")$p.val, 0.05)
})
test_that("compare_models.glm: mis-specified, wls -> p < 0.05", {
    skip_on_cran()
    expect_lt(.cmp_glm(misspec = TRUE, "wls")$p.val, 0.05)
})

# compare_models.lm ----------------------------------------------------------

.cmp_lm <- function(misspec, hunt.style) {
    set.seed(42)
    n <- 500
    dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
    dat$x3 <- dat$x3 + (dat$x1 + dat$x2) / 3
    dat$y  <- 1 + dat$x1 + 2 * dat$x3 + rnorm(n)
    fit.1 <- lm(y ~ x1 + x2 + x3, data = dat)
    fit.0 <- if (misspec) lm(y ~ x1 + x2, data = dat) else lm(y ~ x1 + x3, data = dat)
    compare_models(fit.0, fit.1, hunt.style = hunt.style)
}

test_that("compare_models.lm: well-specified, optimal -> p > 0.05", {
    skip_on_cran()
    expect_gt(.cmp_lm(misspec = FALSE, "optimal")$p.val, 0.05)
})
test_that("compare_models.lm: well-specified, wls -> p > 0.05", {
    skip_on_cran()
    expect_gt(.cmp_lm(misspec = FALSE, "wls")$p.val, 0.05)
})
test_that("compare_models.lm: mis-specified, optimal -> p < 0.05", {
    skip_on_cran()
    expect_lt(.cmp_lm(misspec = TRUE, "optimal")$p.val, 0.05)
})
test_that("compare_models.lm: mis-specified, wls -> p < 0.05", {
    skip_on_cran()
    expect_lt(.cmp_lm(misspec = TRUE, "wls")$p.val, 0.05)
})

# compare_models.gam ---------------------------------------------------------

.cmp_gam <- function(misspec, hunt.style) {
    set.seed(42)
    dat <- mgcv::gamSim(eg = 1, n = 500, dist = "normal",
                        scale = 1, verbose = FALSE)
    dat <- dat[, 1:5]
    fit.1 <- mgcv::gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
    fit.0 <- if (misspec) {
        mgcv::gam(y ~ s(x0) + s(x1) + s(x3), data = dat)
    } else {
        mgcv::gam(y ~ s(x0) + s(x1) + s(x2), data = dat)
    }
    compare_models(fit.0, fit.1, hunt.style = hunt.style)
}

test_that("compare_models.gam: well-specified, optimal -> p > 0.05", {
    skip_on_cran()
    expect_gt(.cmp_gam(misspec = FALSE, "optimal")$p.val, 0.05)
})
test_that("compare_models.gam: well-specified, wls -> p > 0.05", {
    skip_on_cran()
    expect_gt(.cmp_gam(misspec = FALSE, "wls")$p.val, 0.05)
})
test_that("compare_models.gam: mis-specified, optimal -> p < 0.05", {
    skip_on_cran()
    expect_lt(.cmp_gam(misspec = TRUE, "optimal")$p.val, 0.05)
})
test_that("compare_models.gam: mis-specified, wls -> p < 0.05", {
    skip_on_cran()
    expect_lt(.cmp_gam(misspec = TRUE, "wls")$p.val, 0.05)
})
