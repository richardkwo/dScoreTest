# Helpers ---------------------------------------------------------------------

.make_pair_data <- function(n = 200, seed = 1) {
    set.seed(seed)
    dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
    dat$y <- 1 + dat$x1 + dat$x3 + rnorm(n)
    dat
}

# Validation paths (cheap, no grf) --------------------------------------------

test_that("compare_models.default errors on unsupported class", {
    obj <- structure(list(), class = "weirdclass")
    expect_error(compare_models(obj), "not implemented")
})

test_that("compare_models.glm rejects non-nested fit.1", {
    dat <- .make_pair_data()
    fit.0   <- glm(y ~ x1 + x2,      family = gaussian(), data = dat)
    fit.bad <- glm(y ~ x2 + x3,      family = gaussian(), data = dat)  # drops x1
    expect_error(compare_models(fit.0, fit.bad), "supermodel")
})

test_that("compare_models.glm rejects non-glm fit.1", {
    dat <- .make_pair_data()
    fit.0   <- glm(y ~ x1 + x2,      family = gaussian(), data = dat)
    fit.bad <- list()
    expect_error(compare_models(fit.0, fit.bad))
})

test_that("compare_models.gam rejects non-nested fit.1", {
    set.seed(1)
    dat <- mgcv::gamSim(eg = 1, n = 100, dist = "normal",
                        scale = 2, verbose = FALSE)
    fit.0   <- mgcv::gam(y ~ s(x0) + s(x1),         data = dat)
    fit.bad <- mgcv::gam(y ~ s(x1) + s(x2),         data = dat)  # drops x0
    expect_error(compare_models(fit.0, fit.bad), "supermodel")
})

# End-to-end smoke tests (slow; skipped on CRAN) ------------------------------

test_that("compare_models.glm runs end-to-end and returns a dScoreTest", {
    skip_on_cran()
    dat <- .make_pair_data()
    fit.0 <- glm(y ~ x1 + x3,      family = gaussian(), data = dat)
    fit.1 <- glm(y ~ x1 + x2 + x3, family = gaussian(), data = dat)
    res <- compare_models(fit.0, fit.1)
    expect_s3_class(res, "dScoreTest")
    expect_true(is.finite(res$t.stat))
    expect_true(res$p.val >= 0 && res$p.val <= 1)
})

test_that("compare_models.lm dispatches via compare_models.glm", {
    skip_on_cran()
    dat <- .make_pair_data()
    fit.0 <- lm(y ~ x1 + x3,      data = dat)
    fit.1 <- lm(y ~ x1 + x2 + x3, data = dat)
    res <- compare_models(fit.0, fit.1)
    expect_s3_class(res, "dScoreTest")
    expect_true(is.finite(res$t.stat))
})

test_that("compare_models.gam runs end-to-end on a nested gam pair", {
    skip_on_cran()
    set.seed(1)
    dat <- mgcv::gamSim(eg = 1, n = 200, dist = "normal",
                        scale = 2, verbose = FALSE)
    dat <- dat[, 1:5]
    fit.0 <- mgcv::gam(y ~ s(x0) + s(x1) + s(x2),         data = dat)
    fit.1 <- mgcv::gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
    res <- compare_models(fit.0, fit.1)
    expect_s3_class(res, "dScoreTest")
    expect_true(is.finite(res$t.stat))
})

test_that("compare_models.glm respects hunt.style", {
    skip_on_cran()
    dat <- .make_pair_data()
    fit.0 <- glm(y ~ x1 + x3,      family = gaussian(), data = dat)
    fit.1 <- glm(y ~ x1 + x2 + x3, family = gaussian(), data = dat)
    for (style in c("optimal", "wls")) {
        res <- compare_models(fit.0, fit.1, hunt.style = style)
        expect_s3_class(res, "dScoreTest")
    }
    expect_error(compare_models(fit.0, fit.1, hunt.style = "bogus"))
})
