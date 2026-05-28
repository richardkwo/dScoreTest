#' Compare two models using the debiased score test
#'
#' @param fit.0 A fitted null-model object. Methods are provided for
#'   \code{glm}, \code{lm} and \code{mgcv::gam} fits.
#' @param ... Additional arguments passed to the dispatched method,
#'   notably the alternative-model fit \code{fit.1}.
#'
#' @seealso \code{\link{compare_models.glm}}, \code{\link{compare_models.lm}},
#'   \code{\link{compare_models.gam}},
#'   \code{\link{gof_test}}, \code{\link{dScoreTest}}
#'
#' @export
compare_models <- function(fit.0, ...) {
    UseMethod("compare_models")
}

#' @export
compare_models.default <- function(fit.0, ...) {
    stop("compare_models is not implemented for class '", class(fit.0)[1], "'")
}

# GLM and LM -----

#' Compare two fitted GLM models
#'
#' Debiased score test of the null model \code{fit.0} against the alternative
#' \code{fit.1}. GLM \code{fit.1} is used to hunt for signal that \code{fit.0} 
#' potentially misses.
#'
#' @param fit.0 The null model as a fitted \code{glm} object.
#' @param fit.1 The alternative model as a fitted \code{glm} object. It must be
#'   a supermodel of \code{fit.0}: every column of
#'   \code{stats::model.matrix(fit.0)} must appear (by name) in
#'   \code{stats::model.matrix(fit.1)}. \code{fit.0} and \code{fit.1} must be
#'   fit on the same rows in the same order.
#' @param hunt.style Hunting algorithm with the following options.
#'   \itemize{
#'   \item \code{'optimal'}: optimal hunting (default).
#'      See \code{\link{hunt_optimal}}.
#'   \item \code{'wls'}: a simpler hunting using weighted least squares,
#'      which can be less powerful. See \code{\link{hunt_wls}}.}
#' @param trim.outlier.hunt If \code{TRUE} (default),
#'   extreme values produced by the hunted function will be trimmed using
#'   Tukey's IQR rule.
#' @param splits Numeric vector of length 2 or 3 giving the relative sizes
#'   of the sample splits; rescaled internally to sum to one.
#'   Default is \code{c(0.5, 0.5)}, which splits data into two halves for
#'   hunt and test respectively. Though typically unnecessary in practice,
#'   one can also specify a 3-way split for hunt, debiasing and test
#'   respectively.
#' @param verbose Default \code{FALSE}; information is printed if set to
#'   \code{TRUE}.
#' @param ... Unused; present for S3 generic/method consistency.
#'
#' @export
#'
#' @examples
#'  set.seed(42)
#'  n <- 500
#'  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
#'  dat$x3 <- dat$x3 + (dat$x1 + dat$x2) / 3
#'  dat$y <- 5 * exp(dat$x1 + dat$x3) + rnorm(n) * 3
#'  fit.0 <- glm(y ~ x1 + x3, family = gaussian(link = "log"),
#'               data = dat, start = rep(1, 3))
#'  fit.1 <- glm(y ~ x1 + x2 + x3, family = gaussian(link = "log"),
#'               data = dat, start = rep(1, 4))
#'               
#'  # test fit.0 against fit.1: should not be rejected
#'  compare_models(fit.0, fit.1)
#'  \donttest{compare_models(fit.0, fit.1, hunt.style="wls")}
#'  anova(fit.0, fit.1)
#'  
#'  # test a misspecified null model: should be rejected
#'  fit.00 <- glm(y ~ x2, family = gaussian(link = "log"),
#'               data = dat, start = rep(1, 2))
#'  compare_models(fit.00, fit.1)
#'  \donttest{plot(compare_models(fit.00, fit.1))}
#'  \donttest{compare_models(fit.00, fit.1, hunt.style="wls")}
#'  \donttest{plot(compare_models(fit.00, fit.1, hunt.style="wls"))}
#'  anova(fit.00, fit.1)
#'
compare_models.glm <- function(fit.0, fit.1,
                               hunt.style = "optimal",
                               trim.outlier.hunt = TRUE,
                               splits = c(0.5, 0.5),
                               verbose = FALSE,
                               ...) {
    # extract data and check nesting
    stopifnot(isa(fit.1, c("glm", "lm")))
    X.0 <- stats::model.matrix(fit.0)
    X.1 <- stats::model.matrix(fit.1)
    y   <- fit.0$y
    cols.0 <- match(colnames(X.0), colnames(X.1))
    if (any(is.na(cols.0))) {
        stop("fit.1 must be a supermodel of fit.0: every column of ",
             "model.matrix(fit.0) must appear in model.matrix(fit.1)")
    }
    start.0 <- stats::coef(fit.0)
    hunt.style <- match.arg(hunt.style, c("optimal", "wls"))

    # null-model fitting (restrict X to fit.0's columns)
    fit_method <- function(y, X) {
        stats::glm(y ~ . - 1, family = fit.0$family,
                   data = as.data.frame(cbind(y, X[, cols.0, drop = FALSE])),
                   start = start.0)
    }
    wls_method <- function(y, X, w) {
        stats::glm(y ~ . - 1, family = stats::gaussian(),
                   weights = w,
                   data = as.data.frame(cbind(y, X[, cols.0, drop = FALSE])),
                   start = start.0)
    }

    # score, weight, prediction for the null
    score_fun <- function(fit, y, X) {
        y.hat   <- stats::predict(fit, newdata = as.data.frame(X),
                                  type = "response")
        eta.hat <- stats::predict(fit, newdata = as.data.frame(X),
                                  type = "link")
        # (y.hat - y) / V(mu) * (d mu / d eta)
        v        <- fit$family$variance(y.hat)
        dmu.deta <- fit$family$mu.eta(eta.hat)
        return((y.hat - y) / v * dmu.deta)
    }
    weight_fun <- function(fit, X) {
        y.hat   <- stats::predict(fit, newdata = as.data.frame(X),
                                  type = "response")
        eta.hat <- stats::predict(fit, newdata = as.data.frame(X),
                                  type = "link")
        # (d mu / d eta)^2 / V(mu)
        v        <- fit$family$variance(y.hat)
        dmu.deta <- fit$family$mu.eta(eta.hat)
        return(dmu.deta^2 / v)
    }
    predict_fun <- function(fit, X, ...) {
        stats::predict(fit, newdata = as.data.frame(X), ...)
    }

    # custom hunt: fit.1's linear span via WLS
    hunt_fun <- function(y, X, w, ...) {
        stats::glm(y ~ . - 1, family = stats::gaussian(),
                   weights = w,
                   data = as.data.frame(cbind(y, X)))
    }
    predict_fun_alt <- function(fit, X, ...) {
        stats::predict(fit, newdata = as.data.frame(X), ...)
    }

    # run the debiased score test on the full alternative design
    dScoreTest(y, X.1,
               score_fun, weight_fun, fit_method, wls_method,
               hunt.style = hunt.style,
               hunt.method = "glm",
               hunt_fun = hunt_fun,
               trim.outlier.hunt = trim.outlier.hunt,
               X.cols.hunt = 1:ncol(X.1),
               splits = splits,
               predict_fun = predict_fun,
               predict_fun_alt = predict_fun_alt,
               verbose = verbose)
}

#' Compare two fitted linear models
#'
#' Debiased score test of the null \code{fit.0} against the alternative
#' \code{fit.1}. Both are internally refit as Gaussian-family GLMs and the
#' call is dispatched to \code{\link{compare_models.glm}}.
#'
#' @param fit.0 The null model as a fitted \code{lm} object.
#' @param fit.1 The alternative model as a fitted \code{lm} object. Must be a
#'   supermodel of \code{fit.0}; see \code{\link{compare_models.glm}}.
#' @param ... Additional arguments passed to \code{\link{compare_models.glm}}.
#'
#' @export
#'
#' @examples
#'  set.seed(42)
#'  n <- 500
#'  dat <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
#'  dat$x3 <- dat$x3 + (dat$x1 + dat$x2) / 3
#'  dat$y  <- 1 + dat$x1 + 2 * dat$x3 + rnorm(n)
#'  fit.0 <- lm(y ~ x1 + x3,        data = dat)
#'  fit.1 <- lm(y ~ x1 + x2 + x3,   data = dat)
#'  
#'  # test fit.0 against fit.1: should not be rejected
#'  compare_models(fit.0, fit.1)
#'  anova(fit.0, fit.1)
#'
#'  # misspecified model: should be rejected
#'  fit.00 <- lm(y ~ x1 + x2, data = dat)
#'  compare_models(fit.00, fit.1)
#'  \donttest{plot(compare_models(fit.00, fit.1))}
#'  \donttest{compare_models(fit.00, fit.1, hunt.style="wls")}
#'  anova(fit.00, fit.1)
#'
compare_models.lm <- function(fit.0, fit.1, ...) {
    glm.0 <- stats::glm(stats::formula(fit.0),
                        family = stats::gaussian(),
                        data = stats::model.frame(fit.0))
    glm.1 <- stats::glm(stats::formula(fit.1),
                        family = stats::gaussian(),
                        data = stats::model.frame(fit.1))
    compare_models.glm(glm.0, glm.1, ...)
}

# GAM --------

#' Compare two fitted GAM models
#'
#' Debiased score test of the null \code{fit.0} against the alternative
#' \code{fit.1}, both fitted with \code{mgcv::gam}. GLM \code{fit.1} is 
#' used to hunt for signal that \code{fit.0} potentially misses.
#'
#' @param fit.0 The null model as a fitted \code{mgcv::gam} object.
#' @param fit.1 The alternative model as a fitted \code{mgcv::gam} object.
#'   Must be a supermodel of \code{fit.0}: every predictor variable in
#'   \code{fit.0}'s formula must also appear in \code{fit.1}'s. \code{fit.0}
#'   and \code{fit.1} must be fit on the same rows in the same order.
#' @inheritParams compare_models.glm
#'
#' @details
#' Nesting is checked by predictor-variable name only, not by basis span:
#' two models with the same predictors but different smooth specifications
#' (e.g. \code{s(x, k = 5)} vs \code{s(x, k = 20)}) will pass the check.
#' The test remains meaningful as long as \code{fit.1}'s class contains the
#' relevant alternative directions. Factor predictors, \code{offset()} terms,
#' \code{weights} arguments, and multi-column responses are not supported.
#'
#' @export
#'
#' @examples
#'  set.seed(42)
#'  dat <- mgcv::gamSim(eg=1, n=500, dist="normal", scale=1)
#'  dat <- dat[, 1:5]
#'  
#'  # test fit.0 against fit.1: well-specified (f3 = 0) and should not be rejected
#'  fit.0  <- mgcv::gam(y ~ s(x0) + s(x1) + s(x2), data = dat)
#'  fit.1  <- mgcv::gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
#'  compare_models(fit.0, fit.1)
#'  \donttest{plot(compare_models(fit.0, fit.1))}
#'  anova(fit.0, fit.1)
#'  
#'  # mis-specified model: drops f2 and should be rejected
#'  fit.00 <- mgcv::gam(y ~ s(x0) + s(x1) + s(x3), data = dat)
#'  compare_models(fit.00, fit.1)
#'  \donttest{plot(compare_models(fit.00, fit.1))}
#'  \donttest{compare_models(fit.00, fit.1, hunt.style="wls")}
#'  anova(fit.00, fit.1)
#'
compare_models.gam <- function(fit.0, fit.1,
                               hunt.style = "optimal",
                               trim.outlier.hunt = TRUE,
                               splits = c(0.5, 0.5),
                               verbose = FALSE,
                               ...) {
    stopifnot(isa(fit.1, c("gam", "glm", "lm")))
    hunt.style <- match.arg(hunt.style, c("optimal", "wls"))

    # extract data
    formula.0 <- stats::formula(fit.0)
    formula.1 <- stats::formula(fit.1)
    family.0  <- fit.0$family
    mf.0      <- stats::model.frame(fit.0)
    mf.1      <- stats::model.frame(fit.1)

    # response (from the alternative)
    resp.idx <- attr(stats::terms(fit.1), "response")
    if (resp.idx == 0) {
        stop("compare_models.gam requires a formula with a response variable")
    }
    y.name <- names(mf.1)[resp.idx]
    y      <- mf.1[[resp.idx]]
    if (!is.null(dim(y))) {
        stop("multi-column responses (e.g. cbind(succ, fail) ~ ...) ",
             "are not supported")
    }

    # offsets/weights guard for both
    if (any(c("(offset)", "(weights)") %in% c(names(mf.0), names(mf.1)))) {
        stop("compare_models.gam does not support offset() terms or a ",
             "weights argument in the fitted models")
    }

    # nesting check by predictor variable name
    pred.0 <- setdiff(names(mf.0), y.name)
    pred.1 <- setdiff(names(mf.1), y.name)
    if (!all(pred.0 %in% pred.1)) {
        stop("fit.1 must be a supermodel of fit.0: every predictor of ",
             "fit.0 must appear in fit.1")
    }

    # X = fit.1's predictors, used for both null refits and the hunt
    keep <- seq_along(mf.1) != resp.idx
    X    <- as.matrix(mf.1[, keep, drop = FALSE])

    # null-model fitting (formula.0); rebind formula env so `w` resolves
    fit_method <- function(y, X) {
        data <- as.data.frame(X)
        data[[y.name]] <- y
        fm <- formula.0
        environment(fm) <- environment()
        mgcv::gam(fm, family = family.0, data = data)
    }
    wls_method <- function(y, X, w) {
        data <- as.data.frame(X)
        data[[y.name]] <- y
        fm <- formula.0
        environment(fm) <- environment()
        mgcv::gam(fm, family = stats::gaussian(),
                  weights = w, data = data)
    }

    # score, weight, prediction for the null
    # NB: predict.gam returns a 1-D array; as.numeric() coerces to a vector
    # so downstream (notably grf inside hunt_optimal) accepts it.
    score_fun <- function(fit, y, X) {
        y.hat   <- as.numeric(mgcv::predict.gam(fit, newdata = as.data.frame(X),
                                                type = "response"))
        eta.hat <- as.numeric(mgcv::predict.gam(fit, newdata = as.data.frame(X),
                                                type = "link"))
        v        <- fit$family$variance(y.hat)
        dmu.deta <- fit$family$mu.eta(eta.hat)
        return((y.hat - y) / v * dmu.deta)
    }
    weight_fun <- function(fit, X) {
        y.hat   <- as.numeric(mgcv::predict.gam(fit, newdata = as.data.frame(X),
                                                type = "response"))
        eta.hat <- as.numeric(mgcv::predict.gam(fit, newdata = as.data.frame(X),
                                                type = "link"))
        v        <- fit$family$variance(y.hat)
        dmu.deta <- fit$family$mu.eta(eta.hat)
        return(dmu.deta^2 / v)
    }
    predict_fun <- function(fit, X, ...) {
        as.numeric(mgcv::predict.gam(fit, newdata = as.data.frame(X), ...))
    }

    # custom hunt: fit.1's gam class via WLS Gaussian gam
    hunt_fun <- function(y, X, w, ...) {
        data <- as.data.frame(X)
        data[[y.name]] <- y
        fm <- formula.1
        environment(fm) <- environment()
        mgcv::gam(fm, family = stats::gaussian(),
                  weights = w, data = data)
    }
    predict_fun_alt <- function(fit, X, ...) {
        as.numeric(mgcv::predict.gam(fit, newdata = as.data.frame(X), ...))
    }

    dScoreTest(y, X,
               score_fun, weight_fun, fit_method, wls_method,
               hunt.style = hunt.style,
               hunt.method = "gam",
               hunt_fun = hunt_fun,
               trim.outlier.hunt = trim.outlier.hunt,
               X.cols.hunt = 1:ncol(X),
               splits = splits,
               predict_fun = predict_fun,
               predict_fun_alt = predict_fun_alt,
               verbose = verbose)
}
