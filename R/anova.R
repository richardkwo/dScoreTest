#' Compare two models using the debiased score test
#' 
#' @seealso \code{\link{compare_models.glm}}, \code{\link{compare_models.lm}}
#' 
#' @export
compare_models <- function(object, ...) {
    UseMethod("compare_models")
}

#' @export
compare_models.default <- function(object, ...) {
    stop("compare_models is not implemented for class '", class(object)[1], "'")
}

# GLM and LM -----

#' Compare two fitted GLM models
#'
#' Debiased score test of the null model \code{fit.0} against the alternative
#' \code{fit.1}. The hunt is restricted to the linear span of
#' \code{model.matrix(fit.1)}, fit via weighted least squares; the null
#' tangent space is then orthogonalized out per the chosen \code{hunt.style}.
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
#'  compare_models(fit.0, fit.1, hunt.style="wls")
#'  anova(fit.0, fit.1)
#'  
#'  # test a misspecified null model: should be rejected
#'  fit.00 <- glm(y ~ x2, family = gaussian(link = "log"),
#'               data = dat, start = rep(1, 2))
#'  compare_models(fit.00, fit.1)
#'  \donttest{plot(compare_models(fit.00, fit.1))}
#'  compare_models(fit.00, fit.1, hunt.style="wls")
#'  \donttest{plot(compare_models(fit.00, fit.1, hunt.style="wls"))}
#'  anova(fit.00, fit.1)
#'
compare_models.glm <- function(fit.0, fit.1,
                               hunt.style = "optimal",
                               trim.outlier.hunt = TRUE,
                               splits = c(0.5, 0.5),
                               verbose = FALSE) {
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
#'  compare_models(fit.00, fit.1, hunt.style="wls")
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
