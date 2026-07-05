#' Vanilla hunting
#'
#' Hunt by fitting residuals on X. 
#'
#' @param fit_hunt_method Function with signature \code{fit_hunt_method(y, X, ...)}
#'   that returns a fitted \emph{alternative model} \eqn{\hat{g} \in \mathcal{G}}.
#'   For a fitted \code{g}, it must support \code{predict_fun_hunt(g, X)} for 
#'   evaluation. 
#' @param resids Residuals (i.e., negative scores) of length n from the null model.
#' @param X Covariates of dim n x p.
#' @param X.cols Subset of covariates to hunt. (Default: \code{1:ncol(X)})
#' @param trim.outlier If \code{TRUE}, outliers in \eqn{\hat{h}(X)} will be
#'   trimmed from the hunted \eqn{\hat{h}} using Tukey's IQR rule.
#' @param arg.fit_hunt_method Named list of additional arguments passed to
#'   \code{fit_hunt_method} (default to \code{NULL}).
#' @param predict_fun_hunt Function with signature \code{predict_fun_hunt(fit, X)}
#'   returning a numeric vector of predictions from the alternative-model fit.
#'   Default \code{stats::predict}.
#'
#' @return An object of class \code{"hunt"}, a list with elements:
#'   \describe{
#'     \item{\code{hunt.fit}}{The fitted hunt model produced by
#'       \code{fit_hunt_method}.}
#'     \item{\code{trim.bounds}}{The Tukey IQR trimming bounds, or
#'       \code{c(-Inf, Inf)} when \code{trim.outlier = FALSE}.}
#'     \item{\code{predict_fun_hunt}}{The prediction function for
#'       \code{hunt.fit}, as supplied.}
#'     \item{\code{X.cols}}{The columns of \code{X} used for the hunt, as
#'       supplied.}
#'     \item{\code{h}}{A function with signature \code{h(X)} giving the hunted
#'       signal.}
#'   }
#' @export
hunt_vanilla <- function(fit_hunt_method,
                        resids, X, X.cols=1:ncol(X),
                        trim.outlier=TRUE,
                        arg.fit_hunt_method=NULL,
                        predict_fun_hunt=stats::predict) {
    stopifnot(
        "num of obs. in resids and X do not match" = length(resids) == nrow(X),
        "X must be a matrix, array or data frame" = length(dim(X)) == 2
    )
    hunt.fit <- do.call(fit_hunt_method,
                        c(list(resids, X[, X.cols, drop = FALSE]),
                          arg.fit_hunt_method))
    if (trim.outlier) {
        trim.bounds <- tukey_iqr_bounds(predict_fun_hunt(hunt.fit, X[, X.cols, drop = FALSE]))
    } else {
        trim.bounds <- c(-Inf, Inf)
    }
    h <- function(.X) {
        .y <- predict_fun_hunt(hunt.fit, .X[, X.cols, drop = FALSE])
        return(pmin(pmax(.y, trim.bounds[1]), trim.bounds[2]))
    }
    out <- list(hunt.fit = hunt.fit, trim.bounds = trim.bounds,
                predict_fun_hunt = predict_fun_hunt, X.cols = X.cols, h = h)
    class(out) <- "hunt"
    return(out)
}

#' Weighted-least-squares hunting
#'
#' Hunt by fitting residuals on X, trained by solving a weighted least squares. 
#'
#' @param wls_hunt_method Function with signature \code{wls_hunt_method(y, X, w, ...)} that
#'   returns a fitted \emph{alternative model} \eqn{\hat{g} \in \mathcal{G}} by
#'   minimizing \eqn{\sum_i w_i (y_i - g(x_i))^2}. The returned object must
#'   support \code{predict_fun_hunt(g, X)} for evaluation. 
#' @param resids Residuals (i.e., negative scores) of length n from the null model.
#' @param X Covariates of dim n x p.
#' @param X.cols Subset of covariates to hunt. (Default: \code{1:ncol(X)})
#' @param trim.outlier If \code{TRUE}, outliers in \eqn{\hat{h}(X)} will be
#'   trimmed from the hunted \eqn{\hat{h}} using Tukey's IQR rule.
#' @param arg.wls_hunt_method Named list of additional arguments passed to
#'   \code{wls_hunt_method} (default to \code{NULL}).
#' @param predict_fun_hunt Function with signature \code{predict_fun_hunt(fit, X)}
#'   returning a numeric vector of predictions from the alternative-model fit.
#'   Default \code{stats::predict}.
#'
#' @return An object of class \code{"hunt"}, a list with elements:
#'   \describe{
#'     \item{\code{hunt.fit}}{The fitted hunt model produced by
#'       \code{wls_hunt_method}.}
#'     \item{\code{trim.bounds}}{The Tukey IQR trimming bounds, or
#'       \code{c(-Inf, Inf)} when \code{trim.outlier = FALSE}.}
#'     \item{\code{predict_fun_hunt}}{The prediction function for
#'       \code{hunt.fit}, as supplied.}
#'     \item{\code{X.cols}}{The columns of \code{X} used for the hunt, as
#'       supplied.}
#'     \item{\code{h}}{A function with signature \code{h(X)} giving the hunted
#'       signal.}
#'   }
#' @export
hunt_wls <- function(wls_hunt_method,
                         resids, X, X.cols=1:ncol(X),
                         trim.outlier=TRUE,
                         arg.wls_hunt_method=NULL,
                         predict_fun_hunt=stats::predict) {
    stopifnot(
        "num of obs. in resids and X do not match" = length(resids) == nrow(X),
        "X must be a matrix, array or data frame" = length(dim(X)) == 2
    )
    # Floor near-zero residuals to avoid division issues
    .idx <- abs(resids) < 1e-9
    resids[.idx] <- sign(resids[.idx]) * 1e-9
    # fit
    hunt.fit <- do.call(wls_hunt_method,
                        c(list(1/resids, X[, X.cols, drop = FALSE], resids^2),
                          arg.wls_hunt_method))
    if (trim.outlier) {
        trim.bounds <- tukey_iqr_bounds(predict_fun_hunt(hunt.fit, X[, X.cols, drop = FALSE]))
    } else {
        trim.bounds <- c(-Inf, Inf)
    }
    h <- function(.X) {
        .y <- predict_fun_hunt(hunt.fit, .X[, X.cols, drop = FALSE])
        return(pmin(pmax(.y, trim.bounds[1]), trim.bounds[2]))
    }
    out <- list(hunt.fit = hunt.fit, trim.bounds = trim.bounds,
                predict_fun_hunt = predict_fun_hunt, X.cols = X.cols, h = h)
    class(out) <- "hunt"
    return(out)
}


#' Optimal hunting
#'
#' Hunt by fitting residuals on X optimally, trained by solving a weighted 
#' least squares. The hunted function is also pre-debiased. 
#' 
#'
#' @param wls_hunt_method Function with signature \code{wls_hunt_method(y, X, w, ...)} that
#'   returns a fitted \emph{alternative model} \eqn{\hat{g} \in \mathcal{G}} by
#'   minimizing \eqn{\sum_i w_i (y_i - g(x_i))^2}. The returned object must
#'   support \code{predict_fun_hunt(g, X)} for evaluation.
#' @param wls_method Function with signature \code{wls_method(y, X, w, ...)}
#'   that fits the \emph{null model} \eqn{\hat{f} \in \mathcal{F}} with weighted least
#'   squares, i.e., minimizing \eqn{\sum_i w_i (f(x_i) - y_i)^2}. The returned
#'   object must support \code{predict_fun(f, X)} for evaluation. 
#' @param score_fun Function with signature \code{score_fun(fit, y, X)}
#'   returning a vector of scores \eqn{l'(\hat{f}(x_i), y_i)}.
#' @param weight_fun Function with signature \code{weight_fun(fit, X)}
#'   that computes the weight \eqn{\mathbb{E}[l''(\hat{f}(x_i), y_i) | x_i]} for each
#'   row \eqn{x_i} of X.
#' @param fit Fitted null model. Must support \code{predict_fun(fit, X)}.
#' @param y Response vector of length n.
#' @param X Covariates of dim n x p.
#' @param X.cols Subset of covariates to hunt. Default \code{1:ncol(X)}.
#' @param binary.y Set to \code{TRUE} only if y is binary (default: \code{FALSE}).
#'   This only affects how the variance function is estimated. When
#'   \code{TRUE}, \code{predict_fun(fit, X)} must return the conditional
#'   probability \eqn{P(y = 1 | x)}.
#' @param trim.outlier If \code{TRUE}, outliers in \eqn{\hat{h}(X)} will be
#'   trimmed from the hunted \eqn{\hat{h}} using Tukey's IQR rule.
#' @param arg.wls_hunt_method Named list of additional arguments passed to
#'   \code{wls_hunt_method} (default to \code{NULL}).
#' @param arg.wls_method Named list of additional arguments passed to
#'   \code{wls_method} (default to \code{NULL}).
#' @param predict_fun Function with signature \code{predict_fun(fit, X)}
#'   returning a numeric vector of predictions from null-model fits. 
#'   Default \code{stats::predict}.
#' @param predict_fun_hunt Function with signature \code{predict_fun_hunt(fit, X)}
#'   returning a numeric vector of predictions from the alternative-model fit.
#'   Default \code{stats::predict}.
#' 
#' @details
#' If y is binary, then set \code{binary.y=TRUE}. Meanwhile, 
#' \code{predict_fun(fit, X, type="response")} must output the predicted probabilities. 
#' 
#' @return An object of class \code{"hunt"}, a list with elements:
#'   \describe{
#'     \item{\code{hunt.fit}}{The fitted hunt model produced by
#'       \code{wls_hunt_method}.}
#'     \item{\code{trim.bounds}}{The Tukey IQR trimming bounds, or
#'       \code{c(-Inf, Inf)} when \code{trim.outlier = FALSE}.}
#'     \item{\code{predict_fun_hunt}}{The prediction function for
#'       \code{hunt.fit}, as supplied.}
#'     \item{\code{X.cols}}{The columns of \code{X} used for the hunt, as
#'       supplied.}
#'     \item{\code{h}}{A function with signature \code{h(X)} giving the
#'       orthogonalised hunted signal.}
#'   }
#'
#' @export
hunt_optimal <- function(wls_hunt_method, wls_method,
                     score_fun, weight_fun, fit, y, X, X.cols=1:ncol(X),
                     binary.y=FALSE,
                     trim.outlier=TRUE,
                     arg.wls_hunt_method=NULL, arg.wls_method=NULL,
                     predict_fun=stats::predict,
                     predict_fun_hunt=stats::predict) {
    stopifnot(
        "num of obs. in y and X do not match" = length(y) == nrow(X),
        "X must be a matrix, array or data frame" = length(dim(X)) == 2
    )
    # get resids
    resids <- -1 * score_fun(fit, y, X)
    # floor near-zero residuals to avoid division issues
    .idx <- abs(resids) < 1e-9
    resids[.idx] <- sign(resids[.idx]) * 1e-9
    # fit
    hunt.fit <- do.call(wls_hunt_method,
                 c(list(1/resids, X[, X.cols, drop = FALSE], resids^2),
                   arg.wls_hunt_method))
    # estimate the true variance function E[(l')^2 | x]
    if (binary.y) {
        var_fun <- function(.X) {
            p <- predict_fun(fit, .X, type="response")
            if (min(p) < 0 || max(p) > 1) {
                stop("predict_fun(fit,X,type='response') does not produce probabilities")
            }
            return(pmax(p * (1 - p), 1e-8))
        }
    } else {
        v.fit <- grf::regression_forest(X[, X.cols, drop = FALSE], resids^2,
                                        honesty = FALSE, tune.parameters = "all")
        var_fun <- function(.X) {
            v <- stats::predict(v.fit, .X[, X.cols, drop = FALSE])$predictions
            return(pmax(v, 1e-8))
        }
    }
    # get the null model projection with WLS
    w <- weight_fun(fit, X)
    v <- var_fun(X)
    h.raw <- predict_fun_hunt(hunt.fit, X[, X.cols, drop = FALSE])
    proj.fit <- do.call(wls_method,
                        c(list(v / w * h.raw, X[, X.cols, drop = FALSE], w^2 / v),
                          arg.wls_method))

    if (trim.outlier) {
        h.ortho <- h.raw - w / v * predict_fun(proj.fit, X[, X.cols, drop = FALSE])
        trim.bounds <- tukey_iqr_bounds(h.ortho)
    } else {
        trim.bounds <- c(-Inf, Inf)
    }
    h <- function(.X) {
        .w <- weight_fun(fit, .X)
        .v <- var_fun(.X)
        .h <- predict_fun_hunt(hunt.fit, .X[, X.cols, drop = FALSE])
        .proj <- predict_fun(proj.fit, .X[, X.cols, drop = FALSE])
        .y <- .h - .w / .v * .proj
        return(pmin(pmax(.y, trim.bounds[1]), trim.bounds[2]))
    }
    out <- list(hunt.fit = hunt.fit, trim.bounds = trim.bounds,
                predict_fun_hunt = predict_fun_hunt, X.cols = X.cols, h = h)
    class(out) <- "hunt"
    return(out)
}

# grf hunt functions -------
fit_hunt_method_grf <- function(y, X, ...) {
    wls_hunt_method_grf(y, X, NULL, ...)
}

wls_hunt_method_grf <- function(y, X, w, ...) {
    hunt.fit <- grf::regression_forest(X, y, sample.weights = w, ...)
    # prevent degenerating to a constant function 
    if (all(grf::variable_importance(hunt.fit)[,1] == 0)) {
        params <- hunt.fit$tunable.params
        # increase sample.fraction and refit
        params$sample.fraction <- max(params$sample.fraction, 
                                      params$min.node.size / length(y) * 2)
        hunt.fit <- grf::regression_forest(X, y, sample.weights = w, 
                              sample.fraction=params$sample.fraction,
                              mtry=params$mtry,
                              min.node.size=params$min.node.size,
                              alpha=params$alpha,
                              imbalance.penalty=params$imbalance.penalty, 
                              honesty=FALSE)
    } 
    return(hunt.fit)
}
arg.fit_hunt_method_grf <- list(honesty=FALSE, tune.parameters="all")
arg.wls_hunt_method_grf <- list(honesty=FALSE, tune.parameters="all")
predict_fun_hunt_grf <- function(.fit, .X) {
    stats::predict(.fit, .X)$predictions
}