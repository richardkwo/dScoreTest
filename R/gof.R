#' @export
gof_test <- function(object, ...) {
    UseMethod("gof_test")
}

#' @export
gof_test.default <- function(object, ...) {
    stop("dScoreTest is not implemented for class '", class(object)[1], "'")
}

# GLM --------

#' Goodness-of-fit test for GLM
#'
#' Debiased score test for goodness of fit of GLM.
#' 
#' @param object Fitted glm object. 
#' @param hunt.style Hunting algorithm with the following options.  
#'   \itemize{
#'   \item \code{'optimal'}: optimal hunting (default). 
#'      See \code{\link{hunt_optimal}}.
#'   \item \code{'wls'}: a simpler hunting using weighted least squares, 
#'      which can be less powerful. See \code{\link{hunt_wls}}.
#'   \item \code{'vanilla'}: a basic hunting; not 
#'      recommended unless unable to fit an alternative model with weighted 
#'      least squares.  See \code{\link{hunt_vanilla}}.}
#' @param hunt.method Built-in method for hunting. Currently available:
#'   \itemize{
#'   \item \code{'grf'}: regression forest from package \code{grf}.
#'   }
#'   When this is set to any other value, arguments \code{hunt_fun},
#'   \code{arg.hunt_fun} and \code{predict_fun_alt} are used to specify a
#'   customized hunting method.
#' @param hunt_fun Default \code{NULL}. 
#'   When \code{hunt.method} is not set to a built-in method, 
#'   this is a customized function for hunting. When \code{hunt.style} is 
#'   \code{'optimal'} or \code{'wls'}, this function must have signature 
#'   \code{hunt_fun(y, X, w, ...)} that returns a fitted 
#'   \emph{alternative model} \eqn{\hat{g} \in \mathcal{G}} via weighted least 
#'   squares, i.e., by minimizing \eqn{\sum_i w_i (y_i - g(x_i))^2}; 
#'   otherwise, for \code{'vanilla'} hunting,
#'   this function must have signature \code{hunt_fun(y, X, ...)} that 
#'   returns an \emph{alternative model} fitted in any fashion. 
#'   The returned object \code{g} must support \code{predict_fun_alt(g, X)} 
#'   for evaluation.
#' @param trim.outlier.hunt If \code{TRUE} (default), 
#'   extreme values produced by the hunted function will be trimmed using Tukey's 
#'   IQR rule. 
#' @param X.cols.exclude Columns in \code{stats::model.matrix(object)} to be 
#'   excluded when hunting for alternative signal. Default \code{NULL}.
#' @param splits Numeric vector of length 2 or 3 giving the relative sizes
#'   of the sample splits; rescaled internally to sum to one.
#'   Default is \code{c(0.5, 0.5)}, which splits data into two halves for
#'   hunt and test respectively. Though typically unnecessary in practice,
#'   one can also specify a 3-way split for hunt, debiasing and test respectively.
#' @param arg.hunt_fun Extra arguments (default \code{NULL}) passed to the 
#'   customized \code{hunt.fun}.
#' @param predict_fun_alt When a customized \code{hunt.fun} is used, this is a 
#' function with signature \code{predict_fun_alt(fit, X)} returning a numeric 
#' vector of predictions from a fitted alternative model produced by 
#' \code{hunt_fun()}.
#' @param verbose Default \code{FALSE}; information is printed if set to 
#'   \code{TRUE}.
#' 
#' @export
#' 
#' @examples
#'  set.seed(42)
#'  n <- 500
#'  X <- matrix(rnorm(n * 3), nrow = n)
#'  # log(E[y]) ~ X well-specified
#'  y0 <- 5 * exp(X[,1] + X[,3]) + rnorm(n) * 3
#'  fit.0 <- glm(y0 ~ X, family = gaussian(link = "log"), start=rep(1,4))
#'  gof_test(fit.0)
#'  # log(E[y]) ~ X misspecified
#'  y1 <- y0 + exp(6 * cos(X[,1]/6)^2) / sqrt(n)
#'  fit.1 <- glm(y1 ~ X, family = gaussian(link = "log"), start=rep(1,4))
#'  gof_test(fit.1)
#' 
gof_test.glm <- function(object, 
                         hunt.style = "optimal",
                         hunt.method = "grf", 
                         hunt_fun = NULL,
                         trim.outlier.hunt=TRUE,
                         X.cols.exclude=NULL,
                         splits=c(0.5, 0.5),
                         arg.hunt_fun=NULL,
                         predict_fun_alt=NULL, 
                         verbose=FALSE) {
    # extract data
    X <- stats::model.matrix(object)
    y <- object$y
    start <- stats::coef(object)

    # fit and wls
    fit_method <- function(y, X) {
        stats::glm(y ~ . - 1, family = object$family,
            data = as.data.frame(cbind(y, X)),
            start = start)
    }

    wls_method <- function(y, X, w) {
        stats::glm(y ~ . - 1, family = stats::gaussian(),
                   weights=w,
                   data = as.data.frame(cbind(y, X)),
                   start = start)
    }
    
    # score and weight
    score_fun <- function(fit, y, X) {
        y.hat <- stats::predict(fit, newdata = as.data.frame(X), 
                                type="response")
        eta.hat <- stats::predict(fit, newdata = as.data.frame(X), 
                                  type="link")
        # (y.hat - y) / V(mu) * (d mu / d eta)
        v <- fit$family$variance(y.hat)
        dmu.deta <- fit$family$mu.eta(eta.hat)
        return((y.hat - y) / v * dmu.deta)
    }

    weight_fun <- function(fit, X) {
        y.hat <- stats::predict(fit, newdata = as.data.frame(X), 
                                type="response")
        eta.hat <- stats::predict(fit, newdata = as.data.frame(X), 
                                  type="link")
        # (d mu / d eta)^2 / V(mu)
        v <- fit$family$variance(y.hat)
        dmu.deta <- fit$family$mu.eta(eta.hat)
        return(dmu.deta^2 / v)
    }
    
    # prediction 
    predict_fun <- function(fit, X, ...) {
        stats::predict(fit, newdata = as.data.frame(X), ...)
    }
    
    # test
    X.cols.hunt <- 1:ncol(X)
    if (!is.null(X.cols.exclude)) {
        if (is.character(X.cols.exclude)) {
            X.cols.exclude <- match(X.cols.exclude, colnames(X))}
        X.cols.hunt <- setdiff(X.cols.hunt, X.cols.exclude)
    }
    dScoreTest(y, X, 
               score_fun, weight_fun, fit_method, wls_method,
               hunt.style=hunt.style, 
               hunt.method=hunt.method, 
               hunt_fun=hunt_fun, 
               trim.outlier.hunt=trim.outlier.hunt, 
               X.cols.hunt=X.cols.hunt,
               splits=splits, 
               arg.hunt_fun=arg.hunt_fun,
               predict_fun=predict_fun,
               predict_fun_alt=predict_fun_alt,
               verbose=verbose)
}

# LM --------

#' Goodness-of-fit test for a linear model
#'
#' Debiased score test for goodness of fit of an \code{lm}. Internally
#' refits the model as a Gaussian-family GLM and dispatches to
#' \code{\link{gof_test.glm}}.
#'
#' @param object Fitted \code{lm} object.
#' @param ... Additional arguments passed to \code{\link{gof_test.glm}}.
#'
#' @export
#' 
#' @examples
#'  set.seed(42)
#'  n <- 500
#'  X <- matrix(rnorm(n * 3), nrow = n)
#'  X[,3] <- X[,3] + X[,1] + X[,2] / 2
#'  y0 <- 1 + X %*% c(1,1,2) + rnorm(n)  # well-specified
#'  fit.0 <- lm(y0 ~ X)
#'  test.0 <- gof_test(fit.0)
#'  \donttest{plot(test.0)}
#'  y1 <- y0 + cos(X[,1])  # mis-specified
#'  fit.1 <- lm(y1 ~ X)
#'  test.1 <- gof_test(fit.1)
#'  \donttest{plot(test.1)}
#'  
gof_test.lm <- function(object, ...) {
    glm_obj <- stats::glm(stats::formula(object),
                          family = stats::gaussian(),
                          data = stats::model.frame(object))
    gof_test.glm(glm_obj, ...)
}