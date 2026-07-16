#' Fit the conditional treatment effect (CATE) function 
#'
#' Returns a fitted CATE \eqn{\tau(Z)} where covariates \eqn{Z_S} has no effect 
#' modification given the remaining covariates. It employs R-loss and cross fitting 
#' to estimate
#' \deqn{\mathbb{E}[Y|T,Z] = \mu_0(Z) + T \cdot \tau(Z), \quad \mu_0(Z) := \mathbb{E}[Y \mid T=0, Z].}
#'
#' @param y Numeric response vector of length n.
#' @param X Matrix \eqn{X = [T, Z]} of dimension n x (p+1), where the first column 
#'      is the binary treatment T and the remaining are the covariates Z.
#' @param w Non-negative numeric weight vector of length n. Defaults to
#'      \code{rep(1, nrow(X))}. Applied when fitting the CATE (with weights
#'      \eqn{\tilde{T}^2 w}) and the control mean (with weights \eqn{w}). When 
#'      it is not constant, this amounts to fitting CATE (or rather \eqn{\mathbb{E}[Y | T, Z]})
#'      using weighted least squares. 
#' @param S A subset of \code{1:(ncol(X)-1)} such that \code{X[,-1][,S]} gives
#'      the covariates S which have no effect modification conditionally. 
#'      When \code{S=1:(ncol(X)-1)}, CATE must be a constant; when \code{S=NULL}, 
#'      CATE is unrestricted. 
#' @param folds.crossfit An integer for the number of folds in cross fitting
#'      using the R-loss to estimate CATE. When it is 1, no cross fitting is used.
#' 
#' @return An object of class \code{"CATE"}: 
#'      \describe{
#'      \item{\code{control_mean_fun}}{\eqn{\mu_0(Z) = \mathbb{E}[Y | T=0, Z]}}
#'      \item{\code{CATE_fun}}{\eqn{\tau(Z) = \mathbb{E}[Y | T=1, Z] - \mathbb{E}[Y | T=0, Z]}, 
#'      which only depends on \eqn{Z} through \eqn{Z_{S^c}}. }
#'      \item{\code{S}}{\code{S} as specified.}
#'      \item{\code{p}}{Number of covariates, which equals \code{ncol(Z)}.}
#'      }
#'      
#' @export
fit_CATE <- function(y, X, w=rep(1, nrow(X)), 
                     S=1:(ncol(X)-1), 
                     folds.crossfit=5) {
    # check input
    X <- as.matrix(X)
    p <- ncol(X) - 1
    stopifnot(
        "y must be numeric"                       = is.numeric(y),
        "number of obs. in X and y do not match." = length(y) == nrow(X),
        "w must be non-negative and of length nrow(X)" =
            length(w) == nrow(X) && all(w >= 0),
        "X must have at least a treatment and one covariate column" = p >= 1,
        "the first column of X (treatment T) must be binary 0/1" =
            all(X[, 1] %in% c(0, 1)),
        "folds.crossfit must be a positive integer" =
            length(folds.crossfit) == 1 && folds.crossfit >= 1 &&
            folds.crossfit == round(folds.crossfit)
    )
    if (!is.null(S)) {
        S <- as.integer(S)
        stopifnot(
            "S must be a subset of 1:(ncol(X)-1)" =
                all(S %in% seq_len(p)) && !anyDuplicated(S)
        )
    }
    if (folds.crossfit == 1) {
        warning("folds.crossfit = 1: nuisances are fit and evaluated ",
                "in-sample without cross-fitting, which overfits and can bias ",
                "the CATE. Use folds.crossfit >= 2.")
    }
    Tr <- X[, 1]
    Z  <- X[, -1, drop = FALSE]
    n  <- length(y)

    # cross-fitted nuisances Ytilde = Y - m(Z), Ttilde = T - e(Z), where
    # m(Z) = E[Y | Z] (regression forest) and e(Z) = E[T | Z] (probability
    # forest). With K = folds.crossfit folds, each observation is residualized
    # by nuisances trained on the other folds; when folds.crossfit == 1 the
    # nuisances are fit and evaluated in-sample (no cross-fitting).
    fold.id <- if (folds.crossfit == 1) {
        rep(1L, n)
    } else {
        sample(rep_len(seq_len(folds.crossfit), n))
    }
    Y.tilde <- numeric(n)
    T.tilde <- numeric(n)
    for (k in seq_len(folds.crossfit)) {
        if (folds.crossfit == 1) {
            train <- seq_len(n)
            test  <- seq_len(n)
        } else {
            train <- which(fold.id != k)
            test  <- which(fold.id == k)
        }
        m.fit <- grf::regression_forest(Z[train, , drop = FALSE], y[train],
                                        honesty = FALSE,
                                        tune.parameters = "all")
        e.fit <- grf::probability_forest(Z[train, , drop = FALSE],
                                         as.factor(Tr[train]), 
                                         honesty=FALSE)
        m.pred <- stats::predict(m.fit, Z[test, , drop = FALSE])$predictions
        e.pred <- stats::predict(e.fit, Z[test, , drop = FALSE])$predictions[, "1"]
        Y.tilde[test] <- y[test] - m.pred
        T.tilde[test] <- Tr[test] - e.pred
    }

    # fit CATE via the R-loss
    Sc <- setdiff(seq_len(p), S)
    if (length(Sc) == 0) {
        # S is the full set: CATE is constant, the R-learner weighted mean
        tau.const <- sum(w * T.tilde * Y.tilde) / sum(w * T.tilde^2)
        CATE_fun <- function(Z) {
            rep(tau.const, nrow(as.matrix(Z)))
        }
    } else {
        # regress Ytilde / Ttilde on Z_{S^c} with weights Ttilde^2 * w
        pseudo   <- Y.tilde / T.tilde
        cate.fit <- grf::regression_forest(Z[, Sc, drop = FALSE], pseudo,
                                           sample.weights = T.tilde^2 * w,
                                           honesty = FALSE,
                                           tune.parameters = "all")
        CATE_fun <- function(Z) {
            Z.new <- as.matrix(Z)
            stats::predict(cate.fit, Z.new[, Sc, drop = FALSE])$predictions
        }
    }

    # fit control_mean_fun mu_0(Z) = E[Y | T=0, Z] in sample:
    # regress Ycheck = Y - T * CATE(Z) on Z
    Y.check   <- y - Tr * CATE_fun(Z)
    mu0.fit   <- grf::regression_forest(Z, Y.check,
                                        sample.weights = w,
                                        honesty = FALSE,
                                        tune.parameters = "all")
    control_mean_fun <- function(Z) {
        stats::predict(mu0.fit, as.matrix(Z))$predictions
    }

    out <- list(control_mean_fun = control_mean_fun,
                CATE_fun = CATE_fun,
                S = S,
                p = p)
    class(out) <- "CATE"
    return(out)
}

#' Predict from a fitted CATE model
#'
#' Computes the outcome mean \eqn{\mathbb{E}[Y | T, Z] = \mu_0(Z) + T\,\tau(Z)}
#' for a \code{"CATE"} object returned by \code{\link{fit_CATE}}.
#'
#' @param object A \code{"CATE"} object from \code{\link{fit_CATE}}.
#' @param X Matrix \eqn{X = [T, Z]} of dimension m x (p+1), where the first
#'      column is the binary treatment T and the remaining are the covariates Z.
#' @param ... Unused, for S3 consistency.
#'
#' @return Numeric vector of length \code{nrow(X)} giving
#'      \eqn{\mu_0(Z) + T\,\tau(Z)}.
#'
#' @seealso \code{\link{fit_CATE}}
#' @export
predict.CATE <- function(object, X, ...) {
    X <- as.matrix(X)
    stopifnot(
        "X must have ncol equal to object$p + 1" = ncol(X) == object$p + 1
    )
    Tr <- X[, 1]
    Z  <- X[, -1, drop = FALSE]
    object$control_mean_fun(Z) + Tr * object$CATE_fun(Z)
}

#' Score for a fitted CATE model
#'
#' @param fit A fitted model object.
#' @param y Numeric response vector.
#' @param X Covariate matrix.
#' @param ... Additional arguments passed to methods.
#'
#' @return Numeric vector of scores.
#'
#' @seealso \code{\link{score_fun.CATE}}
#' @export
score_fun <- function(fit, y, X, ...) {
    UseMethod("score_fun")
}

#' @describeIn score_fun Score for a \code{"CATE"} object: the residual
#'   \eqn{\hat{\mathbb{E}}[Y | T, Z] - Y}, i.e. predicted value minus \code{y}.
#' @export
score_fun.CATE <- function(fit, y, X, ...) {
    predict.CATE(fit, X) - y
}

#' Weight for a fitted CATE model
#'
#' @param fit A fitted model object.
#' @param X Covariate matrix.
#' @param ... Additional arguments passed to methods.
#'
#' @return Numeric vector of weights.
#'
#' @seealso \code{\link{weight_fun.CATE}}
#' @export
weight_fun <- function(fit, X, ...) {
    UseMethod("weight_fun")
}

#' @describeIn weight_fun Weight for a \code{"CATE"} object: constant 1 for
#'   each row of \code{X}, since the outcome model has a squared-error loss.
#' @export
weight_fun.CATE <- function(fit, X, ...) {
    rep(1, nrow(X))
}

# Test -----

#' Test for detecting heterogeneity in the conditional treatment effect (CATE)
#'
#' For binary treatment \eqn{T}, covariates \eqn{Z} and outcome \eqn{Y}, it tests 
#' the hypothesis 
#' \deqn{H_0^c(S): \tau(Z) \text{ only depends on } Z \text{ through } Z_{S^c},}
#' where \eqn{S} is a subset of covariates and \eqn{\tau(Z)=\mathbb{E}[Y | T=1, Z] - \mathbb{E}[Y | T=0, Z]}
#' is the CATE. When \eqn{S} is the full set, this amounts to testing \eqn{\tau(Z)} 
#' is a constant; when \eqn{S} is a proper subset, this amounts to testing \eqn{Z_S}
#' has no further effect modification while holding \eqn{Z_{S^c}} fixed. 
#'
#' @param y Numeric response vector of length n.
#' @param Tr Binary (0/1) treatment vector of length n.
#' @param Z Numeric covariate matrix of dimension n x p.
#' @param S A subset of \code{1:ncol(Z)} giving the covariates whose effect
#'   modification is tested (they have none under the null, given the rest).
#'   Default \code{1:ncol(Z)} tests for \emph{any} heterogeneity (constant
#'   CATE under the null); \code{S = NULL} leaves the CATE unrestricted.
#' @param hunt.style One of \code{"optimal"} (default), \code{"wls"} or
#'   \code{"vanilla"}, selecting the hunting algorithm in
#'   \code{\link{dScoreTest}}. The hunted alternative is a \code{grf} outcome
#'   model fitted separately for the treated and control arms (T-learner).
#' @param folds.crossfit Number of cross-fitting folds passed to
#'   \code{\link{fit_CATE}} when estimating the CATE. Default 5.
#' @param trim.outlier.hunt,splits,verbose Passed through to
#'   \code{\link{dScoreTest}}; see there for details.
#' @param arg.hunt_grf Arguments passed to [grf::regression_forest()] 
#'   for hunting. 
#'
#' @inherit dScoreTest return
#'
#' @seealso \code{\link{fit_CATE}}, \code{\link{dScoreTest}}
#' @export
#'
#' @examples
#'  set.seed(1)
#'  n <- 600
#'  Z  <- matrix(rnorm(n * 3), n, 3)
#'  Tr <- rbinom(n, 1, plogis(Z[, 1]))
#'  y  <- Z[, 2] + Z[, 3] + Tr * (1 + Z[, 1]) + rnorm(n)  # CATE varies with Z1
#'  \donttest{
#'  # allow modification by Z1 (S = {2,3}): well-specified, should not reject
#'  hte_test_conditional(y, Tr, Z, S = c(2, 3))
#'  # forbid all modification (constant CATE): misspecified, should reject
#'  hte_test_conditional(y, Tr, Z, S = 1:3)
#'  }
hte_test_conditional <- function(y, Tr, Z, S=1:ncol(Z),
                                 hunt.style = "optimal",
                                 folds.crossfit=5,
                                 trim.outlier.hunt = TRUE,
                                 splits = c(0.5, 0.5),
                                 arg.hunt_grf = list(honesty=FALSE, tune.parameters="all"),
                                 verbose = FALSE) {
    hunt.style <- match.arg(hunt.style, c("optimal", "wls", "vanilla"))
    stopifnot("Tr must be binary" = setequal(unique(Tr), c(0,1)))
    # assemble the design X = [T, Z]
    Z <- as.matrix(Z)
    if (is.null(colnames(Z))) {
        colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    }
    if (is.null(S)) {
        warning("S = NULL: You are testing a vacuous hypothesis that is always TRUE.")
    } else {
        S <- as.integer(S)
        stopifnot(
            "S must be a subset of 1:ncol(Z)" =
                all(S %in% seq_len(ncol(Z))) && !anyDuplicated(S)
        )
        if (verbose && setequal(S, 1:ncol(Z))) {
            message("S = full set: You are testing the hypothesis that CATE(Z) = const.")
        }
    }
    X <- cbind(T = Tr, Z)

    # null-model fitters: fit_CATE with unit weights (fit_method) and with the
    # debiasing weights (wls_method). S is supplied via arg.fit_method and 
    # arg.wls_method.
    fit_method <- function(y, X, ...) {
        fit_CATE(y, X, ...)
    }
    wls_method <- function(y, X, w, ...) {
        fit_CATE(y, X, w = w, ...)
    }

    # Customized arm-specific grf hunt (fitted separately for T=0 and T=1).
    # dScoreTest routes hunt_fun to wls_hunt_method for optimal/wls hunting and
    # to fit_hunt_method for vanilla, so hunt_fun must match hunt.style.
    if (hunt.style == "vanilla") {
        hunt_fun <- fit_hunt_method_grf_hte_conditional
    } else {
        hunt_fun <- wls_hunt_method_grf_hte_conditional
    }

    dScoreTest(y, X,
               score_fun = score_fun.CATE,
               weight_fun = weight_fun.CATE,
               fit_method = fit_method,
               wls_method = wls_method,
               hunt.style = hunt.style,
               hunt.method = "grf.hte",
               hunt_fun = hunt_fun,
               debias.method= "hte.conditional", 
               debias_fun = debias_hte_conditional,
               trim.outlier.hunt = trim.outlier.hunt,
               splits = splits,
               arg.fit_method = list(S = S, folds.crossfit=folds.crossfit),
               arg.wls_method = list(S = S, folds.crossfit=folds.crossfit),
               arg.hunt_fun = arg.hunt_grf,
               predict_fun = predict.CATE,
               predict_fun_hunt = predict_fun_hunt_grf_hte_conditional,
               verbose = verbose)
}

# Hunt with grf T-learner --------
predict_fun_hunt_grf_hte_conditional <- function(.fit, X) {
    Tr.new <- X[, 1]
    Z.new  <- X[, -1, drop = FALSE]
    pred.0 <- stats::predict(.fit$mu_0.grf, Z.new)$predictions
    pred.1 <- stats::predict(.fit$mu_1.grf, Z.new)$predictions
    return(pred.0 + Tr.new * (pred.1 - pred.0))
}

wls_hunt_method_grf_hte_conditional <- function(y, X, w, ...) {
    Tr <- X[, 1]
    Z  <- X[, -1, drop = FALSE]
    if (!any(Tr == 0) || !any(Tr == 1)) {
        stop("hunt sample must contain both treated (T=1) and control (T=0) ",
             "units to fit the arm-specific hunt models.")
    }
    mu_0.grf <- wls_hunt_method_grf(y[Tr==0], Z[Tr==0, , drop = FALSE], w[Tr==0], ...)
    mu_1.grf <- wls_hunt_method_grf(y[Tr==1], Z[Tr==1, , drop = FALSE], w[Tr==1], ...)
    return(list(mu_0.grf = mu_0.grf, mu_1.grf = mu_1.grf))
}

fit_hunt_method_grf_hte_conditional <- function(y, X, ...) {
    wls_hunt_method_grf_hte_conditional(y, X, rep(1, nrow(X)), ...)
}

# Customized debiasing ------
#' Customized debiasing for [hte_test_conditional()]
#'
#' For a hunt \eqn{\hat{h}(t,z)}, it can 
#' be rewritten as \deqn{\hat{h}(t,z) = \hat{h}_0(z) + t \cdot \hat{h}_{\Delta}(z),} 
#' where \eqn{\hat{h}_{\Delta}(z) := \hat{h}(1,z) - \hat{h}(0,z)}. 
#' Let the projection of \eqn{\hat{h}(t,z)} onto the null space be
#' \deqn{m_{\hat{h}}(t,z) = m_0(z) + t \cdot m_{\Delta}(z_{S^c}).}
#' This function returns the debiased hunt function \eqn{\hat{h} - \hat{m}_{\hat{h}}}.
#' 
#' @details Function \eqn{\hat{m}_{\Delta}} is fitted by minimizing 
#' \deqn{\sum_i (T_i - \hat{e}(Z_i))^2 \, (\hat{h}_{\Delta}(Z_i) - m_{\Delta}(Z_{i,S^c}))^2,}
#' where \eqn{e(z):= \mathbb{E}(T \mid Z=z)}.
#' The debiased hunt function is given by
#' \deqn{(\hat{h} - \hat{m}_{\hat{h}})(t,z) = (t - \hat{e}(z))\,(\hat{h}_{\Delta}(z) - \hat{m}_{\Delta}(z_{S^c})).}
#'
#' @param h.hat Object of class \code{hunt} produced by [hunt_optimal()], 
#'      [hunt_wls()] or [hunt_vanilla()], where \code{h.hat$h(X)} is \eqn{\hat{h}(t, Z)} for 
#'      \eqn{X=[t, Z]}.
#' @param X.debias,fit.debias,predict_fun,weight_fun,wls_method,arg.wls_method See [dScoreTest()] for details.
#'      Argument \code{arg.wls_method} must contain field \code{S} that defines the model space.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{m.h.fit}}{The null model fitted (over all columns of \code{X})
#'     to project and debias \eqn{\hat{h}}.}
#'     \item{\code{h}}{The debiased hunt function \eqn{\hat{h} - \hat{m}_{\hat{h}}}
#'     with signature \code{h(X)}.}
#'   }
#' @keywords internal
debias_hte_conditional <- function(h.hat, X.debias, fit.debias,
                                   predict_fun, weight_fun, 
                                   wls_method, arg.wls_method) {
    # debias the CATE, fit with constant weight (square loss)
    Tr.debias <- X.debias[, 1]
    Z.debias  <- X.debias[, -1, drop = FALSE]
    n.debias <- nrow(X.debias)
    p <- ncol(Z.debias)
    S <- arg.wls_method$S
    stopifnot("Must have both T=1 and T=0 in the debiasing sample" = 
                  (sum(Tr.debias) > 0 && sum(Tr.debias) < n.debias))
    
    e.fit <- grf::probability_forest(Z.debias, as.factor(Tr.debias),
                                     honesty=FALSE)
    e.pred <- stats::predict(e.fit, Z.debias)$predictions[, "1"]
    # fit m_{\Delta} (i.e., CATE) on Z_{Sc}
    pred.h.delta <- h.hat$h(cbind(rep(1, n.debias), Z.debias)) - 
        h.hat$h(cbind(rep(0, n.debias), Z.debias))
    w <- (Tr.debias - e.pred)^2
    Sc <- setdiff(seq_len(p), S)
    if (length(Sc) == 0) {
        # constant function
        m.delta.fit <- sum(pred.h.delta * w) / sum(w)
        CATE_fun <- function(Z) {
            rep(m.delta.fit, nrow(as.matrix(Z)))
        }
    } else {
        # function of Sc
        m.delta.fit <- grf::regression_forest(Z.debias[, Sc, drop=FALSE], 
                                              pred.h.delta,
                                              sample.weights = w,
                                              honesty = FALSE,
                                              tune.parameters = "all")
        CATE_fun <- function(Z) {
            stats::predict(m.delta.fit, Z[, Sc, drop=FALSE])$predictions
        }
    }
    m.CATE.fit <- list(control_mean_fun=NULL, CATE_fun=CATE_fun)
    # debiased h = T * (h.CATE - m.h.CATE)
    h <- function(X.new) {
        Tr.new <- X.new[,1]
        Z.new <- X.new[, -1, drop = FALSE]
        n.new <- nrow(X.new)
        pred.h.CATE.new <- h.hat$h(cbind(rep(1, n.new), Z.new)) - h.hat$h(cbind(rep(0, n.new), Z.new))
        pred.m.h.CATE.new <- m.CATE.fit$CATE_fun(Z.new)
        pred.e.new <- stats::predict(e.fit, Z.new)$predictions[, "1"]
        return((Tr.new - pred.e.new) * (pred.h.CATE.new - pred.m.h.CATE.new))
    }
    return(list(m.h.fit = m.CATE.fit, h=h))
}