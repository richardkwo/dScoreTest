#' Debiased Score Test for Quantile Generalised Additive Models
#'
#' Tests goodness of fit of a quantile GAM fitted via \code{qgam::qgam}.
#' The score function is the quantile check function derivative, and the
#' projection uses a conditional density estimate (via a quantile forest)
#' to compute the appropriate weights.
#'
#' @param y Numeric response vector of length n.
#' @param X Numeric covariate matrix or data frame of dimension n x p.
#'   If \code{formula = NULL}, columns are renamed to \code{X1, X2, ...}
#'   internally. If a formula is supplied, column names must match those
#'   used in the formula.
#' @param formula A \code{\link[stats]{formula}} or character string
#'   specifying the null qGAM model. If \code{NULL} (default), a full
#'   additive model \code{y ~ s(X1, k=k) + s(X2, k=k) + ...} is used.
#' @param qu Numeric quantile level in (0, 1). Default \code{0.5}.
#' @param k Integer. Basis dimension for smooth terms when
#'   \code{formula = NULL}. Default \code{20}.
#' @param method Smoothing parameter estimation method passed to
#'   \code{qgam::qgam} via \code{argGam}. Default \code{"REML"}.
#' @param nfolds Integer. Number of cross-validation folds used for
#'   conditional density estimation in the projection step. Default \code{5}.
#' @param hunt.style Character. Hunting style passed to
#'   \code{\link{hunt.nonpara}}. One of \code{"opt.WLS"} (default),
#'   \code{"WLS"}, or \code{"vanilla"}.
#' @param in.sample Logical. Passed to \code{\link{debiased.scoretest}}.
#'   Default \code{TRUE}.
#' @param return.hunted Logical. If \code{TRUE}, returns a list with
#'   \code{p.value} and \code{h.hat}. Default \code{FALSE}.
#' @param multi.split Logical. Passed to \code{\link{debiased.scoretest}}.
#'   Default \code{FALSE}.
#' @param n.splits,B,agg.fun Multi-split parameters. See
#'   \code{\link{debiased.scoretest}}.
#' @param ... Additional arguments passed to \code{grf::regression_forest}
#'   within the hunting step.
#'
#' @return A p-value (numeric scalar), or if \code{return.hunted = TRUE}
#'   and \code{multi.split = FALSE}, a list with \code{p.value} and
#'   \code{h.hat}.
#'
#' @seealso \code{\link{debiased.scoretest}}, \code{\link{hunt.nonpara}},
#'   \code{\link{dst.gam}}
#'
#' @export
dst.qgam <- function(y, X,
                     formula = NULL,
                     qu = 0.5,
                     k = 20,
                     method = "REML",
                     nfolds = 5,
                     hunt.style = "opt.WLS",
                     in.sample = TRUE,
                     return.hunted = FALSE,
                     multi.split = FALSE,
                     n.splits = 5,
                     B = 30,
                     agg.fun = function(x) mean(x, na.rm = TRUE),
                     ...) {

  reg.method <- function(y, X) {
    .qgam_fit(y, X, formula = formula, qu = qu, k = k, method = method)
  }

  score.fun <- .qgam_score

  proj.method <- function(fit, h.hat, y, X) {
    .qgam_proj(fit, h.hat, y, X, formula = formula, k = k, nfolds = nfolds)
  }

  # weight.fit for opt.WLS: conditional density f(q_qu(X) | X)
  weight.fit <- function(y, X, fit) {
    fitted.quantile <- as.numeric(stats::predict(fit, newdata = as.data.frame(X),
                                                 type = "response"))
    dens <- .qgam_cond_density(y, X, nfolds = nfolds)
    function(Xnew) {
      Xnew.df        <- as.data.frame(Xnew)
      fitted.q.new   <- as.numeric(stats::predict(fit, newdata = Xnew.df, type = "response"))
      as.numeric(dens$fhat(fitted.q.new, Xnew.df))
    }
  }

  # proj.hunt.method: weighted GAM projection for opt.WLS inside hunting
  proj.hunt.method <- function(h.vals, X, weights) {
    X <- as.data.frame(X)
    if (is.null(formula)) colnames(X) <- paste0("X", seq_len(ncol(X)))
    pfit <- .gam_fit(y = h.vals, X = X, formula = formula,
                     family = "gaussian", link = NULL,
                     k = k, method = method, weights = weights)
    function(Xnew) {
      Xnew <- as.data.frame(Xnew)
      if (is.null(formula)) colnames(Xnew) <- paste0("X", seq_len(ncol(Xnew)))
      as.numeric(stats::predict(pfit, newdata = Xnew, type = "response"))
    }
  }

  hunt.method <- function(resids, X, y, fit) {
    hunt.nonpara(resids, X, y, fit,
                 hunt.style        = hunt.style,
                 weight.fit        = weight.fit,
                 proj.hunt.method  = proj.hunt.method,
                 ...)
  }

  debiased.scoretest(y, X,
                     score.fun     = score.fun,
                     reg.method    = reg.method,
                     hunt.method   = hunt.method,
                     proj.method   = proj.method,
                     in.sample     = in.sample,
                     return.hunted = return.hunted,
                     multi.split   = multi.split,
                     n.splits      = n.splits,
                     B             = B,
                     agg.fun       = agg.fun)
}


# ---- internal helpers -------------------------------------------------------

# Fit a qGAM
.qgam_fit <- function(y, X, formula = NULL, qu = 0.5, k = 20, method = "REML") {
  X <- as.data.frame(X)
  y <- as.numeric(y)

  if (is.null(formula)) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
    cov.names   <- colnames(X)
    formula     <- stats::as.formula(
      paste("y ~", paste0("s(", cov.names, ", k = ", k, ")", collapse = " + "))
    )
  } else {
    formula <- stats::as.formula(formula, env = environment())
  }

  qgam::qgam(formula, data = data.frame(y, X),
             qu = qu, argGam = list(method = method))
}

# Score: indicator of y <= fitted quantile, shifted by qu
.qgam_score <- function(y, X, fit) {
  qu <- fit$family$getQu()
  as.numeric(y - stats::predict(fit, newdata = as.data.frame(X),
                                type = "response") <= 0) - qu
}

# proj.method: project h.hat onto null qGAM space using conditional density weights
.qgam_proj <- function(fit, h.hat, y, X, formula = NULL, k = 20, nfolds = 5) {
  X      <- as.data.frame(X)
  h.vals <- h.hat(X)

  fitted.q <- as.numeric(stats::predict(fit, newdata = X, type = "response"))
  dens     <- .qgam_cond_density(y, X, nfolds = nfolds)
  weights  <- as.numeric(dens$fhat(fitted.q, X))
  weights  <- pmax(weights, 1e-8)

  proj.fit <- .gam_fit(y = h.vals, X = X, formula = formula,
                       family = "gaussian", link = NULL,
                       k = k, method = "REML", weights = weights)

  function(Xnew) {
    Xnew <- as.data.frame(Xnew)
    as.numeric(stats::predict(proj.fit, newdata = Xnew, type = "response"))
  }
}

# Conditional density estimator using a quantile forest + kernel smoothing.
# Returns list(fhat, h.opt) where fhat(y0, X0) evaluates f(y0 | X0).
.qgam_cond_density <- function(y, X, nfolds = 5,
                                h.mult = seq(0.05, 2, by = 0.05)) {
  if (!is.data.frame(X)) X <- as.data.frame(X)
  y <- as.numeric(y)
  n <- length(y)

  qf <- grf::quantile_forest(X, y)
  W  <- as.matrix(grf::get_forest_weights(qf))

  sigma.y <- stats::sd(y)
  h.star  <- 1.06 * sigma.y * n^(-1/5)
  h.grid  <- h.star * h.mult

  folds <- sample(rep(seq_len(nfolds), length.out = n))

  loglik.scores <- vapply(h.grid, function(h) {
    ll.fold <- vapply(seq_len(nfolds), function(fold) {
      train <- which(folds != fold)
      test  <- which(folds == fold)
      W.tt  <- W[test, train, drop = FALSE]
      Y.tr  <- y[train]
      Y.te  <- y[test]
      Y.mat <- matrix(Y.tr, nrow = length(test), ncol = length(train), byrow = TRUE)
      fhat  <- pmax(rowSums(W.tt * stats::dnorm((Y.te - Y.mat) / h) / h), 1e-12)
      mean(log(fhat))
    }, numeric(1))
    mean(ll.fold)
  }, numeric(1))

  h.opt <- h.grid[which.max(loglik.scores)]

  fhat <- function(y0, X0) {
    if (!is.data.frame(X0)) X0 <- as.data.frame(X0)
    W.new <- as.matrix(grf::get_forest_weights(qf, newdata = X0))
    Y.mat <- matrix(y, nrow = nrow(W.new), ncol = ncol(W.new), byrow = TRUE)
    y.mat <- matrix(y0, nrow = nrow(W.new), ncol = ncol(W.new), byrow = FALSE)
    rowSums(W.new * stats::dnorm((y.mat - Y.mat) / h.opt) / h.opt)
  }

  list(fhat = fhat, h.opt = h.opt)
}
