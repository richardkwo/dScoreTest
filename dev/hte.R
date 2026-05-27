#' Debiased Score Test for Heterogeneous Treatment Effects
#'
#' Tests whether treatment effects are heterogeneous across covariates using
#' the Debiased Score Test. The null model is a partially linear model (PLM)
#' fitted via double machine learning (DML) with GRF. The test hunts
#' non-parametrically for heterogeneity using a T-learner over all covariates.
#'
#' @param y Numeric response vector of length n.
#' @param Z Numeric covariate matrix or data frame of dimension n x p
#'   (not including the treatment variable).
#' @param D Numeric treatment vector of length n. Can be binary or continuous.
#' @param het.vars Integer vector of column indices of \code{Z} to
#'   \emph{test} for heterogeneity. The null model allows heterogeneity in
#'   all columns of \code{Z} \emph{except} those in \code{het.vars}; the
#'   hunt then searches non-parametrically for heterogeneity across all
#'   covariates. \code{NULL} (default) tests all variables, with a constant
#'   treatment effect as the null.
#' @param binary.treatment Logical. If \code{TRUE}, indicates that \code{D}
#'   is binary. Used in \code{opt.WLS} to estimate the variance function via
#'   \code{fit(Z) * (1 - fit(Z))} rather than a separate forest. Default
#'   \code{FALSE}.
#' @param n.folds Integer. Number of cross-fitting folds for DML. Default
#'   \code{2}.
#' @param hunt.style Character. One of \code{"opt.WLS"} (default),
#'   \code{"WLS"}, or \code{"vanilla"}.
#' @param in.sample Logical. Passed to \code{\link{debiased.scoretest}}.
#'   Default \code{TRUE}.
#' @param return.hunted Logical. If \code{TRUE}, returns a list with
#'   \code{p.value} and \code{h.hat}. Default \code{FALSE}.
#' @param multi.split Logical. Passed to \code{\link{debiased.scoretest}}.
#'   Default \code{FALSE}.
#' @param n.splits,B,agg.fun Multi-split parameters. See
#'   \code{\link{debiased.scoretest}}.
#' @param arg.dml Named list of additional arguments passed to
#'   \code{grf::boosted_regression_forest} within the DML fitting step.
#'   Default \code{list()}.
#' @param arg.hunt Named list of additional arguments passed to
#'   \code{grf::regression_forest} within the hunting step.
#'   Default \code{list()}.
#'
#' @return A p-value (numeric scalar), or if \code{return.hunted = TRUE}
#'   and \code{multi.split = FALSE}, a list with \code{p.value} and
#'   \code{h.hat}.
#'
#' @seealso \code{\link{debiased.scoretest}}, \code{\link{hunt.nonpara}}
#'
#' @export
dst.hte <- function(y, Z, D,
                    het.vars = NULL,
                    binary.treatment = FALSE,
                    n.folds = 2,
                    hunt.style = "opt.WLS",
                    in.sample = TRUE,
                    return.hunted = FALSE,
                    multi.split = FALSE,
                    n.splits = 5,
                    B = 30,
                    agg.fun = function(x) mean(x, na.rm = TRUE),
                    arg.dml = list(),
                    arg.hunt = list()) {

  # Pack D and Z into a single matrix X; treatment is always column 1
  X <- .hte_combine(D, Z)

  reg.method <- function(y, X) {
    do.call(.hte_fit, c(list(y = y, X = X, het.vars = het.vars,
                             n.folds = n.folds), arg.dml))
  }

  score.fun <- .hte_score

  proj.method <- function(fit, h.hat, y, X) {
    do.call(.hte_proj, c(list(fit = fit, h.hat = h.hat, y = y, X = X,
                              het.vars = het.vars, n.folds = n.folds), arg.dml))
  }

  hunt.method <- function(resids, X, y, fit) {
    do.call(.hte_hunt, c(list(resids = resids, X = X, y = y, fit = fit,
                              hunt.style = hunt.style,
                              binary.treatment = binary.treatment),
                         arg.hunt))
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

# Pack D (treatment) and Z (covariates) into one matrix, D in column 1
.hte_combine <- function(D, Z) {
  Z <- as.data.frame(Z)
  as.matrix(cbind(D = as.numeric(D), Z))
}

# Unpack the combined matrix X into treatment D and covariates Z
.hte_split <- function(X) {
  X <- as.data.frame(X)
  list(D = as.numeric(X[, 1]),
       Z = X[, -1, drop = FALSE])
}

# Fit a PLM via DML.
# Null allows heterogeneity in the complement of het.vars;
# het.vars are the variables being tested.
.hte_fit <- function(y, X, het.vars = NULL, n.folds = 2, weights = NULL, ...) {
  s   <- .hte_split(X)
  n.z <- ncol(s$Z)
  het.selection <- if (is.null(het.vars)) NULL else {
    comp <- setdiff(seq_len(n.z), het.vars)
    if (length(comp) == 0) NULL else comp
  }
  model <- .plm_grf(y, s$D, s$Z, het.selection = het.selection,
                    n.folds = n.folds, weights = weights, ...)
  list(model = model, het.vars = het.vars)
}

# Score: residuals from the PLM
.hte_score <- function(y, X, fit) {
  s <- .hte_split(X)
  y - fit$model$predict(s$D, s$Z)
}

# proj.method: project h.hat onto the PLM space
.hte_proj <- function(fit, h.hat, y, X, het.vars = NULL, n.folds = 2, weights = NULL, ...) {
  h.vals   <- h.hat(X)
  proj.fit <- .hte_fit(h.vals, X, het.vars = het.vars,
                       n.folds = n.folds, weights = weights, ...)$model
  function(Xnew) {
    sn <- .hte_split(Xnew)
    proj.fit$predict(sn$D, sn$Z)
  }
}

# Hunt for heterogeneous treatment effects using a T-learner over all Z
.hte_hunt <- function(resids, X, y, fit,
                      hunt.style = "opt.WLS",
                      binary.treatment = FALSE,
                      trim.outlier = TRUE,
                      ...) {

  stopifnot(hunt.style %in% c("opt.WLS", "WLS", "vanilla"))

  resids <- ifelse(abs(resids) < 1e-9, ifelse(resids >= 0, 1e-9, -1e-9), resids)

  s <- .hte_split(X)
  D <- s$D
  Z <- s$Z  # hunt over all covariates

  # T-learner: fits separate forests for treated/control, predicts mu0 + D*(mu1-mu0)
  .t_learner <- function(y.train, D.train, Z.train) {
    rf1 <- grf::regression_forest(Z.train[D.train == 1, , drop = FALSE],
                                  y.train[D.train == 1], ...)
    rf0 <- grf::regression_forest(Z.train[D.train == 0, , drop = FALSE],
                                  y.train[D.train == 0], ...)
    function(D.new, Z.new) {
      Z.new <- as.data.frame(Z.new)
      mu1   <- predict(rf1, Z.new)$predictions
      mu0   <- predict(rf0, Z.new)$predictions
      mu0 + as.numeric(D.new) * (mu1 - mu0)
    }
  }

  .t_learner_wls <- function(y.train, D.train, Z.train, w.train) {
    rf1 <- grf::regression_forest(Z.train[D.train == 1, , drop = FALSE],
                                  y.train[D.train == 1],
                                  sample.weights = w.train[D.train == 1], ...)
    rf0 <- grf::regression_forest(Z.train[D.train == 0, , drop = FALSE],
                                  y.train[D.train == 0],
                                  sample.weights = w.train[D.train == 0], ...)
    function(D.new, Z.new) {
      Z.new <- as.data.frame(Z.new)
      mu1   <- predict(rf1, Z.new)$predictions
      mu0   <- predict(rf0, Z.new)$predictions
      mu0 + as.numeric(D.new) * (mu1 - mu0)
    }
  }

  pred_on_X <- function(pred.fun, X) {
    s <- .hte_split(X)
    pred.fun(s$D, s$Z)
  }

  if (hunt.style == "vanilla") {
    pred.fun <- .t_learner(resids, D, Z)
    rf.pred  <- function(X) pred_on_X(pred.fun, X)

  } else if (hunt.style == "WLS") {
    pred.fun <- .t_learner_wls(1 / resids, D, Z, resids^2)
    rf.pred  <- function(X) pred_on_X(pred.fun, X)

  } else { # opt.WLS
    pred.base.fun <- .t_learner_wls(1 / resids, D, Z, resids^2)
    rf.pred.base  <- function(X) pred_on_X(pred.base.fun, X)

    # Variance function
    if (binary.treatment) {
      v.fun <- function(X) {
        s <- .hte_split(X)
        p <- pmax(fit$model$propensity(s$Z), 1e-8)
        pmax(p * (1 - p), 1e-8)
      }
    } else {
      v.fit.fun <- .t_learner(resids^2, D, Z)
      v.fun.raw <- function(X) pred_on_X(v.fit.fun, X)
      v.fun     <- function(X) pmax(v.fun.raw(X), 1e-8)
    }

    # proj.hunt: project onto PLM space with weights w.tilde
    proj.hunt <- function(h.vals, X, weights) {
      proj.fit <- .hte_fit(h.vals, X, weights = weights)$model
      function(Xnew) {
        sn <- .hte_split(Xnew)
        proj.fit$predict(sn$D, sn$Z)
      }
    }

    w.tilde <- 1 / v.fun(X)
    h.tilde <- function(x) rf.pred.base(x) * v.fun(x)
    m.h.hat <- proj.hunt(h.tilde(X), X = X, weights = w.tilde)
    rf.pred  <- function(x) rf.pred.base(x) - m.h.hat(x) / v.fun(x)
  }

  if (!trim.outlier) return(rf.pred)

  # Tukey IQR clipping on training predictions
  train.vals  <- rf.pred(X)
  q1          <- stats::quantile(train.vals, 0.25)
  q3          <- stats::quantile(train.vals, 0.75)
  iqr         <- q3 - q1
  lower.bound <- q1 - 1.5 * iqr
  upper.bound <- q3 + 1.5 * iqr

  function(x) {
    vals <- rf.pred(x)
    vals[vals < lower.bound] <- lower.bound
    vals[vals > upper.bound] <- upper.bound
    vals
  }
}


# ---- DML fitting (PLM via cross-fitted boosted forests) ---------------------

.plm_grf <- function(Y, D, Z,
                     het.selection = NULL,
                     n.folds = 2,
                     weights = NULL,
                     crossfit = TRUE,
                     ...) {

  if (!is.data.frame(Z)) Z <- as.data.frame(Z)
  D <- as.numeric(D)
  Y <- as.numeric(Y)
  n <- length(Y)
  if (is.null(weights)) weights <- rep(1, n)

  Y_tilde  <- numeric(n)
  D_tilde  <- numeric(n)
  f_models <- vector("list", n.folds)

  if (crossfit) {
    folds <- sample(rep(seq_len(n.folds), length.out = n))

    for (k in seq_len(n.folds)) {
      train_idx <- which(folds == k)
      test_idx  <- setdiff(seq_len(n), train_idx)

      g_model <- grf::boosted_regression_forest(
        Z[train_idx, , drop = FALSE], Y[train_idx],
        sample.weights = weights[train_idx], ...)
      m_model <- grf::boosted_regression_forest(
        Z[train_idx, , drop = FALSE], D[train_idx],
        sample.weights = weights[train_idx], ...)

      Y_tilde[test_idx] <- Y[test_idx] -
        predict(g_model, Z[test_idx, , drop = FALSE])$predictions
      D_tilde[test_idx] <- D[test_idx] -
        predict(m_model, Z[test_idx, , drop = FALSE])$predictions

      rm(g_model, m_model); gc()

      if (is.null(het.selection)) {
        f_val <- stats::weighted.mean(Y_tilde[test_idx] * D_tilde[test_idx],
                                      weights[test_idx]) /
          stats::weighted.mean(D_tilde[test_idx]^2, weights[test_idx])
        f_models[[k]] <- list(constant = TRUE, value = f_val)
      } else {
        Z_S       <- Z[test_idx, het.selection, drop = FALSE]
        D_tilde_k <- D_tilde[test_idx]
        D_tilde_k <- ifelse(abs(D_tilde_k) < 1e-6,
                            ifelse(D_tilde_k >= 0, 1e-6, -1e-6), D_tilde_k)
        f_models[[k]] <- grf::boosted_regression_forest(
          X              = Z_S,
          Y              = Y_tilde[test_idx] / D_tilde_k,
          sample.weights = weights[test_idx] * D_tilde_k^2, ...)
        gc()
      }
    }

    f_hat_all <- numeric(n)
    for (k in seq_len(n.folds)) {
      test_idx <- which(folds == k)
      f_hat_all[test_idx] <- if (is.null(het.selection)) {
        f_models[[k]]$value
      } else {
        predict(f_models[[k]],
                Z[test_idx, het.selection, drop = FALSE])$predictions
      }
    }
    g_refit <- grf::boosted_regression_forest(Z, Y - D * f_hat_all,
                                              sample.weights = weights, ...)

  } else {
    g_model <- grf::boosted_regression_forest(Z, Y, sample.weights = weights, ...)
    m_model <- grf::boosted_regression_forest(Z, D, sample.weights = weights, ...)
    g_hat   <- predict(g_model, Z)$predictions
    m_hat   <- predict(m_model, Z)$predictions
    Y_tilde <- Y - g_hat
    D_tilde <- D - m_hat

    if (is.null(het.selection)) {
      f_val <- stats::weighted.mean(Y_tilde * D_tilde, weights) /
        stats::weighted.mean(D_tilde^2, weights)
      f_models[[1]] <- list(constant = TRUE, value = f_val)
    } else {
      Z_S     <- Z[, het.selection, drop = FALSE]
      D_tilde <- ifelse(abs(D_tilde) < 1e-6,
                        ifelse(D_tilde >= 0, 1e-6, -1e-6), D_tilde)
      f_models[[1]] <- grf::boosted_regression_forest(
        X              = Z_S,
        Y              = Y_tilde / D_tilde,
        sample.weights = weights * D_tilde^2, ...)
    }
  }

  # Prediction function for new (D, Z)
  predict_fun <- function(D_new, Z_new) {
    if (!is.data.frame(Z_new)) Z_new <- as.data.frame(Z_new)

    if (is.null(het.selection)) {
      f_hat <- if (crossfit)
        mean(sapply(f_models, function(m) m$value))
      else
        f_models[[1]]$value
      f_hat <- rep(f_hat, nrow(Z_new))
    } else {
      Z_Snew <- Z_new[, het.selection, drop = FALSE]
      preds  <- if (crossfit)
        sapply(f_models, function(m) predict(m, Z_Snew)$predictions)
      else
        predict(f_models[[1]], Z_Snew)$predictions
      f_hat <- if (is.matrix(preds)) rowMeans(preds) else preds
    }

    if (crossfit) {
      g_hat <- predict(g_refit, Z_new)$predictions
      as.numeric(g_hat + D_new * f_hat)
    } else {
      g_hat <- predict(g_model, Z_new)$predictions
      m_hat <- predict(m_model, Z_new)$predictions
      as.numeric(g_hat + (D_new - m_hat) * f_hat)
    }
  }

  # Propensity score: E[D|Z]. For crossfit, refit on all data for prediction.
  if (crossfit) {
    m_refit <- grf::boosted_regression_forest(Z, D, sample.weights = weights, ...)
    propensity_fun <- function(Z_new) {
      if (!is.data.frame(Z_new)) Z_new <- as.data.frame(Z_new)
      predict(m_refit, Z_new)$predictions
    }
  } else {
    propensity_fun <- function(Z_new) {
      if (!is.data.frame(Z_new)) Z_new <- as.data.frame(Z_new)
      predict(m_model, Z_new)$predictions
    }
  }

  list(predict = predict_fun, propensity = propensity_fun)
}
