#' Debiased Score Test for Generalised Additive Models
#'
#' Tests goodness of fit of a GAM by hunting non-parametrically for
#' directions of misspecification. The null model is a GAM fitted via
#' \code{mgcv::gam} (or \code{mgcv::bam}). Internally wires the
#' \code{reg.method}, \code{score.fun}, \code{proj.method}, and
#' \code{hunt.method} arguments of \code{\link{debiased.scoretest}}.
#'
#' @param y Numeric response vector of length n.
#' @param X Numeric covariate matrix or data frame of dimension n x p.
#'   If \code{formula = NULL}, columns are renamed to \code{X1, X2, ...}
#'   internally. If a formula is supplied, column names in \code{X} must
#'   match the variable names used in the formula.
#' @param formula A \code{\link[stats]{formula}} or character string
#'   specifying the null GAM model. If \code{NULL} (default), a full
#'   additive model \code{y ~ s(X1, k=k) + s(X2, k=k) + ...} is used.
#' @param family Character name of the GAM family, e.g. \code{"gaussian"},
#'   \code{"binomial"}, \code{"poisson"}. Default \code{"gaussian"}.
#' @param link Character name of the link function. If \code{NULL}
#'   (default), the canonical link for \code{family} is used.
#' @param k Integer. Basis dimension for smooth terms when
#'   \code{formula = NULL}. Default \code{20}.
#' @param method Smoothing parameter estimation method passed to
#'   \code{mgcv::gam}. Default \code{"REML"}.
#' @param bam Logical. If \code{TRUE}, use \code{mgcv::bam} instead of
#'   \code{mgcv::gam} (recommended for large n). Default \code{FALSE}.
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
#'   \code{\link{dst.gam.nested}}
#'
#' @export
dst.gam <- function(y, X,
                    formula = NULL,
                    family = "gaussian",
                    link = NULL,
                    k = 20,
                    method = "REML",
                    bam = FALSE,
                    hunt.style = "opt.WLS",
                    in.sample = TRUE,
                    return.hunted = FALSE,
                    multi.split = FALSE,
                    n.splits = 5,
                    B = 30,
                    agg.fun = function(x) mean(x, na.rm = TRUE),
                    ...) {

  reg.method <- function(y, X) {
    .gam_fit(y, X, formula = formula, family = family, link = link,
             k = k, method = method, bam = bam)
  }

  score.fun <- .gam_score

  proj.method <- function(fit, h.hat, y, X) {
    .gam_proj(fit, h.hat, y, X, formula = formula, k = k, method = method, bam = bam)
  }

  hunt.method <- function(resids, X, y, fit) {
    hunt.nonpara(resids, X, y, fit,
                 hunt.style = hunt.style,
                 weight.fit = .gam_weight_fit,
                 proj.hunt.method = .gam_proj_hunt(formula, k, method, bam),
                 ...)
  }

  debiased.scoretest(y, X,
                     score.fun    = score.fun,
                     reg.method   = reg.method,
                     hunt.method  = hunt.method,
                     proj.method  = proj.method,
                     in.sample    = in.sample,
                     return.hunted = return.hunted,
                     multi.split  = multi.split,
                     n.splits     = n.splits,
                     B            = B,
                     agg.fun      = agg.fun)
}


#' Debiased Score Test for Nested GAM Comparison
#'
#' Tests whether a larger GAM (\code{formula.alt}) fits significantly better
#' than a smaller null GAM (\code{formula.null}), analogous to
#' \code{anova(gam0, gam1)} but using the Debiased Score Test.
#'
#' If \code{formula.alt} is supplied, the hunting direction is the
#' difference in fitted values between the alternative and null models
#' (parametric hunt). If \code{formula.alt = NULL}, non-parametric GRF
#' hunting is used instead.
#'
#' @param y Numeric response vector of length n.
#' @param X Numeric covariate matrix or data frame of dimension n x p.
#'   Column names must match those used in \code{formula.null} and
#'   \code{formula.alt}.
#' @param formula.null Formula or character string for the null (smaller) GAM.
#' @param formula.alt Formula or character string for the alternative (larger)
#'   GAM. If \code{NULL} (default), non-parametric hunting is used.
#' @param family Character name of the GAM family. Default \code{"gaussian"}.
#' @param link Character link function name. Default \code{NULL} (canonical).
#' @param k Integer basis dimension for smooth terms. Default \code{20}.
#' @param method Smoothing parameter estimation method. Default \code{"REML"}.
#' @param bam Logical. Use \code{mgcv::bam} instead of \code{mgcv::gam}.
#'   Default \code{FALSE}.
#' @param hunt.style Character. One of \code{"opt.WLS"} (default),
#'   \code{"WLS"}, or \code{"vanilla"}. When \code{formula.alt} is supplied,
#'   the same weighting logic applies but the alternative GAM is used as the
#'   learner instead of GRF.
#' @param in.sample Logical. Passed to \code{\link{debiased.scoretest}}.
#'   Default \code{TRUE}.
#' @param return.hunted Logical. Return \code{h.hat} alongside p-value.
#'   Default \code{FALSE}.
#' @param multi.split Logical. Use multi-split aggregation. Default \code{FALSE}.
#' @param n.splits,B,agg.fun Multi-split parameters. See
#'   \code{\link{debiased.scoretest}}.
#' @param ... Additional arguments passed to \code{grf::regression_forest}
#'   when \code{formula.alt = NULL}, or to \code{mgcv::gam}/\code{mgcv::bam}
#'   when \code{formula.alt} is supplied.
#'
#' @return A p-value (numeric scalar), or if \code{return.hunted = TRUE}
#'   and \code{multi.split = FALSE}, a list with \code{p.value} and
#'   \code{h.hat}.
#'
#' @seealso \code{\link{dst.gam}}, \code{\link{debiased.scoretest}}
#'
#' @export
dst.gam.nested <- function(y, X,
                           formula.null,
                           formula.alt = NULL,
                           family = "gaussian",
                           link = NULL,
                           k = 20,
                           method = "REML",
                           bam = FALSE,
                           hunt.style = "opt.WLS",
                           in.sample = TRUE,
                           return.hunted = FALSE,
                           multi.split = FALSE,
                           n.splits = 5,
                           B = 30,
                           agg.fun = function(x) mean(x, na.rm = TRUE),
                           ...) {

  reg.method <- function(y, X) {
    .gam_fit(y, X, formula = formula.null, family = family, link = link,
             k = k, method = method, bam = bam)
  }

  score.fun <- .gam_score

  proj.method <- function(fit, h.hat, y, X) {
    .gam_proj(fit, h.hat, y, X, formula = formula.null, k = k,
              method = method, bam = bam)
  }

  if (!is.null(formula.alt)) {
    # Parametric hunt: use formula.alt as the learner instead of GRF,
    # applying the same WLS / opt.WLS / vanilla weighting logic
    hunt.method <- function(resids, X, y, fit) {
      .gam_hunt(resids, X, y, fit,
                formula.alt  = formula.alt,
                formula.null = formula.null,
                hunt.style   = hunt.style,
                k = k, method = method, bam = bam)
    }
  } else {
    # Non-parametric hunt via GRF
    hunt.method <- function(resids, X, y, fit) {
      hunt.nonpara(resids, X, y, fit,
                   hunt.style = hunt.style,
                   weight.fit = .gam_weight_fit,
                   proj.hunt.method = .gam_proj_hunt(formula.null, k, method, bam),
                   ...)
    }
  }

  debiased.scoretest(y, X,
                     score.fun    = score.fun,
                     reg.method   = reg.method,
                     hunt.method  = hunt.method,
                     proj.method  = proj.method,
                     in.sample    = in.sample,
                     return.hunted = return.hunted,
                     multi.split  = multi.split,
                     n.splits     = n.splits,
                     B            = B,
                     agg.fun      = agg.fun)
}


# ---- internal helpers -------------------------------------------------------

# Fit a GAM or BAM
.gam_fit <- function(y, X, formula = NULL, family = "gaussian", link = NULL,
                     k = 20, method = "REML", weights = NULL, bam = FALSE) {
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

  fam.obj <- .gam_family(family, link)

  if (bam) {
    mgcv::bam(formula, data = data.frame(y, X), family = fam.obj,
              weights = weights, discrete = TRUE)
  } else {
    mgcv::gam(formula, data = data.frame(y, X), family = fam.obj,
              method = method, weights = weights)
  }
}

# Build a family object from name + optional link
.gam_family <- function(family, link) {
  fam.fun <- match.fun(family)
  if (is.null(link)) fam.fun() else fam.fun(link = link)
}

# Score function: response-scale residuals
.gam_score <- function(y, X, fit) {
  as.numeric(y) - as.numeric(stats::predict(fit, newdata = as.data.frame(X),
                                             type = "response"))
}

# proj.method: project h.hat onto the null GAM space with GLM weights
.gam_proj <- function(fit, h.hat, y, X, formula = NULL, k = 20,
                      method = "REML", bam = FALSE) {
  X       <- as.data.frame(X)
  h.vals  <- h.hat(X)

  lnk         <- fit$family$link
  weight.fun  <- stats::make.link(lnk)$mu.eta
  fitted.link <- stats::predict(fit, newdata = X, type = "link")
  weights     <- as.numeric(weight.fun(fitted.link))

  proj.fit <- .gam_fit(y = h.vals, X = X, formula = formula,
                       family = "gaussian", link = NULL,
                       k = k, method = method, weights = weights, bam = bam)

  function(Xnew) {
    as.numeric(stats::predict(proj.fit, newdata = as.data.frame(Xnew),
                              type = "response"))
  }
}

# weight.fit for opt.WLS: returns link derivative at fitted values
.gam_weight_fit <- function(y, X, fit) {
  function(Xnew) {
    lnk    <- fit$family$link
    wfun   <- stats::make.link(lnk)$mu.eta
    as.numeric(wfun(stats::predict(fit, newdata = as.data.frame(Xnew),
                                   type = "link")))
  }
}

# Parametric hunt using formula.alt as the learner (mirrors hunt.nonpara logic)
# Hunting is always gaussian — we are regressing scores, not the original response.
.gam_hunt <- function(resids, X, y, fit,
                      formula.alt, formula.null,
                      hunt.style = "opt.WLS",
                      k = 20, method = "REML", bam = FALSE) {

  stopifnot(hunt.style %in% c("opt.WLS", "WLS", "vanilla"))

  resids <- ifelse(abs(resids) < 1e-9, sign(resids) * 1e-9, resids)
  X      <- as.data.frame(X)

  .gam_pred_fun <- function(gam.fit) {
    function(Xnew) as.numeric(stats::predict(gam.fit, newdata = as.data.frame(Xnew),
                                              type = "response"))
  }

  if (hunt.style == "vanilla") {
    gam.fit <- .gam_fit(resids, X, formula = formula.alt, family = "gaussian",
                        link = NULL, k = k, method = method, bam = bam)
    gam.pred <- .gam_pred_fun(gam.fit)

  } else if (hunt.style == "WLS") {
    gam.fit  <- .gam_fit(1 / resids, X, formula = formula.alt, family = "gaussian",
                         link = NULL, k = k, method = method,
                         weights = resids^2, bam = bam)
    gam.pred <- .gam_pred_fun(gam.fit)

  } else { # opt.WLS
    # Base learner: predict 1/resids weighted by resids^2
    gam.fit.base <- .gam_fit(1 / resids, X, formula = formula.alt,
                             family = "gaussian", link = NULL,
                             k = k, method = method, weights = resids^2, bam = bam)
    gam.pred.base <- .gam_pred_fun(gam.fit.base)

    # Variance function: estimated via a separate fit on resids^2.
    # Clipped to positive to avoid negative weights downstream.
    v.fit     <- .gam_fit(resids^2, X, formula = formula.alt, family = "gaussian",
                          link = NULL, k = k, method = method, bam = bam)
    v.fun.raw <- .gam_pred_fun(v.fit)
    v.fun     <- function(x) pmax(v.fun.raw(x), 1e-8)

    # Weight function from null model link derivative
    weight.fun <- .gam_weight_fit(y, X, fit)

    w.tilde <- (weight.fun(X)^2) / v.fun(X)
    h.tilde <- function(x) gam.pred.base(x) * v.fun(x) / weight.fun(x)

    # Project h.tilde onto null model space
    proj.hunt <- .gam_proj_hunt(formula.null, k, method, bam)
    m.h.hat   <- proj.hunt(h.tilde(X), X = X, weights = w.tilde)

    raw.pred <- function(x) {
      gam.pred.base(x) - m.h.hat(x) * weight.fun(x) / v.fun(x)
    }

    # Reproject back onto formula.alt space — the scaling w/v can push the
    # result outside the model space
    raw.vals  <- raw.pred(X)
    refit     <- .gam_fit(raw.vals, X, formula = formula.alt, family = "gaussian",
                          link = NULL, k = k, method = method, bam = bam)
    gam.pred  <- function(x) {
      as.numeric(stats::predict(refit, newdata = as.data.frame(x), type = "response"))
    }
  }

  # Tukey IQR clipping
  train.vals  <- gam.pred(X)
  q1          <- stats::quantile(train.vals, 0.25)
  q3          <- stats::quantile(train.vals, 0.75)
  iqr         <- q3 - q1
  lower.bound <- q1 - 1.5 * iqr
  upper.bound <- q3 + 1.5 * iqr

  function(x) {
    vals <- gam.pred(x)
    vals[vals < lower.bound] <- lower.bound
    vals[vals > upper.bound] <- upper.bound
    vals
  }
}


# proj.hunt.method factory for use inside hunt.nonpara (opt.WLS)
.gam_proj_hunt <- function(formula, k, method, bam) {
  function(h.vals, X, weights) {
    X <- as.data.frame(X)
    if (is.null(formula)) colnames(X) <- paste0("X", seq_len(ncol(X)))
    pfit <- .gam_fit(y = h.vals, X = X, formula = formula,
                     family = "gaussian", link = NULL,
                     k = k, method = method, weights = weights, bam = bam)
    function(Xnew) {
      Xnew <- as.data.frame(Xnew)
      if (is.null(formula)) colnames(Xnew) <- paste0("X", seq_len(ncol(Xnew)))
      as.numeric(stats::predict(pfit, newdata = Xnew, type = "response"))
    }
  }
}
