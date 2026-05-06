#' Debiased Score Test
#'
#' Tests goodness of fit of a null model by hunting for a direction of
#' misspecification using a non-parametric method, then forming a debiased
#' test statistic via sample splitting.
#'
#' @param y Numeric response vector of length n.
#' @param X Numeric covariate matrix of dimension n x p.
#' @param score_fun Function with signature \code{score_fun(fit, y, X)}
#'   returning a vector of scores (i.e., residuals) \eqn{l'(\hat{f}(x_i), y_i)}.
#' @param weight_fun Function with signature \code{weight_fun(fit, X)} 
#'   that computes the weight \eqn{\mathbb{E}[l''(\hat{f}(x_i), y_i) | x_i]} for each 
#'   row \eqn{x_i} of X. 
#' @param fit_method Function with signature \code{fit_method(y, X, ...)}
#'   that returns a fitted null model \eqn{\hat{f} \in \mathcal{F}} by
#'   minimizing the loss \eqn{\sum_i l(f(x_i), y_i)}. For a fitted \code{f},
#'   it must support \code{predict_fun(f, X)} for evaluation.
#' @param wls_method Function with signature \code{wls_method(y, X, w, ...)}
#'   that fits the null model \eqn{\hat{f} \in \mathcal{F}} with weighted least
#'   squares, i.e., minimizing \eqn{\sum_i w_i (f(x_i) - y_i)^2}. For a fitted
#'   \code{f}, it must support \code{predict_fun(f, X)} for evaluation.
#' @param hunt.style Hunting algorithm with the following options.  
#'   \itemize{
#'   \item \code{'optimal'}: optimal hunting (default); 
#'   \item \code{'wls'}: a simpler hunting using weighted least squares, 
#'   which can be less powerful; 
#'   \item \code{'vanilla'}: an even simpler hunting (not recommended).}
#' @param hunt.method Built-in method for hunting. Currently available:
#'   \itemize{
#'   \item \code{'grf'}: regression forest from package \code{grf}.
#'   }
#'   When this is set to any other value, arguments \code{hunt_fun} and 
#'   \code{predict_fun_alt} must be set properly to supply a customized 
#'   hunting method. 
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
#' @param trim.outlier.hunt If \code{TRUE} (default), outliers produced by 
#'   the hunted function will be trimmed using Tukey's IQR rule. 
#' @param splits Numeric vector of length 2 or 3 giving the relative sizes
#'   of the sample splits; rescaled internally to sum to one.
#'   Default is \code{c(0.5, 0.5)}, which splits data into two halves for
#'   hunt and test respectively. Though typically unnecessary in practice,
#'   one can also specify a 3-way split for hunt, debiasing and test respectively.
#' @param arg.fit_method Named list of additional arguments passed to
#'   \code{fit_method} (default to \code{NULL}).
#' @param arg.wls_method Named list of additional arguments passed to
#'   \code{wls_method} (default to \code{NULL}).
#' @param predict_fun Function with signature \code{predict_fun(fit, X)}
#'   returning a numeric vector of predictions from a fitted null model, which 
#'   is produced by \code{fit_method()} and \code{wls_method()}. Default
#'   \code{stats::predict}. When y is binary, it must also support signature 
#'   \code{predict_fun(fit, X, type='response')} for returning probabilities.
#' @param predict_fun_alt Default \code{NULL}. 
#'   When \code{hunt.method} is not set to a built-in method, this is a function
#'   with signature \code{predict_fun_alt(fit, X)} 
#'   returning a numeric vector of predictions from a fitted alternative model 
#'   produced by \code{hunt_fun()}. 
#'
#' @return An object of class \code{"dScoreTest"}.
#'
#' @seealso \code{\link{new_dScoreTest}}
#'
#' @export
dScoreTest <- function(y, X,
                       score_fun, weight_fun,
                       fit_method, wls_method,
                       hunt.style = "optimal",
                       hunt.method = "grf", hunt_fun = NULL,
                       trim.outlier.hunt=TRUE,
                       splits=c(0.5, 0.5),
                       arg.fit_method=NULL, arg.wls_method=NULL, arg.hunt_fun=NULL,
                       predict_fun=stats::predict,
                       predict_fun_alt=NULL) {
    # check input
    stopifnot(
        "score_fun must be a function" = is.function(score_fun),
        "weight_fun must be a function" = is.function(weight_fun),
        "fit_method must be a function" = is.function(fit_method),
        "wls_method must be a function" = is.function(wls_method),
        "number of obs. in X and y do not match." = length(y) == nrow(X),
        "X must be a matrix, array or data frame" = length(dim(X)) == 2,
        "splits must be of length 2 or 3"         = length(splits) %in% c(2, 3),
        "splits must be positive"                 = all(splits > 0)
    )
    hunt.style <- match.arg(hunt.style, c("optimal", "wls", "vanilla"))
    if (hunt.method=="grf") {
        message("Using grf::regression_forest() for hunting.")
        predict_fun_alt <- predict_fun_alt_grf
        if (hunt.style=="vanilla") {
            hunt_fun <- fit_alt_method_grf
            arg.hunt_fun <- arg.fit_alt_method_grf
        } else {
            hunt_fun <- wls_alt_method_grf
            arg.hunt_fun <- arg.wls_alt_method_grf
        }
    } else {
        message("Hunting using the supplied hunt_fun().")
        stopifnot(
            "hunt_fun is not a function" = is.function(hunt_fun),
            "predict_fun_alt is not a function" = is.function(predict_fun_alt)
        )
    }
    # run constructor
    n <- length(y)
    splits <- splits / sum(splits)
    if (length(splits) == 2) {
        idx.hunt <- sort(sample(n, n * splits[1]))
        idx.debias <- setdiff(1:n, idx.hunt)
        idx.test <- idx.debias
    } else {
        idx.hunt <- sort(sample(1:n, n * splits[1]))
        idx.debias <- sort(sample(setdiff(1:n, idx.hunt), n * splits[2]))
        idx.test <- setdiff(1:n, c(idx.hunt, idx.debias))
    }
    if (hunt.style=="optimal") {
        binary.y <- all(y %in% c(0,1))
        if (binary.y) {
            message("y is binary:\nIt is assumed that predict_fun(fit, X, type='response') produces probabilities.")  
        }
        score.test <- new_dScoreTest(y, X,
                         idx.hunt, idx.debias, idx.test,
                         score_fun, weight_fun,
                         fit_method, wls_method,
                         hunt.style = "optimal",
                         hunt.method = hunt.method,
                         fit_alt_method=NULL, wls_alt_method=hunt_fun,
                         binary.y=binary.y,
                         trim.outlier.hunt=trim.outlier.hunt,
                         predict_fun=predict_fun,
                         predict_fun_alt=predict_fun_alt,
                         arg.fit_method=arg.fit_method,
                         arg.wls_method=arg.wls_method,
                         arg.fit_alt_method=NULL, 
                         arg.wls_alt_method=arg.hunt_fun)
    } else if (hunt.style=="wls") {
        score.test <- new_dScoreTest(y, X,
                         idx.hunt, idx.debias, idx.test,
                         score_fun, weight_fun,
                         fit_method, wls_method,
                         hunt.style = "wls",
                         hunt.method = hunt.method,
                         fit_alt_method=NULL, wls_alt_method=hunt_fun,
                         trim.outlier.hunt=trim.outlier.hunt,
                         predict_fun=predict_fun,
                         predict_fun_alt=predict_fun_alt,
                         arg.fit_method=arg.fit_method,
                         arg.wls_method=arg.wls_method,
                         arg.fit_alt_method=NULL, 
                         arg.wls_alt_method=arg.hunt_fun)
    } else if (hunt.style=="vanilla") {
        score.test <- new_dScoreTest(y, X,
                         idx.hunt, idx.debias, idx.test,
                         score_fun, weight_fun,
                         fit_method, wls_method,
                         hunt.style = "vanilla",
                         hunt.method = hunt.method,
                         fit_alt_method=hunt_fun, wls_alt_method=NULL,
                         trim.outlier.hunt=trim.outlier.hunt,
                         predict_fun=predict_fun,
                         predict_fun_alt=predict_fun_alt,
                         arg.fit_method=arg.fit_method,
                         arg.wls_method=arg.wls_method,
                         arg.fit_alt_method=arg.hunt_fun, 
                         arg.wls_alt_method=NULL)
    }
    return(score.test)
}


#' Constructor for the debiased score test
#'
#' Internal worker that builds a \code{dScoreTest} object from a fixed
#' three-way (or two-way) sample split. Called by \code{\link{dScoreTest}}.
#' Splits, hunting algorithm, and predict semantics are taken as fully
#' resolved arguments — no defaults are inferred from the data.
#'
#' @param y Numeric response vector of length n.
#' @param X Numeric covariate matrix of dimension n x p.
#' @param idx.hunt Integer indices into \code{1:n} for the hunting subsample.
#' @param idx.debias Integer indices into \code{1:n} for the debiasing
#'   subsample (refitting the null model and projecting the hunted direction).
#' @param idx.test Integer indices into \code{1:n} for the test subsample
#'   on which the test statistic is evaluated. May coincide with
#'   \code{idx.debias} (two-way split) or be disjoint (three-way split).
#' @param score_fun Function with signature \code{score_fun(fit, y, X)}
#'   returning a vector of scores \eqn{l'(\hat{f}(x_i), y_i)}.
#' @param weight_fun Function with signature \code{weight_fun(fit, X)}
#'   returning the weight \eqn{\mathbb{E}[l''(\hat{f}(x_i), y_i) | x_i]} for
#'   each row of \code{X}.
#' @param fit_method Function with signature \code{fit_method(y, X, ...)}
#'   returning a fitted null model. The returned object must support
#'   \code{predict_fun(fit, X)}.
#' @param wls_method Function with signature \code{wls_method(y, X, w, ...)}
#'   that fits the null model by weighted least squares. The returned object
#'   must support \code{predict_fun(fit, X)}.
#' @param hunt.style One of \code{"optimal"} (default), \code{"wls"}, or
#'   \code{"vanilla"}. Selects which \code{hunt_*} routine is used.
#' @param hunt.method String for the hunting method.
#' @param fit_alt_method Required when \code{hunt.style = "vanilla"};
#'   ignored otherwise. Function with signature
#'   \code{fit_alt_method(y, X, ...)} returning a fitted alternative model
#'   that supports \code{predict_fun_alt(g, X)}.
#' @param wls_alt_method Required when \code{hunt.style \%in\% c("optimal", "wls")};
#'   ignored otherwise. Function with signature
#'   \code{wls_alt_method(y, X, w, ...)} returning a fitted alternative model
#'   that supports \code{predict_fun_alt(g, X)}.
#' @param X.cols.hunt Integer or name vector selecting which columns of
#'   \code{X} drive the hunt. Default: all columns.
#' @param binary.y Logical. When \code{TRUE}, the optimal hunter computes
#'   \eqn{\mathrm{Var}(l' | x)} from \code{predict_fun(fit, X, type='response')}
#'   assuming a Bernoulli response. Only consulted by
#'   \code{hunt.style = "optimal"}.
#' @param trim.outlier.hunt Logical. Passed to the chosen \code{hunt_*}
#'   routine as \code{trim.outlier}. Default \code{TRUE}.
#' @param predict_fun Function with signature \code{predict_fun(fit, X)}
#'   returning predictions from null-model fits. Default \code{stats::predict}.
#' @param predict_fun_alt Function with signature \code{predict_fun_alt(fit, X)}
#'   returning predictions from alt-model fits. Default \code{stats::predict}.
#' @param arg.fit_method,arg.wls_method,arg.fit_alt_method,arg.wls_alt_method
#'   Named lists of additional arguments forwarded to the corresponding
#'   fitter via \code{do.call}. Default \code{NULL}.
#'
#' @return A list of class \code{"dScoreTest"} with elements:
#'   \describe{
#'     \item{\code{t.stat}}{Debiased test statistic
#'       \eqn{\sqrt{n_{\mathrm{test}}}\,\bar{L}/\hat{\sigma}_L}.}
#'     \item{\code{p.val}}{One-sided p-value (right tail of the standard
#'       normal).}
#'     \item{\code{resids}}{Score residuals on the test subsample.}
#'     \item{\code{h}}{Orthogonalised hunted direction on the test subsample.}
#'     \item{\code{h.raw}}{Hunted direction before the outer debias projection.}
#'     \item{\code{Data}}{List with \code{X}, \code{y}, and the three index
#'       vectors.}
#'     \item{\code{Call}}{Named list of methods,
#'       \code{hunt.style}, \code{hunt.method}, both predict functions, 
#'       and the four \code{arg.*} lists.}
#'   }
#'
#' @seealso \code{\link{dScoreTest}}
#'
#' @export
new_dScoreTest <- function(y, X,
                           idx.hunt, idx.debias, idx.test,
                           score_fun, weight_fun,
                           fit_method, wls_method,
                           hunt.style = "optimal",
                           hunt.method = "customized",
                           fit_alt_method=NULL, wls_alt_method=NULL,
                           X.cols.hunt=1:ncol(X), binary.y=FALSE,
                           trim.outlier.hunt = TRUE,
                           predict_fun = stats::predict, 
                           predict_fun_alt = stats::predict,
                           arg.fit_method = NULL, arg.wls_method=NULL,
                           arg.fit_alt_method = NULL, arg.wls_alt_method=NULL) {
    
    # (fit and) hunt with data in idx.hunt
    fit.hunt <- do.call(fit_method, 
                        c(list(y[idx.hunt], X[idx.hunt,,drop=FALSE]), 
                          arg.fit_method))
    stopifnot(hunt.style %in% c("optimal", "wls", "vanilla"))
    if (hunt.style == "optimal") {
        stopifnot(is.function(wls_alt_method))
        h_hat <- hunt_optimal(wls_alt_method, wls_method, score_fun, weight_fun,
            fit.hunt, y[idx.hunt], X[idx.hunt,,drop=FALSE],
            X.cols = X.cols.hunt,
            binary.y = binary.y,
            trim.outlier = trim.outlier.hunt,
            arg.wls_alt_method = arg.wls_alt_method,
            arg.wls_method = arg.wls_method,
            predict_fun = predict_fun,
            predict_fun_alt = predict_fun_alt)
    } else if (hunt.style == "wls") {
        stopifnot(is.function(wls_alt_method))
        resids.hunt <- score_fun(fit.hunt, y[idx.hunt], X[idx.hunt,,drop=FALSE])
        h_hat <- hunt_wls(wls_alt_method, resids.hunt, X[idx.hunt,,drop=FALSE],
                          X.cols = X.cols.hunt,
                          trim.outlier = trim.outlier.hunt,
                          arg.wls_alt_method = arg.wls_alt_method,
                          predict_fun_alt = predict_fun_alt)
    } else {
        stopifnot(is.function(fit_alt_method))
        resids.hunt <- score_fun(fit.hunt, y[idx.hunt], X[idx.hunt,,drop=FALSE])
        h_hat <- hunt_vanilla(fit_alt_method, resids.hunt, X[idx.hunt,,drop=FALSE],
                          X.cols = X.cols.hunt,
                          trim.outlier = trim.outlier.hunt,
                          arg.fit_alt_method = arg.fit_alt_method,
                          predict_fun_alt = predict_fun_alt)
    }
    # (refit and) debias h_hat
    fit.debias <- do.call(fit_method, 
                          c(list(y[idx.debias], X[idx.debias,,drop=FALSE]), 
                            arg.fit_method))
    w.debias <- weight_fun(fit.debias, X[idx.debias,,drop=FALSE])
    m.h.fit <- do.call(wls_method,
                          c(list(h_hat(X[idx.debias,,drop=FALSE]),
                                 X[idx.debias,,drop=FALSE],
                                 w.debias), arg.wls_method))
    # evaluate the test 
    h.test.raw <- h_hat(X[idx.test,,drop=FALSE])
    h.test <-  h.test.raw - predict_fun(m.h.fit, X[idx.test,,drop=FALSE])
    resids.test <- score_fun(fit.debias, y[idx.test], X[idx.test,,drop=FALSE])
    
    L.test <- resids.test * h.test
    L.test.raw <- resids.test * h.test.raw
    
    t.stat <- sum(L.test) / sqrt(length(L.test) * stats::var(L.test))
    p.val <- stats::pnorm(t.stat, lower.tail = FALSE)
    out <- list(t.stat=t.stat, p.val=p.val,
                resids=resids.test, h=h.test, h.raw=h.test.raw)
    out$Data <- list(X=X, y=y,
                     idx.hunt=idx.hunt, idx.debias=idx.debias, idx.test=idx.test)
    out$Call <- list(score_fun = score_fun, weight_fun = weight_fun,
                        fit_method = fit_method, wls_method = wls_method,
                        hunt.style = hunt.style,
                        hunt.method = hunt.method, 
                        fit_alt_method = fit_alt_method,
                        wls_alt_method = wls_alt_method,
                        predict_fun = predict_fun,
                        predict_fun_alt = predict_fun_alt,
                        arg.fit_method = arg.fit_method,
                        arg.wls_method = arg.wls_method,
                        arg.fit_alt_method = arg.fit_alt_method,
                        arg.wls_alt_method = arg.wls_alt_method)
    class(out) <- "dScoreTest"
    return(out)
}

#' Print the score test
#' @export
print.dScoreTest <- function(x, ...) {
    cat("Debiased score test\n")
    cat(sprintf("(hunt.style = %s, hunt.method = %s)\n\n", 
                x$Call$hunt.style, x$Call$hunt.method))
    if (setequal(x$Data$idx.debias, x$Data$idx.test)) {
        cat(sprintf("n = %d, two-way split: hunt = %d, debias & test = %d\n",
                    length(x$Data$y),
                    length(x$Data$idx.hunt),
                    length(x$Data$idx.test)))
    } else {
        cat(sprintf("n = %d, three-way split: hunt = %d, debias = %d, test = %d\n",
                    length(x$Data$y),
                    length(x$Data$idx.hunt),
                    length(x$Data$idx.debias),
                    length(x$Data$idx.test)))
    }
    cat(sprintf("T = %.4f, p-value = %g\n", x$t.stat, x$p.val))
    invisible(x)
}

#' Plot the score test
#' @export
plot.dScoreTest <- function(x, ...) {
    with(x, {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par))
        par(mfrow = c(1, 2),
            mar = c(4, 4, 1, 0.5),
            oma = c(0, 0, 0, 0))
        # left plot
        plot(resids, h, pch=20, col="blue", cex=0.5, 
             xlab="resid", ylab="hunt", ylim=range(c(h, h.raw)))
        points(resids, h.raw, pch=20, col="grey", cex=0.5, 
             xlab="resid", ylab="hunt")
        .idx.up <- which(h > h.raw)
        .idx.down <- which(h < h.raw)
        if (length(.idx.up) > 0) {
            segments(x0=resids[.idx.up], y0=h[.idx.up], y1=h.raw[.idx.up], 
                     col="red4", lwd=0.5)
        }
        if (length(.idx.down) > 0) {
            segments(x0=resids[.idx.down], y0=h[.idx.down], y1=h.raw[.idx.down], 
                     col="green4", lwd=0.5)
        }
        # right plot
        L.debiased <- resids * h
        L.raw <- resids * h.raw
        .ord <- order(L.debiased)
        L.debiased <- L.debiased[.ord]
        L.raw <- L.raw[.ord]
        .idx.up <- which(L.debiased / sd(L.debiased) > L.raw / sd(L.raw))
        .idx.down <- which(L.debiased / sd(L.debiased) < L.raw / sd(L.raw))
        plot(L.debiased / sd(L.debiased), pch=18, col="blue", cex=0.7,
             xlab = "index (ordered)", ylab="L / sd(L)", 
             ylim=range(c(L.debiased / sd(L.debiased), L.raw / sd(L.raw))))
        points(L.raw / sd(L.raw), col="grey", pch=18, cex=0.7)
        if (length(.idx.up) > 0) {
            segments(x0=.idx.up, 
                     y0=L.debiased[.idx.up] / sd(L.debiased), 
                     y1=L.raw[.idx.up] / sd(L.raw), 
                     col="red4", lwd=0.5)
        }
        if (length(.idx.down) > 0) {
            segments(x0=.idx.down, 
                     y0=L.debiased[.idx.down] / sd(L.debiased), 
                     y1=L.raw[.idx.down] / sd(L.raw), 
                     col="green4", lwd=0.5)
        }
        legend("bottomright", c("debiased", "raw"), pch=c(18,18), 
               col=c("blue", "grey"))
    })
}

