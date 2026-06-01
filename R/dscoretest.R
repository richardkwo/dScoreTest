#' Debiased Score Test
#'
#' @param y Numeric response vector of length n.
#' @param X Numeric covariate matrix of dimension n x p.
#' @param score_fun Function with signature \code{score_fun(fit, y, X)}
#'   returning a vector of scores \eqn{l'(\hat{f}(x_i), y_i)}, which can be viewed
#'   as negative residuals.
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
#' @param trim.outlier.hunt If \code{TRUE} (default), 
#' extreme values produced by the hunted function will be trimmed using Tukey's 
#' IQR rule. 
#' @param X.cols.hunt Integer vector selecting which columns of
#'   \code{X} drive the hunt. Default \code{1:ncol(X)}. This is modified only 
#'   in special settings, e.g., when there is an offset in the null model. 
#' @param splits Numeric vector of length 2 or 3 giving the relative sizes
#'   of the sample splits; rescaled internally to sum to one.
#'   Default is \code{c(0.5, 0.5)}, which splits data into two halves for
#'   hunt and test respectively. Though typically unnecessary in practice,
#'   one can also specify a 3-way split for hunt, debiasing and test respectively.
#' @param arg.fit_method Named list of additional arguments passed to
#'   \code{fit_method} (default to \code{NULL}).
#' @param arg.wls_method Named list of additional arguments passed to
#'   \code{wls_method} (default to \code{NULL}).
#' @param arg.hunt_fun Extra arguments (default \code{NULL}) passed to the 
#'   customized \code{hunt.fun}.
#' @param predict_fun Function with signature \code{predict_fun(fit, X)}
#'   returning a numeric vector of predictions from a fitted null model, which 
#'   is produced by \code{fit_method()} and \code{wls_method()}. 
#'   Note that if \code{fit} is \eqn{\hat{f}}, this function should return 
#'   \eqn{\hat{f}(X)}. Default \code{stats::predict}. 
#'   When y is binary, it must also support signature 
#'   \code{predict_fun(fit, X, type='response')} for returning probabilities.
#' @param predict_fun_alt Default \code{NULL}. 
#'   When \code{hunt.method} is not set to a built-in method, this is a function
#'   with signature \code{predict_fun_alt(fit, X)} 
#'   returning a numeric vector of predictions from a fitted alternative model 
#'   produced by \code{hunt_fun()}. 
#' @param verbose Default \code{FALSE}; information is printed if set to 
#'   \code{TRUE}.
#'
#' @return An object of class \code{"dScoreTest"}: a list whose key elements
#'   are the debiased test statistic \code{t.stat} and the one-sided p-value
#'   \code{p.val} (right tail of the standard normal), along with the test-set
#'   score residuals, the hunted direction, and the call. It has
#'   \code{\link[=print.dScoreTest]{print}},
#'   \code{\link[=summary.dScoreTest]{summary}} and
#'   \code{\link[=plot.dScoreTest]{plot}} methods.
#'
#' @seealso \code{\link{plot.dScoreTest}}, \code{\link{summary.dScoreTest}},
#'   \code{\link{hunt_optimal}}, \code{\link{hunt_wls}}, 
#'   \code{\link{hunt_vanilla}}, \code{\link{new_dScoreTest}}
#' @export
dScoreTest <- function(y, X,
                       score_fun, weight_fun,
                       fit_method, wls_method,
                       hunt.style = "optimal",
                       hunt.method = "grf", hunt_fun = NULL,
                       trim.outlier.hunt=TRUE,
                       X.cols.hunt=1:ncol(X),
                       splits=c(0.5, 0.5),
                       arg.fit_method=NULL, arg.wls_method=NULL, arg.hunt_fun=NULL,
                       predict_fun=stats::predict,
                       predict_fun_alt=NULL, 
                       verbose=FALSE) {
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
        if (verbose) {
            message("Using grf::regression_forest() for hunting.")
        }
        predict_fun_alt <- predict_fun_alt_grf
        if (hunt.style=="vanilla") {
            hunt_fun <- fit_alt_method_grf
            arg.hunt_fun <- arg.fit_alt_method_grf
        } else {
            hunt_fun <- wls_alt_method_grf
            arg.hunt_fun <- arg.wls_alt_method_grf
        }
    } else {
        if (verbose) {
            message("Hunting using the supplied hunt_fun().")
        }
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
        if (binary.y && verbose) {
            message("y is binary:\nIt is assumed that predict_fun(fit, X, type='response') produces probabilities.")  
        }
        score.test <- new_dScoreTest(y, X,
                         idx.hunt, idx.debias, idx.test,
                         score_fun, weight_fun,
                         fit_method, wls_method,
                         hunt.style = "optimal",
                         hunt.method = hunt.method,
                         fit_alt_method=NULL, wls_alt_method=hunt_fun,
                         X.cols.hunt=X.cols.hunt, binary.y=binary.y,
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
                         X.cols.hunt=X.cols.hunt, 
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
                         X.cols.hunt=X.cols.hunt, 
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
#'   routine as \code{trim.outlier}. If \code{TRUE} (default), extreme values 
#'   in the hunted function will be removed using Tukey's IQR rule. 
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
#'     \item{\code{hunted_fun}}{The hunted function that can be applied to X.}      
#'     \item{\code{Data}}{List with \code{X}, \code{y}, and the three index
#'       vectors.}
#'     \item{\code{Call}}{Named list of methods,
#'       \code{hunt.style}, \code{hunt.method}, both predict functions, 
#'       and the four \code{arg.*} lists.}
#'   }
#'
#' @seealso \code{\link{dScoreTest}}, \code{\link{hunt_optimal}}, 
#'   \code{\link{hunt_wls}}, \code{\link{hunt_vanilla}}
#'
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
        # resids is (-score)
        resids.hunt <- -1 * score_fun(fit.hunt, y[idx.hunt], X[idx.hunt,,drop=FALSE])
        h_hat <- hunt_wls(wls_alt_method, resids.hunt, X[idx.hunt,,drop=FALSE],
                          X.cols = X.cols.hunt,
                          trim.outlier = trim.outlier.hunt,
                          arg.wls_alt_method = arg.wls_alt_method,
                          predict_fun_alt = predict_fun_alt)
    } else if (hunt.style == "vanilla") {
        stopifnot(is.function(fit_alt_method))
        # resids is (-score)
        resids.hunt <- -1 * score_fun(fit.hunt, y[idx.hunt], X[idx.hunt,,drop=FALSE])
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
    resids.test <- -1 * score_fun(fit.debias, y[idx.test], X[idx.test,,drop=FALSE])
    
    L.test <- resids.test * h.test
    L.test.raw <- resids.test * h.test.raw
    
    if (any(is.na(L.test))) {
        warning("L contains NAs. Ignoring them.")
    }
    if (stats::var(L.test, na.rm=TRUE)==0) {
        warning("L is constant. Either hunt or debiasing misbehaved.")
    }
    t.stat <- sum(L.test, na.rm=TRUE) / 
        sqrt(sum(!is.na(L.test)) * stats::var(L.test, na.rm=TRUE))
    
    p.val <- stats::pnorm(t.stat, lower.tail = FALSE)
    out <- list(t.stat=t.stat, p.val=p.val,
                resids=resids.test, h=h.test, h.raw=h.test.raw, 
                L=L.test, L.raw=L.test.raw, 
                hunted_fun=h_hat)
    out$Data <- list(X=X, y=y,
                     idx.hunt=idx.hunt, idx.debias=idx.debias, idx.test=idx.test)
    out$Call <- list(score_fun = score_fun, weight_fun = weight_fun,
                        fit_method = fit_method, wls_method = wls_method,
                        hunt.style = hunt.style,
                        hunt.method = hunt.method, 
                        fit_alt_method = fit_alt_method,
                        wls_alt_method = wls_alt_method,
                        X.cols.hunt = X.cols.hunt,
                        predict_fun = predict_fun,
                        predict_fun_alt = predict_fun_alt,
                        arg.fit_method = arg.fit_method,
                        arg.wls_method = arg.wls_method,
                        arg.fit_alt_method = arg.fit_alt_method,
                        arg.wls_alt_method = arg.wls_alt_method)
    class(out) <- "dScoreTest"
    return(out)
}

# generics for the class -------
#' Print the score test
#'
#' @param x A \code{dScoreTest} object.
#' @param ... Unused, for S3 consistency.
#'
#' @return The input \code{x}, invisibly. Called for the side effect of
#'   printing a summary of the test to the console.
#'
#' @export
print.dScoreTest <- function(x, ...) {
    cat("Debiased score test: \n")
    cat(sprintf("y ~ X, with X consists of %s.\n", 
                paste(colnames(x$Data$X)[x$Call$X.cols.hunt], collapse=", ")))
    cat(sprintf("(hunt.style = %s, hunt.method = %s)\n", 
                x$Call$hunt.style, x$Call$hunt.method))
    if (setequal(x$Data$idx.debias, x$Data$idx.test)) {
        cat(sprintf("n = %d, two-way split: hunt = %d, debias & test = %d\n\n",
                    length(x$Data$y),
                    length(x$Data$idx.hunt),
                    length(x$Data$idx.test)))
    } else {
        cat(sprintf("n = %d, three-way split: hunt = %d, debias = %d, test = %d\n\n",
                    length(x$Data$y),
                    length(x$Data$idx.hunt),
                    length(x$Data$idx.debias),
                    length(x$Data$idx.test)))
    }
    cat(sprintf("T = %.4f, p-value = %g\n", x$t.stat, x$p.val))
    invisible(x)
}

#' Plot the score test
#'
#' Diagonotic plots for the test:
#' \enumerate{
#'   \item Histogram of \eqn{\{L_i\}}, where \eqn{L_i = resid_i \times h_i}. 
#'   
#'   \item \eqn{\{L_i\}} against the index \eqn{i}, where \eqn{i} refers to the 
#'     \eqn{i}-th observation in the full dataset. Only those \eqn{i}'s in the 
#'     test split are drawn. The mean is drawn as a horizontal line. 
#'     Extremes values under the null can result in bad normal approximation.
#'     In this case, consider setting \code{trim.outlier.hunt=TRUE}.
#'     
#'   \item Residuals (negative scores) versus the hunted signal. 
#'     A horizontal segment is drawn between each pair of raw hunted
#'     signal and the debiased hunted signal. If debiased gets higher, colored
#'     in red; otherwise colored in green. A regression line (blue) with a large 
#'     positive slope indicates the model is misspecified. 
#'     
#'   \item Normalized \eqn{\{L_i\}} drawn in order.
#' }
#'
#' @param x A \code{dScoreTest} object.
#' @param ... Further graphical parameters passed to underlying plotting
#'   functions.
#'
#' @return No return value; called for its side effect of producing the
#'   diagnostic plots described above.
#'
#' @importFrom graphics par hist abline plot points segments legend
#'
#' @export
plot.dScoreTest <- function(x, ...) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(1, 2),
        mar = c(4, 4, 1, 0.5),
        oma = c(0, 0, 0, 0))
    with(x, {
        # normalize
        L.norm <- L / sd(L)
        L.raw.norm <- L.raw / sd(L.raw)
        # 1st plot
        hist(L.norm, breaks=20, xlab="L / sd(L)", 
             main=sprintf("mean = %.2f", mean(L.norm)))
        abline(v=mean(L.norm), col="red", lwd=1.5)
        # 2nd plot
        plot(Data$idx.test, L, pch=20, cex=0.6, type="p", xlab="index", ylab="L")
        abline(h=0, lty=2)
        abline(h=mean(L), col="red", lwd=1.5)
        # 3rd plot
        plot(h, resids, pch=20, col="blue", cex=0.5, 
             xlab="hunt", ylab="resids", xlim=range(c(h, h.raw)))
        points(h.raw, resids, pch=20, col="grey", cex=0.5)
        abline(h=0, lty=2)
        .idx.up <- which(h > h.raw)
        .idx.down <- which(h < h.raw)
        if (length(.idx.up) > 0) {
            segments(x0=h[.idx.up], x1=h.raw[.idx.up], y0=resids[.idx.up], 
                     col="red4", lwd=0.5)
        }
        if (length(.idx.down) > 0) {
            segments(x0=h[.idx.down], x1=h.raw[.idx.down], y0=resids[.idx.down], 
                     col="green4", lwd=0.5)
        }
        abline(a=0, b=cov(h, resids) / var(h), col="blue", lwd=2)
        # 4th plot
        .ord <- order(L.norm)
        L.norm <- L.norm[.ord]
        L.raw.norm <- L.raw.norm[.ord]
        .idx.up <- which(L.norm > L.raw.norm)
        .idx.down <- which(L.norm < L.raw.norm)
        plot(L.norm, pch=18, col="blue", cex=0.7,
             xlab = "index (ordered)", ylab="L / sd(L)", 
             ylim=range(c(L.norm, L.raw.norm)))
        points(L.raw.norm, col="grey", pch=18, cex=0.7)
        abline(h=0, lty=2)
        if (length(.idx.up) > 0) {
            segments(x0=.idx.up, 
                     y0=L.norm[.idx.up], 
                     y1=L.raw.norm[.idx.up], 
                     col="red4", lwd=0.5)
        }
        if (length(.idx.down) > 0) {
            segments(x0=.idx.down, 
                     y0=L.norm[.idx.down], 
                     y1=L.raw.norm[.idx.down], 
                     col="green4", lwd=0.5)
        }
        legend("bottomright", c("debiased", "raw"), pch=c(18,18),
               col=c("blue", "grey"))
    })
}

#' Summary of the score test
#'
#' Reports the headline statistic and p-value, the sample-split sizes, a
#' raw-vs-debiased comparison of the test statistic (so the effect of the
#' outer-projection debiasing step is visible), and \code{\link[base]{summary}}
#' digests of the diagnostic vectors \code{L}, \code{L.raw}, \code{h},
#' \code{h.raw} and \code{resids}.
#'
#' @param object A \code{dScoreTest} object.
#' @param ... Unused, for S3 consistency.
#'
#' @return A list of class \code{"summary.dScoreTest"}.
#'
#' @export
summary.dScoreTest <- function(object, ...) {
    twoway <- setequal(object$Data$idx.debias, object$Data$idx.test)
    # t-stat that would have resulted without the outer debias projection
    L.raw    <- object$L.raw
    L.raw.ok <- !is.na(L.raw)
    t.raw    <- if (stats::var(L.raw, na.rm = TRUE) == 0) {
        NA_real_
    } else {
        sum(L.raw, na.rm = TRUE) /
            sqrt(sum(L.raw.ok) * stats::var(L.raw, na.rm = TRUE))
    }
    out <- list(
        t.stat        = object$t.stat,
        p.val         = object$p.val,
        t.raw         = t.raw,
        p.raw         = stats::pnorm(t.raw, lower.tail = FALSE),
        n             = length(object$Data$y),
        n.hunt        = length(object$Data$idx.hunt),
        n.debias      = length(object$Data$idx.debias),
        n.test        = length(object$Data$idx.test),
        twoway        = twoway,
        hunt.style    = object$Call$hunt.style,
        hunt.method   = object$Call$hunt.method,
        na.L          = sum(is.na(object$L)),
        L.summary     = summary(object$L),
        L.raw.summary = summary(object$L.raw),
        h.summary     = summary(object$h),
        h.raw.summary = summary(object$h.raw),
        resids.summary = summary(object$resids)
    )
    class(out) <- "summary.dScoreTest"
    out
}

#' @param x A \code{summary.dScoreTest} object.
#' @rdname summary.dScoreTest
print.summary.dScoreTest <- function(x, ...) {
    cat("Debiased score test\n")
    cat(sprintf("(hunt.style = %s, hunt.method = %s)\n",
                x$hunt.style, x$hunt.method))
    if (x$twoway) {
        cat(sprintf("n = %d, two-way split: hunt = %d, debias & test = %d\n",
                    x$n, x$n.hunt, x$n.test))
    } else {
        cat(sprintf("n = %d, three-way split: hunt = %d, debias = %d, test = %d\n",
                    x$n, x$n.hunt, x$n.debias, x$n.test))
    }
    cat(sprintf("\n  Debiased:  T = %8.4f,  p = %g\n", x$t.stat, x$p.val))
    cat(sprintf("  Raw:       T = %8.4f,  p = %g  (may contain bias)\n",
                x$t.raw, x$p.raw))
    if (x$na.L > 0) {
        cat(sprintf("\nNAs in L: %d (excluded from T)\n", x$na.L))
    }
    cat("\nL = resids * h (debiased):\n");        print(x$L.summary)
    cat("\nL.raw = resids * h.raw:\n");           print(x$L.raw.summary)
    cat("\nh (debiased hunted direction):\n");    print(x$h.summary)
    cat("\nh.raw (before outer debias):\n");      print(x$h.raw.summary)
    cat("\nresids (score residuals on test):\n"); print(x$resids.summary)
    invisible(x)
}
