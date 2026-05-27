#' Debiased score test: goodness-of-fit test and model comparison
#'
#' Test whether a parametric (e.g., glm) or a semiparametric (e.g., GAM) model 
#' is well-specified. The test is a debiased (Neyman-orthogonalized) score test 
#' computed via sample splitting:
#' on a held-out hunt sample, a flexible auxiliary fit is hunted for a
#' direction in which the null model's score is non-zero; on a held-out test
#' sample, that direction's score is evaluated and standardized under the
#' null. The orthogonalization absorbs plug-in bias from estimating the
#' direction, so the resulting test statistic is asymptotically standard
#' normal under the null without requiring a parametric form for the
#' alternative.
#' 
#' For most scenarios, use one of these methods instead:
#' 
#' \itemize{
#'   \item Use \code{\link{gof_test}} to test whether a fitted model is
#'   well-specified against a nonparametric alternative. S3 methods are
#'   provided for \code{glm} (\code{\link{gof_test.glm}}), \code{lm}
#'   (\code{\link{gof_test.lm}}) and \code{mgcv::gam}
#'   (\code{\link{gof_test.gam}}).
#'
#' \item Use \code{\link{compare_models}} to test a null model \code{fit.0}
#'   against an alternative supermodel \code{fit.1} in the same model class. 
#'   Similar to \code{\link{stats:anova}}, method can be used to conduct a 
#'   significance test of one or more predictors. 
#'   In contrast with \code{\link{gof_test}}, this method targets the alternative 
#'   \code{fit.1}.
#'   S3 methods are provided for \code{glm} (\code{\link{compare_models.glm}}),
#'   \code{lm} (\code{\link{compare_models.lm}}) and \code{mgcv::gam}
#'   (\code{\link{compare_models.gam}}).}
#'
#' Use \code{\link{dScoreTest}} directly for full control over the score,
#'   weight, refit and hunt routines: this is the underlying engine that the
#'   S3 methods wrap. 
#'
#' @docType package
#' @name dScoreTest
"_PACKAGE"
