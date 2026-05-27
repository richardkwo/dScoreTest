# Temporary test script for dst.gam and dst.gam.nested
# Run from the dst package root: source("test_gam.R")

# ---- Packages ---------------------------------------------------------------
library(mgcv)
library(grf)
library(MASS)

# ---- Source package files ---------------------------------------------------
source("R/dst.R")
source("R/hunt.R")
source("R/gam.R")

# ---- Data generation --------------------------------------------------------
# Gaussian: null model is y ~ sin(2*X1), alternative adds X1*X3
# X is AR(1) correlated

gen_data <- function(n = 1000, p = 3, tau = 0, family = "gaussian") {
  Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
  X     <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)

  eta <- sin(X[, 1]) + (tau / sqrt(n)) * X[, 1] * X[, 3]

  if (family == "gaussian") {
    y <- eta + rnorm(n)
  } else if (family == "binomial") {
    y <- rbinom(n, 1, plogis(eta))
  }

  list(y = y, X = X)
}

# ---- Helpers ----------------------------------------------------------------
pval_label <- function(p) {
  if (p < 0.01) "** (p < 0.01)"
  else if (p < 0.05) "* (p < 0.05)"
  else sprintf("(p = %.3f, not significant)", p)
}

cat_result <- function(label, p) {
  cat(sprintf("%-55s p = %.4f  %s\n", label, p, pval_label(p)))
}

# =============================================================================
cat("\n====  dst.gam  (gaussian)  ====\n\n")

d0 <- gen_data(n = 1000, p = 5, tau = 0)
d5 <- gen_data(n = 1000, p = 5, tau = 5)


p1 <- dst.gam(d0$y, d0$X, family = "gaussian")
cat_result("null (tau=0)", p1)

p2 <- dst.gam(d5$y, d5$X, family = "gaussian")
cat_result("alternative (tau=5)", p2)

p3 <- dst.gam(d5$y, d5$X, family = "gaussian", hunt.style = "WLS")
cat_result("alternative (tau=5, WLS)", p3)

p4 <- dst.gam(d5$y, d5$X, family = "gaussian", hunt.style = "vanilla")
cat_result("alternative (tau=5, vanilla)", p4)

# =============================================================================
cat("\n====  dst.gam  (binomial)  ====\n\n")

db0 <- gen_data(n = 1500, p = 5, tau = 0,  family = "binomial")
db5 <- gen_data(n = 1500, p = 5, tau = 10, family = "binomial")


p5 <- dst.gam(db0$y, db0$X, family = "binomial")
cat_result("null (tau=0)", p5)

p6 <- dst.gam(db5$y, db5$X, family = "binomial")
cat_result("alternative (tau=10)", p6)

# =============================================================================
cat("\n====  dst.gam  (return.hunted)  ====\n\n")

res <- dst.gam(d5$y, d5$X, family = "gaussian", return.hunted = TRUE)
cat(sprintf("p-value: %.4f\n", res$p.value))
cat(sprintf("h.hat is a function: %s\n", is.function(res$h.hat)))
cat(sprintf("h.hat predictions on first 3 obs: %s\n",
            paste(round(res$h.hat(d5$X[1:3, ]), 3), collapse = ", ")))


# =============================================================================
cat("\n====  dst.gam.nested  (no formula.alt — non-parametric hunt)  ====\n\n")

formula.null <- "y ~ s(X1, k=10) + s(X2, k=10) + s(X3, k=10) + s(X4, k=10) + s(X5, k=10)"

p9  <- dst.gam.nested(d0$y, d0$X, formula.null = formula.null)
cat_result("null (tau=0)", p9)

p10 <- dst.gam.nested(d5$y, d5$X, formula.null = formula.null)
cat_result("alternative (tau=5)", p10)

# =============================================================================
cat("\n====  dst.gam.nested  (formula.alt supplied — parametric hunt)  ====\n\n")

formula.alt <- "y ~ s(X1, k=10) + s(X2, k=10) + s(X3, k=10) + s(X4, k=10) + s(X5, k=10) + ti(X1, X3)"

p11 <- dst.gam.nested(d0$y, d0$X,
                      formula.null = formula.null,
                      formula.alt  = formula.alt)
cat_result("null (tau=0)", p11)

p12 <- dst.gam.nested(d5$y, d5$X,
                      formula.null = formula.null,
                      formula.alt  = formula.alt)
cat_result("alternative (tau=5)", p12)

p13 <- dst.gam.nested(d5$y, d5$X,
                      formula.null = formula.null,
                      formula.alt  = formula.alt,
                      hunt.style   = "WLS")
cat_result("alternative (tau=5, WLS)", p13)

# =============================================================================
cat("\n====  anova.gam comparison  ====\n\n")

# Helper: robust anova.gam p-value (mirrors lrt.gam logic from simulations)
anova_gam_pval <- function(y, X, formula.null, formula.alt,
                           family = "gaussian", link = NULL, k = 20) {
  fam.obj  <- if (is.null(link)) match.fun(family)() else match.fun(family)(link = link)
  dat      <- data.frame(y = y, X)
  fit.null <- mgcv::gam(stats::as.formula(formula.null), data = dat, family = fam.obj)
  fit.alt  <- mgcv::gam(stats::as.formula(formula.alt),  data = dat, family = fam.obj)
  av       <- anova(fit.null, fit.alt, test = "Chisq")
  df       <- av[2, "Df"]
  dev      <- av[2, "Deviance"]
  if (df > 0) {
    p <- av[2, "Pr(>Chi)"]
    ifelse(is.na(p), 1, p)
  } else if (dev < 0) {
    1   # null model fits better
  } else {
    0   # positive deviance gain but df <= 0: conservative p = 0
  }
}

p_anova_null <- anova_gam_pval(d0$y, d0$X, formula.null, formula.alt)
cat_result("anova.gam null (tau=0)", p_anova_null)

p_anova_alt <- anova_gam_pval(d5$y, d5$X, formula.null, formula.alt)
cat_result("anova.gam alternative (tau=5)", p_anova_alt)

cat("\n-- Side-by-side (parametric hunt) --\n\n")
cat(sprintf("%-40s  DST = %.4f   anova.gam = %.4f\n", "null (tau=0)",    p11, p_anova_null))
cat(sprintf("%-40s  DST = %.4f   anova.gam = %.4f\n", "alternative (tau=5)", p12, p_anova_alt))

# =============================================================================
cat("\n====  Both models misspecified, alt adds nothing over null  ====\n\n")
cat("DGP: y = sin(X1*X2) + noise. X1,X2,X3 independent.\n")
cat("Null: s(X1)+s(X3)  |  Alt: s(X1)+s(X3)+ti(X1,X3)\n")
cat("Alt still ignores X2 -- should NOT reject.\n\n")

n.mis <- 1000
X.mis <- matrix(runif(n.mis * 3, -1, 1), nrow = n.mis,
                dimnames = list(NULL, paste0("X", 1:3)))
y.mis <- (X.mis[, 1] * X.mis[, 2]) + rnorm(n.mis, sd = 0.5)

formula.null.mis <- "y ~ s(X1, k=10) + s(X3, k=10)"
formula.alt.mis  <- "y ~ s(X1, k=10) + s(X3, k=10) + ti(X1, X3, k=5)"

p_mis_dst   <- dst.gam.nested(y.mis, X.mis,
                               formula.null = formula.null.mis,
                               formula.alt  = formula.alt.mis)
p_mis_anova <- anova_gam_pval(y.mis, X.mis, formula.null.mis, formula.alt.mis)

cat_result("DST   (alt adds nothing)", p_mis_dst)
cat_result("anova.gam (alt adds nothing)", p_mis_anova)

cat("\n---- Done ----\n")
