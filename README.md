# `dst` — Debiased Score Tests for Nonparametric Goodness of Fit

A framework for goodness-of-fit testing where you fit a null model, "hunt" non-parametrically for a direction of misspecification (default via GRF), then form a debiased test statistic via sample splitting to recover correct asymptotic size.

## Core engine — `R/dst.R`

- **`debiased.scoretest(y, X, score.fun, reg.method, hunt.method, proj.method, ...)`** (exported) — the generic test. Glue layer: takes user-supplied score / fit / hunt / projection callbacks and returns a one-sided p-value. Supports a 2-way split (`in.sample = TRUE`) or 3-way split, with optional multi-split aggregation via `MultiSplit::test.multisplit`.
- **`.dst.single`** (internal) — runs one split: fit null on aux, hunt direction `h.hat`, refit + project on main half, build statistic `T = sqrt(n) * mean(L) / sd(L)` from `L = resid * (h.hat - m.h.hat)`.

## Hunting — `R/hunt.R`

- **`hunt.nonpara(resids, X, y, fit, ...)`** (exported) — canonical non-parametric hunter using `grf::regression_forest`. Three styles: `"vanilla"` (regress residuals on X), `"WLS"` (regress 1/resid weighted by resid^2), `"opt.WLS"` (additionally estimates the conditional variance and projects out its component for an asymptotically optimal direction). Tukey IQR clipping on outputs.
- `.rf` / `.rf_predict` — internal wrappers that handle categorical encoding (rank by weighted mean response) before calling `grf`.

## GAM tests — `R/gam.R`

- **`dst.gam(y, X, formula, family, link, ...)`** (exported) — wires `mgcv::gam`/`bam` into the generic test with non-parametric hunting.
- **`dst.gam.nested(y, X, formula.null, formula.alt, ...)`** (exported) — `anova`-style nested GAM comparison. If `formula.alt` supplied, hunts parametrically using the alt model as the learner; otherwise falls back to GRF.
- Internals: `.gam_fit`, `.gam_score` (response-scale residuals), `.gam_proj` (weighted GAM regression of `h.hat`), `.gam_weight_fit` (link derivative), `.gam_hunt` (parametric hunt mirroring `hunt.nonpara` logic), `.gam_proj_hunt` factory.

## Quantile GAM test — `R/qgam.R`

- **`dst.qgam(y, X, formula, qu, ...)`** (exported) — uses `qgam::qgam` for the null. Score is the quantile check derivative `1{y <= q.hat} - qu`. Weights/projections rely on a conditional density `f(q.hat(X)|X)`.
- Internals: `.qgam_fit`, `.qgam_score`, `.qgam_proj`, `.qgam_cond_density` — quantile-forest weights x Gaussian kernel with CV-tuned bandwidth (Silverman scaled).

## Heterogeneous treatment effects — `R/hte.R`

- **`dst.hte(y, Z, D, het.vars, binary.treatment, ...)`** (exported) — tests heterogeneity of the treatment effect in a partially linear model. Null permits heterogeneity in `setdiff(all, het.vars)`; hunt searches over all covariates with a T-learner.
- Internals: `.hte_combine`/`.hte_split` (D in column 1), `.hte_fit`, `.hte_score`, `.hte_proj`, `.hte_hunt` (T-learner forests for treated/control with optional WLS / opt.WLS weighting), and `.plm_grf` — DML estimator using `grf::boosted_regression_forest` with optional cross-fitting; returns `predict` and `propensity` closures.

## Pattern across families

Every family-specific exported function (`dst.gam`, `dst.qgam`, `dst.hte`, `dst.gam.nested`) is a thin wrapper that builds the four callbacks (`reg.method`, `score.fun`, `hunt.method`, `proj.method`) and passes them to `debiased.scoretest`. The hunting functions all share a common skeleton: fit a base learner of `1/resids` weighted by `resids^2`, estimate a variance function, project onto the null space, then trim outliers via Tukey's IQR.
