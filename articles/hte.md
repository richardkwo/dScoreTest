# Effect modification in HIV treatment (ACTG175)

``` r

library(dScoreTest)
```

We use
[`hte_test_conditional()`](https://unbiased.co.in/dScoreTest/reference/hte_test_conditional.md)
to look for effect modifiers in treating HIV, using the AIDS Clinical
Trials Group (ACTG) Protocol 175 data (Hammer et al., 1996), available
from the R package `speff2trial`. This is a randomized trial that
compares treatments for HIV among adults with CD4 T-cell counts between
200 and 500 cells per \\\text{mm}^3\\. We consider two treatment arms:
monotherapy with didanosine (\\T=0\\) versus combination therapy with
zidovudine and didanosine (\\T=1\\), which have 561 and 522 patients
respectively. We study the CD4 count at \\20 \pm 5\\ weeks as the
continuous outcome \\Y\\, with higher values indicating better health.

There are 12 baseline covariates. Here, we look at potential effect
modification by *homosexual activity* (binary) and *baseline CD4 count*
(continuous).

``` r

data(ACTG175, package = "speff2trial")

covs <- c(
    "age", "wtkg", "karnof", "cd40", "cd80", "gender",
    "homo", "race", "symptom", "drugs", "str2", "hemo"
)
label_list <- c(
    cd40 = "CD4", homo = "Homosexual activity", age = "Age", wtkg = "Weight",
    str2 = "Antiretroviral history", race = "Race", hemo = "Hemophilia",
    symptom = "Symptomatic", cd80 = "CD8", karnof = "Karnofsky score",
    drugs = "IV drug use", gender = "Gender"
)

dat <- ACTG175[ACTG175$arms %in% c(1, 3), ]
Z  <- dat[, covs]
Tr <- as.integer(dat$arms == 1)
y  <- dat$cd420
```

For a covariate \\Z_i\\, we use `hte_test_conditional(y, Tr, Z, S = i)`
to test the conditional hypothesis \\H_0(\\i\\)\\: the conditional
treatment effect \\ \tau(Z) := \mathbb{E}\[Y(t=1) - Y(t=0) \mid Z\]\\
does not depend on \\Z_i\\, while holding all the remaining covariates
fixed. Although here \\S\\ is chosen to be a singleton, the method is
also applicable where \\S\\ is a group of covariates. The test involves
a random sample split, so we repeat it a number of times and combine the
p-values by their harmonic mean, which controls the type-I error if the
\\\alpha\\ level is small in the sense that \\\text{Type-I error} /
\alpha \to 1\\ holds as \\\alpha \to 0\\ under the null, even when the
\\p\\-values are computed on the same sample and are therefore
dependent.

``` r

hmean <- function(p) length(p) / sum(1 / p)

n.reps <- 10
covariates.of.interest <- c("homo", "cd40")

set.seed(42)
results <- do.call(rbind, lapply(covariates.of.interest, function(cov) {
    cov.idx <- which(covs == cov)
    p.vals <- replicate(n.reps, hte_test_conditional(y, Tr, Z, S = cov.idx)$p.val)
    data.frame(covariate = cov, label = label_list[[cov]],
              n.reps = n.reps, harmonic.p = hmean(p.vals))
}))
rownames(results) <- NULL
results
#>   covariate               label n.reps harmonic.p
#> 1      homo Homosexual activity     10 0.02219521
#> 2      cd40                 CD4     10 0.02535719
```

Both covariates come back with harmonic p-values below the conventional
0.05 threshold, suggesting that both homosexual activity and baseline
CD4 count may modify the effect of combination therapy versus
monotherapy on CD4 count at 20 weeks. This is only a preview with 10
reps per covariate, so the harmonic-mean estimate is still noisy; a
production run would use many more reps for a stable estimate.
