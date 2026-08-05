## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* local macOS (aarch64), R 4.5.2
* win-builder, R-devel (2026-08-04 r90350 ucrt), x86_64-w64-mingw32

## The one note (CRAN incoming feasibility)

* "New submission" -- this is the first submission of this package.

* "Possibly misspelled words in DESCRIPTION: Debiased, Dhawan, Guo, Neyman,
  orthogonalization". These are spelled as intended: "Debiased",
  "orthogonalization" and "Neyman" (as in Neyman orthogonality) are standard
  terms in the semiparametric inference literature, and "Dhawan" and "Guo" are
  author surnames of the cited reference.

## Example timings

Several examples call machine-learning routines (regression forests via 'grf')
on samples of a few hundred observations, which is the smallest size at which
the tests being illustrated behave as documented. Those calls are wrapped in
\donttest{} so that the standard example check remains fast; the surrounding
data generation and model fitting still run.
