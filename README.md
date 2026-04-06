<!-- badges: start -->
[![R-CMD-check](https://github.com/juliusolaifa/mixchisq/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/juliusolaifa/mixchisq/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# mixchisq

`mixchisq` is a lightweight R package for working with finite mixtures of chi-square distributions, including chi-bar-square style laws with an atom at zero.

It provides:

- `dmixchisq()` for densities,
- `pmixchisq()` for cumulative probabilities,
- `qmixchisq()` for numerical inversion,
- `rmixchisq()` for random generation,
- `plot_mixchisq()` for quick visual inspection.

## Why this package?

In constrained likelihood ratio testing, the asymptotic null distribution is often not a single chi-square law. Instead, it can be a weighted mixture such as:

- `0.5 * chi^2_0 + 0.5 * chi^2_1`,
- `0.5 * chi^2_2 + 0.5 * chi^2_3`,
- or a more general chi-bar-square distribution.

This package gives you a compact utility layer for those cases.

## Installation

```r
# local install from source directory
# install.packages("devtools")
devtools::install_local("mixchisq")
```

## Basic usage

```r
library(mixchisq)

# A simple two-component mixture
pmixchisq(5, df = c(2, 3), w = c(0.5, 0.5))
dmixchisq(5, df = c(2, 3), w = c(0.5, 0.5))
qmixchisq(0.95, df = c(2, 3), w = c(0.5, 0.5))
rmixchisq(10, df = c(2, 3), w = c(0.5, 0.5))
```

## Mixtures with an atom at zero

```r
library(mixchisq)

# Chi-bar-square style law
# 0.5 * chi^2_0 + 0.5 * chi^2_1
pmixchisq(c(0, 0.5, 2), df = c(0, 1), w = c(0.5, 0.5))
qmixchisq(c(0.25, 0.5, 0.95), df = c(0, 1), w = c(0.5, 0.5))
set.seed(123)
rmixchisq(20, df = c(0, 1), w = c(0.5, 0.5))
```

## Quick plotting

```r
plot_mixchisq(df = c(2, 3), w = c(0.5, 0.5), type = "cdf")
plot_mixchisq(df = c(0, 1), w = c(0.5, 0.5), type = "cdf")
```
