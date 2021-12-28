# distdiffR

An R package for bivariate two-sample tests of distributional equality.

## Overview

This R package provides a collection of bivariate two-sample tests for distributional equality. The tests make use of statistics between empirical cumulative distribution functions averaged across a series of rotations and or toroidal shifts of the pooled samples. The variety of tests with their respective parameters can be called from the main function `distdiffr()`. It takes as input two bivariate samples (not necessarily the same size) in the form of two-column matrices, referred to as `data1` and `data2`. The data inputs order is arbitrary.

Another version of the test also exists when combining bivariate data from multiple sources into each of the two-samples, respectively, which treats each contribution equally. This can be called by setting the `group` parameter within `distdiffr()` to `TRUE` and providing a vector, say `X` of subject labels, e.g., `distdiffr(groupLabels = X, group = TRUE)`.

## Installation

You can install the development version of dplyr from GitHub.

`# install.packages("devtools")`
`devtools::install_github("EricMcKinney77/distdiffR")`

## Usage

Please see `vignette("distdiffR")` for more examples and greater explanation of the implementation of `distdiffR` test functions.

## Providing Feedback

If you encounter a bug, please file an issue with a minimal reproducible example on [GitHub](https://github.com/EricMcKinney77/distdiffR/issues). Thank you!

---
