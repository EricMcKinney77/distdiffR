# distdiffR

An R package for bivariate two-sample tests of distributional equality.

## Overview

This R package provides a collection of bivariate two-sample tests for distributional equality. The tests make use of statistics between empirical cumulative distribution functions averaged across a series of rotations and or toroidal shifts of the pooled samples. The variety of tests with their respective parameters can be called from the main function `distdiffr()`. It takes as input two bivariate samples in the form of two-column matrices, referred to as `data1` and `data2`.

Another version of the test also exists when combining bivariate data from multiple sources into each of the two-samples, respectively, which treats each contribution equally. This can be called via the `grouped_distdiffr()` function.

## Installation

You can install the development version of dplyr from GitHub.

```
install.packages("devtools")
devtools::install_github("EricMcKinney77/distdiffR")
```

## Usage

```
seedNum <- 123
set.seed(seedNum)

# Test when the null is true (both distributions are equivalent)
data(iris)
# Randomly assign all three species to two samples
iris <- iris[sample.int(nrow(iris)), ]
data1 <- as.matrix(iris[1:75, -(3:5)])
data2 <- as.matrix(iris[76:150, -(3:5)])

# Rotational test
output <- distdiffr(data1,
                    data2,
                    testType = "rotational",
                    seedNum = seedNum)
output$pval

# Toroidal shift test
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    seedNum = seedNum)
output$pval

# Combined rotational and toroidal shift test
output <- distdiffr(data1,
                    data2,
                    testType = "combined",
                    seedNum = seedNum)
output$pval


# Test when the null is false
data(iris)
data1 <- as.matrix(iris[iris[5] == "setosa", -(3:5)])
data2 <- as.matrix(iris[iris[5] == "virginica", -(3:5)])

# Rotational test
output <- distdiffr(data1,
                    data2,
                    testType = "rotational",
                    seedNum = seedNum)
output$pval

# Toroidal shift test
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    seedNum = seedNum)
output$pval

# Combined rotational and toroidal shift test
output <- distdiffr(data1,
                    data2,
                    testType = "combined",
                    seedNum = seedNum)
output$pval
```

Please see `vignette("distdiffR")` for more examples and greater explanation of the implementation of `distdiffR` test functions.

## Providing Feedback

If you encounter a bug, please file an issue with a minimal reproducible example on [GitHub](https://github.com/EricMcKinney77/distdiffR/issues). Thank you!

---
