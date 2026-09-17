# Check that the distdiffR tests are producing reasonable outputs.

library(distdiffR)

cat("Starting distdiffR sanity checks...\n")


seedNum <- 123
set.seed(seedNum)

# Test when the null is true (both distributions are equivalent)
cat("\nTesting cases where Null is True (distributions are equivalent)...\n")
data(iris)
# Randomly assign all three species to two samples
irisPermuted <- iris[sample.int(nrow(iris)), ]
sample1 <- as.matrix(irisPermuted[1:75, -(3:5)])
sample2 <- as.matrix(irisPermuted[76:150, -(3:5)])

# Rotational test
output <- distdiffr(sample1,
    sample2,
    testType = "rotational",
    seedNum = seedNum
)
stopifnot(output$pval == 0.754)


# Order does not matter
output2 <- distdiffr(sample2,
    sample1,
    testType = "rotational",
    seedNum = seedNum
)
stopifnot(output2$pval == output$pval)

# Toroidal shift test with proportions of points
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    propPnts = 0.1,
    seedNum = seedNum
)
stopifnot(output$pval == 0.272)

# Toroidal shift test with thresholds below pooled sample size
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    seedNum = seedNum
)
stopifnot(output$pval == 0.205)

# Toroidal shift test with thresholds above pooled sample size
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    shiftThrshld = 100,
    seedNum = seedNum
)
stopifnot(output$pval == 0.494)

# Toroidal shift test with a number of shifts
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    numShifts = 8,
    seedNum = seedNum
)
stopifnot(output$pval == 0.206)

# Combined rotational and toroidal shift test
output <- distdiffr(sample1,
    sample2,
    testType = "combined",
    seedNum = seedNum
)
stopifnot(output$pval == 0.267)


# Test when the null is false
cat("\nTesting cases where Null is False (distributions differ)...\n")
data(iris)
sample1 <- as.matrix(iris[iris$Species == "setosa", 1:2])
sample2 <- as.matrix(iris[iris$Species == "virginica", 1:2])
pooled_data <- rbind(sample1, sample2)
n1 <- nrow(sample1)
n2 <- nrow(sample2)

# Rotational test
output <- distdiffr(sample1,
    sample2,
    testType = "rotational",
    seedNum = seedNum
)
stopifnot(output$pval == 0.001)

# Toroidal shift test with proportions of points
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    propPnts = 0.1,
    seedNum = seedNum
)
stopifnot(output$pval == 0.001)

# Toroidal shift test with thresholds below pooled sample size
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    seedNum = seedNum
)
stopifnot(output$pval == 0.001)

# Toroidal shift test with thresholds above pooled sample size
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    shiftThrshld = 100,
    seedNum = seedNum
)
stopifnot(output$pval == 0.001)

# Toroidal shift test with a number of shifts
output <- distdiffr(sample1,
    sample2,
    testType = "toroidal",
    numShifts = 8,
    seedNum = seedNum
)
stopifnot(output$pval == 0.001)

# Combined rotational and toroidal shift test
output <- distdiffr(sample1,
    sample2,
    testType = "combined",
    seedNum = seedNum
)
stopifnot(output$pval == 0.001)



# Testing helper functions
cat("\nTesting helper functions (Hashing, BCDF, Toroidal Shifts)...\n")
# Toroidal shifts
output_num <- NumToroShiftData(pooled_data, n1, n2, 5)
stopifnot(length(output_num) == 5)

output_prop <- PropToroShiftData(pooled_data, n1, n2, 0.1)
stopifnot(length(output_prop) == ceiling(0.1 * (n1 + n2)))

# Hashing and BCDF
stopifnot(is.numeric(hashMat(sample1)))
stopifnot(is.numeric(hashMat(sample2)))
stopifnot(is.numeric(bcdf(sample1, sample2)))

cat("\nAll sanity checks passed successfully!\n")

