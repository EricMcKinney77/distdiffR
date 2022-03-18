# Check that the distdiffR tests are producing reasonable outputs.

if (!require(devtools)) install.packages("devtools")
devtools::install_github("https://github.com/EricMcKinney77/distdiffR", build_vignettes = TRUE)
library(distdiffR)

seedNum <- 123
set.seed(seedNum)

# Test when the null is true (both distributions are equivalent)
data(iris)
# Randomly assign all three species to two samples
irisPermuted <- iris[sample.int(nrow(iris)), ]
sample1 <- as.matrix(irisPermuted[1:75, -(3:5)])
sample2 <- as.matrix(irisPermuted[76:150, -(3:5)])

# Rotational test
output <- distdiffr(sample1,
                    sample2,
                    testType = "rotational",
                    seedNum = seedNum)
output$pval == 0.58

# Order does not matter
output2 <- distdiffr(sample2,
                     sample1,
                     testType = "rotational",
                     seedNum = seedNum)
output2$pval == output$pval

# Toroidal shift test with proportions of points
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    propPnts = 0.1,
                    seedNum = seedNum)
output$pval == 0.449

# Toroidal shift test with thresholds below pooled sample size
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    seedNum = seedNum)
output$pval == 0.5

# Toroidal shift test with thresholds above pooled sample size
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    shiftThrshld = 100,
                    seedNum = seedNum)
output$pval == 0.487

# Toroidal shift test with a number of shifts
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    numShifts = 8,
                    seedNum = seedNum)
output$pval == 0.531

# Combined rotational and toroidal shift test
output <- distdiffr(sample1,
                    sample2,
                    testType = "combined",
                    seedNum = seedNum)
output$pval == 0.331


# Test when the null is false
data(iris)
sample1 <- as.matrix(iris[iris$Species == "setosa", 1:2])
sample2 <- as.matrix(iris[iris$Species == "virginica", 1:2])
pooled_data <- rbind(sample1, sample2)
n1 <- nrow(sample1)
n2 <- nrow(sample2)

output <- NumToroShiftData(pooled_data, n1, n2, 25)

sample_labels <- rep(1:2, c(sum(iris[5] == "setosa"),
                            sum(iris[5] == "virginica")))

CalcPsiCWS(as.matrix(pooled_data), sample_labels)
# Rotational test
output <- distdiffr(sample1,
                    sample2,
                    testType = "rotational",
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with proportions of points
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    propPnts = 0.1,
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with thresholds below pooled sample size
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with thresholds above pooled sample size
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    shiftThrshld = 100,
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with a number of shifts
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    numShifts = 8,
                    seedNum = seedNum)
output$pval == 0.001

# Combined rotational and toroidal shift test
output <- distdiffr(sample1,
                    sample2,
                    testType = "combined",
                    seedNum = seedNum)
output$pval == 0.001



# Randomly assign all three species to two samples
data(iris)
sample1 <- as.matrix(iris[iris$Species == "setosa", 1:2])

# Generate five rotated versions of sample1 (every 72 degrees) within 360 degrees.
RotateData(sample1, 5)


data(iris)
sample1 <- as.matrix(iris[iris$Species == "setosa", 1:2])
iris$Species <- rep(1:3, each = 50) # Species will serve as the subject label
irisPermuted <- iris[sample.int(nrow(iris)), ]
sample1 <- as.matrix(irisPermuted[1:75, c(1:2, 5)])
sample2 <- as.matrix(irisPermuted[76:150, c(1:2, 5)])
pooled_data <- rbind(cbind(sample1, 1), cbind(sample2, 2))

CalcGroupPsiCWS(pooled_data[, 1:2], pooled_data[, 4], pooled_data[, 3])


data(iris)
sample1 <- as.matrix(iris[iris$Species == "setosa", 1:2])
sample2 <- as.matrix(iris[iris$Species == "virginica", 1:2])
pooled_data <- rbind(sample1, sample2)
n1 <- nrow(sample1)
n2 <- nrow(sample2)

# Create five
output <- NumToroShiftData(pooled_data, n1, n2, 5)
output
output <- PropToroShiftData(pooled_data, n1, n2, 0.1)
output


data(iris)
sample1 <- as.matrix(iris[iris$Species == "virginica", 1:2])
sample2 <- as.matrix(iris[iris$Species == "versicolor", 1:2])

hashMat(sample1)
hashMat(sample2)

bcdf(sample1, sample2)
