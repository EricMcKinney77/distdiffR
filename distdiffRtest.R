# Check that the distdiffR tests are producing reasonable outputs.

library(devtools)
install_github("https://github.com/EricMcKinney77/distdiffR")
library(distdiffR)

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
output$pval == 0.561

# Order does not matter
output2 <- distdiffr(data2,
                     data1,
                     testType = "rotational",
                     seedNum = seedNum)
output2$pval == output$pval

# Toroidal shift test with proportions of points
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    propPnts = 0.1,
                    seedNum = seedNum)
output$pval == 0.566

# Toroidal shift test with thresholds below pooled sample size
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    seedNum = seedNum)
output$pval == 0.477

# Toroidal shift test with thresholds above pooled sample size
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    shiftThrshld = 100,
                    seedNum = seedNum)
output$pval == 0.502

# Toroidal shift test with a number of shifts
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    numShifts = 8,
                    seedNum = seedNum)
output$pval == 0.648

# Combined rotational and toroidal shift test
output <- distdiffr(data1,
                    data2,
                    testType = "combined",
                    seedNum = seedNum)
output$pval == 0.307


# Test when the null is false
data(iris)
data1 <- as.matrix(iris[iris[5] == "setosa", -(3:5)])
data2 <- as.matrix(iris[iris[5] == "virginica", -(3:5)])

# Rotational test
output <- distdiffr(data1,
                    data2,
                    testType = "rotational",
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with proportions of points
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    propPnts = 0.1,
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with thresholds below pooled sample size
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with thresholds above pooled sample size
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    shiftThrshld = 100,
                    seedNum = seedNum)
output$pval == 0.001

# Toroidal shift test with a number of shifts
output <- distdiffr(data1,
                    data2,
                    testType = "toroidal",
                    numShifts = 8,
                    seedNum = seedNum)
output$pval == 0.001

# Combined rotational and toroidal shift test
output <- distdiffr(data1,
                    data2,
                    testType = "combined",
                    seedNum = seedNum)
output$pval == 0.001
