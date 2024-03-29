% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distdiffr.R
\name{distdiffr}
\alias{distdiffr}
\title{The distdiffR two-sample tests of bivariate distributional equality}
\usage{
distdiffr(
  data1,
  data2,
  testType = "combined",
  numRot = 8,
  propPnts = NULL,
  numShifts = NULL,
  shiftThrshld = 25,
  numPerms = 999,
  psiFun = CalcPsiCWS,
  seedNum = NULL
)
}
\arguments{
\item{data1}{A two column matrix of bivariate observations from one sample.}

\item{data2}{A two column matrix of bivariate observations from another sample.}

\item{testType}{A string indicating the type of test to be used. Must be one of c("rotational", "toroidal", "combined").}

\item{numRot}{An integer number of rotational shifts of the pooled samples.}

\item{propPnts}{A numeric proportion of points to be used as toroidal shift origins. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.}

\item{numShifts}{A numeric integer. The number of points to be used as toroidal shift origins. Must be less than the pooled sample size. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.}

\item{shiftThrshld}{A numeric integer. Used if neither propPnts or numShifts are provided. If the pooled sample size is less than shiftThrshld, every point will be used as a toroidal shift origin. Otherwise, only a random sample of shiftThrshld points will be used.}

\item{numPerms}{An integer number of permutations of the original data.}

\item{psiFun}{A function specifying the Psi statistic calculation.}

\item{seedNum}{An integer random seed value.}
}
\value{
A list including three objects:
(1) the Psi statistic computed on the original data
(2) a vector of Psi statistics computed on the permuted data
(3) the p-value for the test
}
\description{
The \code{distdiffr()} function conducts two-sample permutation tests of
distributional equality based on differences in the bivariate empirical
cumulative density functions (BECDFs). The differences in BECDFs are computed
across a series of rotations, toroidal shifts, or both rotations and toroidal
shifts of the combined data (specified via \code{testType}). The number of rotations
and toroidal shifts may be specified (via \code{numRot} or \code{numShifts},
respectively). The number of toroidal shifts may also be determined by a
proportion of the combined sample size (via \code{propPnts}). However,
\insertCite{mckinney2022extensions;textual}{distdiffR} has shown that
limiting the number of toroidal shifts to ease the computational load of the
test will still provide stable results. Simulations have shown the combined
rotational and toroidal shift test to be the most powerful yet appropriately
conservative test. For more information, see
\insertCite{mckinney2022extensions;textual}{distdiffR} and
\insertCite{mckinney2021extensions;textual}{distdiffR}.
}
\examples{
# Randomly assign all three species to two samples
seedNum <- 123
set.seed(seedNum)

data(iris)
# Randomly assign all three species to two samples
irisPermuted <- iris[sample.int(nrow(iris)), ]
sample1 <- as.matrix(irisPermuted[1:75, 1:2])
sample2 <- as.matrix(irisPermuted[76:150, 1:2])
pooled_data <- rbind(cbind(sample1, 1), cbind(sample2, 2))

# Rotational test
output <- distdiffr(sample1, # Note: Data inputs must be matrices
                    sample2,
                    testType = "rotational",
                    numRot = 8, # Default value
                    seedNum = seedNum)
output$pval

# Toroidal shift test with proportions of points
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    propPnts = 0.1,
                    seedNum = seedNum)
output$pval

# Toroidal shift test with a threshold below pooled sample size
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    shiftThrshld = 25, # Default
                    seedNum = seedNum)
output$pval

# Toroidal shift test with a threshold above pooled sample size
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    shiftThrshld = 200,
                    seedNum = seedNum)
output$pval

# Toroidal shift test with a number of shifts
output <- distdiffr(sample1,
                    sample2,
                    testType = "toroidal",
                    numShifts = 8,
                    seedNum = seedNum)
output$pval

# Combined rotational and toroidal shift test
output <- distdiffr(sample1,
                    sample2,
                    testType = "combined", # Default
                    numRot = 8,            # Default
                    shiftThrshld = 25,     # Default
                    seedNum = seedNum)
output$pval

# Also see browseVignettes(package = "distdiffR")
}
\references{
\insertRef{mckinney2022extensions}{distdiffR}

\insertRef{mckinney2021extensions}{distdiffR}
}
