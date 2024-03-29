% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcGroupPsiCWS.R
\name{CalcGroupPsiCWS}
\alias{CalcGroupPsiCWS}
\title{The Psi CWS statistic for aggregated group data}
\usage{
CalcGroupPsiCWS(data, groups, subjects)
}
\arguments{
\item{data}{A two column matrix of the bivariate pooled samples}

\item{groups}{A numeric vector of sample (or group) labels (use either 1 or 2)}

\item{subjects}{A numeric vector of subject labels}
}
\value{
The Psi CWS statistic for aggregated group data
}
\description{
This statistic computes the complementary weighted squared (CWS) differences
between the averaged subject empirical cumulative distribution functions for
the two samples. For more information, see
\insertCite{mckinney2022extensions;textual}{distdiffR} and
\insertCite{mckinney2021extensions;textual}{distdiffR}.
}
\examples{
# Randomly assign all three species to two samples
data(iris)
iris$Species <- rep(1:3, each = 50) # Species will serve as the subject label
irisPermuted <- iris[sample.int(nrow(iris)), ]
sample1 <- as.matrix(irisPermuted[1:75, c(1:2, 5)])
sample2 <- as.matrix(irisPermuted[76:150, c(1:2, 5)])
pooled_data <- rbind(cbind(sample1, 1), cbind(sample2, 2))

CalcGroupPsiCWS(pooled_data[, 1:2], pooled_data[, 4], pooled_data[, 3])
}
\references{
\insertRef{mckinney2022extensions}{distdiffR}

\insertRef{mckinney2021extensions}{distdiffR}
}
