#' The combined rotational and toroidal shift distdiffR test for aggregated group data
#'
#' The `grouped_distdiffr()` function conducts two-sample permutation tests of
#' distributional equality based on differences in the bivariate empirical
#' cumulative density functions (BECDFs). The differences in BECDFs are computed
#' across a series of rotations, toroidal shifts, or both rotations and toroidal
#' shifts of the combined data (specified via `testType`). The number of rotations
#' and toroidal shifts may be specified (via `numRot` or `numShifts`,
#' respectively). The number of toroidal shifts may also be determined by a
#' proportion of the combined sample size (via `propPnts`).
#'
#' Additionally, `grouped_distdiffr()` assumes multiple sources (subjects) are
#' contributing to each sample. As such, the function weights each sources
#' contribution as to treat each equally within the samples, respectively.
#' However, this test is only currently available with the grouped CWS statistic
#' and employs both rotational and toroidal shifts. For more information, see
#' \insertCite{mckinney2022extensions;textual}{distdiffR} and
#' \insertCite{mckinney2021extensions;textual}{distdiffR}.
#'
#' @param aggdata1 A three column matrix of bivariate observations from one sample with the third column being the numeric subject labels
#' @param aggdata2 A three column matrix of bivariate observations from another sample with the third column being the numeric subject labels
#' @param numRot An integer number of rotational shifts of the pooled samples
#' @param propPnts A numeric proportion of points to be used as toroidal shift origins. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.
#' @param numShifts A numeric integer. The number of points to be used as toroidal shift origins. Must be less than the pooled sample size. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.
#' @param shiftThrshld A numeric integer. Used if neither propPnts or numShifts are provided. If the pooled sample size is less than shiftThrshld, every point will be used as a toroidal shift origin. Otherwise, only a random sample of shiftThrshld points will be used.
#' @param numPerms An integer number of permutations of the original data
#' @param psiFun A function specifying the Psi statistic calculation. Default is the CalcGroupPsiCWS.
#' @param seedNum An integer random seed value
#'
#' @return A list including three objects:
#'     (1) the Psi statistic computed on the original data
#'     (2) a vector of Psi statistics computed on the permuted data
#'     (3) the p-value for the test
#' @importFrom stats runif
#' @importFrom stats median
#' @export
#' @references
#' \insertRef{mckinney2022extensions}{distdiffR}
#'
#' \insertRef{mckinney2021extensions}{distdiffR}
#' @examples
#' # Randomly assign all three species to two samples
#' # The species serve as subject labels within each sample
#' seedNum <- 123
#' set.seed(seedNum)
#'
#' data(iris)
#' irisPermuted <- iris
#' irisPermuted$Species <- rep(1:3, each = 50)
#' irisPermuted <- irisPermuted[sample.int(nrow(irisPermuted)), ]
#' sample1 <- as.matrix(irisPermuted[1:75, c(1:2, 5)])
#' sample2 <- as.matrix(irisPermuted[76:150, c(1:2, 5)])
#'
#' output <- grouped_distdiffr(sample1,
#'                             sample2,
#'                             seedNum = seedNum)
#' output$pval
grouped_distdiffr <- function(aggdata1,
                              aggdata2,
                              numRot = 8,
                              propPnts = NULL,
                              numShifts = NULL,
                              shiftThrshld = 25,
                              numPerms = 999,
                              psiFun = CalcGroupPsiCWS,
                              seedNum = NULL) {
  # aggdata1 and aggdata2 are matrices where the three columns represent X, Y, and subjectNumber.
  subjNums1 <- aggdata1[, 3]
  subjNums2 <- aggdata2[, 3]

  aggdata1 <- aggdata1[, 1:2]
  aggdata2 <- aggdata2[, 1:2]

  # Combines the data from two subjects into one long matrix.
  n1 <- nrow(aggdata1)
  n2 <- nrow(aggdata2)

  # Check for conflict between propPnts and numShifts.
  noPropPnts <- is.null(propPnts)
  noNumShifts <- is.null(numShifts)
  n_pooled <- n1 + n2
  if (!noPropPnts & !noNumShifts) {
    stop("Must provide either propPnts or numShifts, but not both.")
  } else if (noPropPnts & noNumShifts) { # Use either numShifts or set using shiftThrshld.
    numShifts <- ifelse(n_pooled > shiftThrshld, shiftThrshld, n_pooled)
    noNumShifts <- FALSE
  }

  hash1 <- hashMat(aggdata1)
  hash2 <- hashMat(aggdata2)

  if (hash1 >= hash2) {
    data <- rbind(aggdata1, aggdata2)
    groupNums <- rep(1:2, times = c(n1, n2))
    subjNums <- c(subjNums1, subjNums2)
  } else {
    data <- rbind(aggdata2, aggdata1)
    groupNums <- rep(2:1, times = c(n2, n1))
    subjNums <- c(subjNums2, subjNums1)
  }

  set.seed(seedNum)

  ## Center the data around the bivariate median of the combined data sets.
  medians <- apply(data, 2, median)
  data <- sweep(data, 2, medians)

  ## Rotate the data and stores the rotated data frames in a list.
  rotDataList <- RotateData(data, numRot)

  ## Applies toroidal shifts to the data and store the shifted data frames in a list.
  if (!noPropPnts) {
    lstOfRotShiftDataLists <- lapply(rotDataList, PropToroShiftData, n1, n2, propPnts)
  } else if (!noNumShifts) {
    lstOfRotShiftDataLists <- lapply(rotDataList, NumToroShiftData, n1, n2, numShifts)
  }

  ## Calculate psi for the real data
  truePsi <- mean(sapply(lstOfRotShiftDataLists, function(rotDataLst) mean(sapply(rotDataLst, psiFun, groupNums, subjNums))))

  # Calculate psi for all permutations of the rotated and toroidal shifted data.
  permPsi <- rep(0, numPerms)
  for (i in 1:numPerms) {
    permGroups <- sample(groupNums, length(groupNums), replace = FALSE)
    permPsi[i] <- mean(sapply(lstOfRotShiftDataLists, function(rotDataLst) mean(sapply(rotDataLst, psiFun, permGroups, subjNums))))
  }

  list(psiStat = truePsi,
       permPsi = permPsi,
       pval = mean(c((permPsi >= truePsi), 1)))
}
