#' The combined rotational and toroidal shift distDiffR test for aggregated group data
#'
#' @param aggdata1 A three column matrix of bivariate observations from one sample with the third column being the subject labels
#' @param aggdata2 A three column matrix of bivariate observations from another sample with the third column being the subject labels
#' @param numRot An integer number of rotational shifts of the pooled samples
#' @param propPnts A numeric proportion of points to be used as toroidal shift origins
#' @param numShifts A numeric integer. The number of points to be used as toroidal shift origins. Must be less than the pooled sample size.
#' @param numPerms An integer number of permutations of the original data
#' @param psiFun A function specifying the Psi statistic calculation. Default is the CalcGroupPsiRWS.
#' @param seedNum An integer random seed value
#'
#' @return A list including three objects:
#'     (1) the Psi statistic computed on the original data
#'     (2) a vector of Psi statistics computed on the permuted data
#'     (3) the p-value for the test
#' @importFrom stats runif
#' @importFrom stats median
#' @export
AggregatedRotateToroDiffTest <- function(aggdata1, aggdata2, numRot = 8, propPnts = NULL, numShifts = NULL,
                                         numPerms = 99, psiFun = CalcGroupPsiRWS, seedNum = NULL) {
  noPropPnts <- is.null(propPnts)
  noNumShifts <- is.null(numShifts)
  if (noPropPnts & noNumShifts) {
    stop("Must provide propPnts or numShifts.")
  } else if (!noPropPnts & !noNumShifts) {
    stop("Must provide either propPnts or numShifts, but not both.")
  }

  # aggdata1 and aggdata2 are matrices where the three columns represent X, Y, and subjectNumber.
  subjNums1 <- aggdata1[, 3]
  subjNums2 <- aggdata2[, 3]
  subjNums <- c(subjNums1, subjNums2)

  aggdata1 <- aggdata1[, 1:2]
  aggdata2 <- aggdata2[, 1:2]

  # NOTE: Data cleaning must be done before applying this function, e.g., filter(X != 0 & Y != 0)
  set.seed(seedNum)

  ## Combines the data from two groups into one long matrix.
  n1 <- nrow(aggdata1)
  n2 <- nrow(aggdata2)
  groupNums <- rep(1:2, times = c(n1, n2))

  # Check and remove duplicate data values between samples
  aggdata1DuplRwsWaggdata2 <- function(aggdata1, aggdata2) {
    aggdata1strgs <- unlist(sapply(1:n1, function(j) paste(aggdata1[j, ], collapse = '_')))
    aggdata2strgs <- unlist(sapply(1:n2, function(j) paste(aggdata2[j, ], collapse = '_')))
    which(aggdata1strgs %in% intersect(aggdata1strgs, aggdata2strgs))
  }
  duplRwsInd <- aggdata1DuplRwsWaggdata2(aggdata1, aggdata2)

  if (length(duplRwsInd) != 0) {
    print("Some data values were not unique between groups. Adding an insignificant amount of noise.")

    # Determine the number of significant digits that were used to record the data
    decimalplaces <- function(x) {
      if (abs(x - round(x)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', sprintf(fmt = "%f", x)), ".", fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }
    numSigDigits <- max(mapply(decimalplaces, rbind(aggdata1, aggdata2)))
    roundBound <- 5 * 10^(-(numSigDigits + 1))

    # Keep adding small amounts of noise until aggdata1 and aggdata2 have no more duplicate rows
    while (length(duplRwsInd) != 0) {
      aggdata1[duplRwsInd, ] <- aggdata1[duplRwsInd, ] + runif(length(duplRwsInd) * 2, -roundBound, roundBound)
      duplRwsInd <- aggdata1DuplRwsWaggdata2(aggdata1, aggdata2)
    }
  }

  data <- rbind(aggdata1, aggdata2)
  srtDatOrdr <- order(data[, 1], data[, 2])
  data <- data[srtDatOrdr, ]
  subjNums <- subjNums[srtDatOrdr]
  groupNums <- groupNums[srtDatOrdr]

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
