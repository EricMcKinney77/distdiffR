#' The combined rotational and toroidal shift distDiffR test
#'
#' @param data1 a two column matrix of bivariate observations from one sample
#' @param data2 a two column matrix of bivariate observations from another sample
#' @param numRot an integer number of rotational shifts of the pooled samples
#' @param propPnts A numeric proportion of points to be used as toroidal shift origins
#' @param numShifts A numeric integer. The number of points to be used as toroidal shift origins. Must be less than the pooled sample size.
#' @param numPerms an integer number of permutations of the original data
#' @param psiFun a function specifying the Psi statistic calculation
#' @param seedNum an integer random seed value
#'
#' @return A list including three objects:
#'     (1) the Psi statistic computed on the original data
#'     (2) a vector of Psi statistics computed on the permuted data
#'     (3) the p-value for the test
#' @importFrom stats runif
#' @importFrom stats median
#' @export
RotateToroDiffTest <- function(data1, data2, numRot = 8, propPnts = NULL, numShifts = NULL,
                               numPerms = 999, psiFun = CalcPsiRWS, seedNum = NULL) {
  noPropPnts <- is.null(propPnts)
  noNumShifts <- is.null(numShifts)
  if (noPropPnts & noNumShifts) {
    stop("Must provide propPnts or numShifts.")
  } else if (!noPropPnts & !noNumShifts) {
    stop("Must provide either propPnts or numShifts, but not both.")
  }

  # NOTE: Data cleaning must be done before applying this function, e.g., filter(X != 0 & Y != 0)
  set.seed(seedNum)

  ## Combines the data from two subjects into one long matrix.
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  subjects <- rep(1:2, times = c(n1, n2))

  # Check and remove duplicate data values between samples
  data1DuplRwsWdata2 <- function(data1, data2) {
    data1strgs <- unlist(sapply(1:n1, function(j) paste(data1[j, ], collapse = '_')))
    data2strgs <- unlist(sapply(1:n2, function(j) paste(data2[j, ], collapse = '_')))
    which(data1strgs %in% intersect(data1strgs, data2strgs))
  }
  duplRwsInd <- data1DuplRwsWdata2(data1, data2)

  if (length(duplRwsInd) != 0) {
    print("Some data values were not unique between subjects. Adding an insignificant amount of noise.")

    # Determine the number of significant digits that were used to record the data
    decimalplaces <- function(x) {
      if (abs(x - round(x)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', sprintf(fmt = "%f", x)), ".", fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }
    numSigDigits <- max(mapply(decimalplaces, rbind(data1, data2)))
    roundBound <- 5 * 10^(-(numSigDigits + 1))

    # Keep adding small amounts of noise until data1 and data2 have no more duplicate rows
    while (length(duplRwsInd) != 0) {
      data1[duplRwsInd, ] <- data1[duplRwsInd, ] + runif(length(duplRwsInd) * 2, -roundBound, roundBound)
      duplRwsInd <- data1DuplRwsWdata2(data1, data2)
    }
  }

  data <- rbind(data1, data2)
  srtDatOrdr <- order(data[, 1], data[, 2])
  data <- data[srtDatOrdr, ]
  subjects <- subjects[srtDatOrdr]

  ## Center the data around the bivariate median of the combined data sets.
  medians <- apply(data, 2, median)
  data <- sweep(data, 2, medians)

  ## Rotate the data and stores the rotated data frames in a list.
  rotDataList <- RotateData(data, numRot)

  ## Applies toroidal shifts to the data and stores the shifted data frames in a list.
  if (!noPropPnts) {
    lstOfRotShiftDataLists <- lapply(rotDataList, PropToroShiftData, n1, n2, propPnts)
  } else if (!noNumShifts) {
    lstOfRotShiftDataLists <- lapply(rotDataList, NumToroShiftData, n1, n2, numShifts)
  }

  ## Calculate psi for the real data
  truePsi <- mean(sapply(lstOfRotShiftDataLists, function(rotDataLst) mean(sapply(rotDataLst, psiFun, subjects))))

  # Calculate psi for all permutations of the rotated data.
  permPsi <- rep(0, numPerms)
  for (i in 1:numPerms) {
    permSubj <- sample(subjects, length(subjects), replace = FALSE)
    permPsi[i] <- mean(sapply(lstOfRotShiftDataLists, function(rotDataLst) mean(sapply(rotDataLst, psiFun, permSubj))))
  }

  list(psiStat = truePsi,
       permPsi = permPsi,
       pval = mean(c((permPsi >= truePsi), 1)))
}
