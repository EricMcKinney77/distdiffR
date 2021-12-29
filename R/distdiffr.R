#' The distdiffr test
#'
#' @param data1 A two column matrix of bivariate observations from one sample.
#' @param data2 A two column matrix of bivariate observations from another sample.
#' @param testType A string indicating the type of test to be used. Must be one of c("rotational", "toroidal", "combined").
#' @param numRot An integer number of rotational shifts of the pooled samples.
#' @param propPnts A numeric proportion of points to be used as toroidal shift origins. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.
#' @param numShifts A numeric integer. The number of points to be used as toroidal shift origins. Must be less than the pooled sample size. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.
#' @param shiftThrshld A numeric integer. Used if neither propPnts or numShifts are provided. If the pooled sample size is less than shiftThrshld, every point will be used as a toroidal shift origin. Otherwise, only a random sample of shiftThrshld points will be used.
#' @param numPerms An integer number of permutations of the original data.
#' @param psiFun A function specifying the Psi statistic calculation.
#' @param seedNum An integer random seed value.
#'
#' @return A list including three objects:
#'     (1) the Psi statistic computed on the original data
#'     (2) a vector of Psi statistics computed on the permuted data
#'     (3) the p-value for the test
#' @importFrom stats runif
#' @importFrom stats median
#' @export
distdiffr <- function(data1,
                      data2,
                      testType = "combined",
                      numRot = 8,
                      propPnts = NULL,
                      numShifts = NULL,
                      shiftThrshld = 25,
                      numPerms = 999,
                      psiFun = CalcPsiCWS,
                      seedNum = NULL) {
  ## Combines the data from two subjects into one long matrix.
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  subjects <- rep(1:2, times = c(n1, n2))

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

  # Check and remove duplicate data values between samples
  data1DuplRowsWdata2 <- function(data1, data2) {
    data1strgs <- unlist(sapply(1:n1, function(j) paste(data1[j, ], collapse = '_')))
    data2strgs <- unlist(sapply(1:n2, function(j) paste(data2[j, ], collapse = '_')))
    which(data1strgs %in% intersect(data1strgs, data2strgs))
  }
  duplRowsInd <- data1DuplRowsWdata2(data1, data2)

  if (length(duplRowsInd) != 0) {
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
    while (length(duplRowsInd) != 0) {
      data1[duplRowsInd, ] <- data1[duplRowsInd, ] + runif(length(duplRowsInd) * 2, -roundBound, roundBound)
      duplRowsInd <- data1DuplRowsWdata2(data1, data2)
    }
  }

  set.seed(seedNum)

  data <- rbind(data1, data2)
  srtDatOrdr <- order(data[, 1], data[, 2])
  data <- data[srtDatOrdr, ]
  subjects <- subjects[srtDatOrdr]

  ## Center the data around the bivariate median of the combined data sets.
  medians <- apply(data, 2, median)
  data <- sweep(data, 2, medians)

  ## Rotate the data and stores the rotated data frames in a list.
  if (testType != "toroidal") {
    rotDataList <- RotateData(data, numRot)

    if (testType == "combined") {
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

    } else if (testType == "rotational") {
      ## Calculate psi for the real data
      truePsi <- mean(sapply(rotDataList, psiFun, subjects))

      # Calculate psi for all permutations of the rotated data.
      permPsi <- rep(0, numPerms)
      for (i in 1:numPerms) {
        permSubj <- sample(subjects, length(subjects), replace = FALSE)
        permPsi[i] <- mean(sapply(rotDataList, psiFun, permSubj))
      }
    }
  } else if (testType == "toroidal") {
    if (!noPropPnts) {
      shiftDataList <- PropToroShiftData(data, n1, n2, propPnts)
    } else if (!noNumShifts) {
      shiftDataList <- NumToroShiftData(data, n1, n2, numShifts)
    }

    ## Calculate psi for the real data
    truePsi <- mean(sapply(shiftDataList, psiFun, subjects))

    # Calculate psi for all permutations of the rotated data.
    permPsi <- rep(0, numPerms)
    for (i in 1:numPerms) {
      permSubj <- sample(subjects, length(subjects), replace = FALSE)
      permPsi[i] <- mean(sapply(shiftDataList, psiFun, permSubj))
    }
  }

  list(psiStat = truePsi,
       permPsi = permPsi,
       pval = mean(c((permPsi >= truePsi), 1)))
}
