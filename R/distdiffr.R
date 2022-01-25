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

  # Check for conflict between propPnts and numShifts.
  usePropPnts <- !is.null(propPnts)
  useNumShifts <- !is.null(numShifts)
  n_pooled <- n1 + n2
  if (usePropPnts & useNumShifts) {
    stop("Must provide either propPnts or numShifts, but not both.")
  } else if (!usePropPnts & useNumShifts) { # Use numShifts instead of shiftThrshld.
    warning("Using numShifts instead of default shiftThrshld.")
    shiftThrshld <- NULL
    if (n_pooled < numShifts) {
      stop("n_pooled is smaller than numShifts.")
    }
  } else if (usePropPnts & !useNumShifts) { # Use propPnts and instead of shiftThrshld.
    print("Using propPnts instead of default shiftThrshld.")
    shiftThrshld <- NULL
  } else if (!usePropPnts & !useNumShifts) { # Use default shiftThrshld.
    if (n_pooled < shiftThrshld) {
      warning("n_pooled is smaller than shiftThrshld.\nCan only compute n_pooled toroidal shifts.")
      numShifts <- n_pooled
    } else {
      numShifts <- shiftThrshld
    }
    useNumShifts <- TRUE
  }

  hash1 <- hashMat(data1)
  hash2 <- hashMat(data2)

  if (hash1 >= hash2) {
    data <- rbind(data1, data2)
    subjects <- rep(1:2, times = c(n1, n2))
  } else {
    data <- rbind(data2, data1)
    subjects <- rep(2:1, times = c(n2, n1))
  }

  set.seed(seedNum)

  ## Center the data around the bivariate median of the combined data sets.
  medians <- apply(data, 2, median)
  data <- sweep(data, 2, medians)

  ## Rotate the data and stores the rotated data frames in a list.
  if (testType != "toroidal") {
    rotDataList <- RotateData(data, numRot)

    if (testType == "combined") {
      ## Applies toroidal shifts to the data and stores the shifted data frames in a list.
      if (usePropPnts) {
        lstOfRotShiftDataLists <- lapply(rotDataList, PropToroShiftData, n1, n2, propPnts)
      } else if (useNumShifts) {
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
    if (usePropPnts) {
      shiftDataList <- PropToroShiftData(data, n1, n2, propPnts)
    } else if (useNumShifts) {
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
