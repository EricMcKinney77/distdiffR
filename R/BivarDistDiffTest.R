#' The rotational distDiffR test
#'
#' @param data1 a two column matrix of bivariate observations from one sample
#' @param data2 a two column matrix of bivariate observations from another sample
#' @param numRot an integer number of rotational shifts of the pooled samples
#' @param numPerms an integer number of permutations of the original data
#' @param psiFun a function specifying the Psi statistic calculation
#' @param seedNum an integer random seed value
#'
#' @return
#' @importFrom stats runif
#' @importFrom stats median
#' @export
BivarDistDiffTest <- function(data1, data2, numRot = 8, numPerms = 999,
                              psiFun = CalcPsiRWS, seedNum = NULL) {

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

  ## Calculate psi for the real data
  truePsi <- mean(sapply(rotDataList, psiFun, subjects))

  # Calculate psi for all permutations of the rotated data.
  permPsi <- rep(0, numPerms)
  for (i in 1:numPerms) {
    permSubj <- sample(subjects, length(subjects), replace = FALSE)
    permPsi[i] <- mean(sapply(rotDataList, psiFun, permSubj))
  }

  list(psiStat = truePsi,
       permPsi = permPsi,
       pval = mean(c((permPsi >= truePsi), 1)))
}
