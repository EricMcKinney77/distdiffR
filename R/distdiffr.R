#' The distdiffR two-sample tests of bivariate distributional equality
#'
#' The `distdiffr()` function conducts two-sample permutation tests of
#' distributional equality based on differences in the bivariate empirical
#' cumulative density functions (BECDFs). The differences in BECDFs are computed
#' across a series of rotations, toroidal shifts, or both rotations and toroidal
#' shifts of the combined data (specified via `testType`). The number of rotations
#' and toroidal shifts may be specified (via `numRot` or `numShifts`,
#' respectively). The number of toroidal shifts may also be determined by a
#' proportion of the combined sample size (via `propPnts`). However,
#' \insertCite{mckinney2022extensions;textual}{distdiffR} has shown that
#' limiting the number of toroidal shifts to ease the computational load of the
#' test will still provide stable results. Simulations have shown the combined
#' rotational and toroidal shift test to be the most powerful yet appropriately
#' conservative test. For more information, see
#' \insertCite{mckinney2022extensions;textual}{distdiffR} and
#' \insertCite{mckinney2021extensions;textual}{distdiffR}.
#'
#' @param data1 A two column matrix of bivariate observations from one sample.
#' @param data2 A two column matrix of bivariate observations from another sample.
#' @param testType A string indicating the type of test to be used. Must be one of c("rotational", "toroidal", "combined").
#' @param numRot An integer number of rotational shifts of the pooled samples.
#' @param propPnts A numeric proportion of points to be used as toroidal shift origins. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.
#' @param numShifts A numeric integer. The number of points to be used as toroidal shift origins. Must be less than the pooled sample size. Cannot provide both propPnts and numShifts. If neither are provided, shiftThrshld is used.
#' @param shiftThrshld A numeric integer. Used if neither propPnts or numShifts are provided. If the pooled sample size is less than shiftThrshld, every point will be used as a toroidal shift origin. Otherwise, only a random sample of shiftThrshld points will be used.
#' @param numPerms An integer number of permutations of the original data.
#' @param psiStat A string specifying the Psi statistic calculation. Must be one of c("CWA", "DWA", "UWA", "CWS", "DWS", "UWS").
#' @param seedNum An integer random seed value.
#'
#' @return A list including three objects:
#'     (1) the Psi statistic computed on the original data
#'     (2) a vector of Psi statistics computed on the permuted data
#'     (3) the p-value for the test
#' @importFrom stats runif
#' @importFrom stats median
#' @importFrom Rdpack reprompt
#' @export
#' @references
#' \insertRef{mckinney2022extensions}{distdiffR}
#'
#' \insertRef{mckinney2021extensions}{distdiffR}
#' @examples
#' # Randomly assign all three species to two samples
#' seedNum <- 123
#' set.seed(seedNum)
#'
#' data(iris)
#' # Randomly assign all three species to two samples
#' irisPermuted <- iris[sample.int(nrow(iris)), ]
#' sample1 <- as.matrix(irisPermuted[1:75, 1:2])
#' sample2 <- as.matrix(irisPermuted[76:150, 1:2])
#' pooled_data <- rbind(cbind(sample1, 1), cbind(sample2, 2))
#'
#' # Rotational test
#' output <- distdiffr(sample1, # Note: Data inputs must be matrices
#'   sample2,
#'   testType = "rotational",
#'   numRot = 8, # Default value
#'   seedNum = seedNum
#' )
#' output$pval
#'
#' # Toroidal shift test with proportions of points
#' output <- distdiffr(sample1,
#'   sample2,
#'   testType = "toroidal",
#'   propPnts = 0.1,
#'   seedNum = seedNum
#' )
#' output$pval
#'
#' # Toroidal shift test with a threshold below pooled sample size
#' output <- distdiffr(sample1,
#'   sample2,
#'   testType = "toroidal",
#'   shiftThrshld = 25, # Default
#'   seedNum = seedNum
#' )
#' output$pval
#'
#' # Toroidal shift test with a threshold above pooled sample size
#' output <- distdiffr(sample1,
#'   sample2,
#'   testType = "toroidal",
#'   shiftThrshld = 200,
#'   seedNum = seedNum
#' )
#' output$pval
#'
#' # Toroidal shift test with a number of shifts
#' output <- distdiffr(sample1,
#'   sample2,
#'   testType = "toroidal",
#'   numShifts = 8,
#'   seedNum = seedNum
#' )
#' output$pval
#'
#' # Combined rotational and toroidal shift test
#' output <- distdiffr(sample1,
#'   sample2,
#'   testType = "combined", # Default
#'   numRot = 8, # Default
#'   shiftThrshld = 25, # Default
#'   seedNum = seedNum
#' )
#' output$pval
#'
#' # Also see browseVignettes(package = "distdiffR")
distdiffr <- function(data1,
                      data2,
                      testType = "combined",
                      numRot = 8,
                      propPnts = NULL,
                      numShifts = NULL,
                      shiftThrshld = 25,
                      numPerms = 999,
                      psiStat = "CWS",
                      seedNum = NULL) {
  n1 <- nrow(data1)
  n2 <- nrow(data2)

  # Determine numShifts
  usePropPnts <- !is.null(propPnts)
  useNumShifts <- !is.null(numShifts)
  n_pooled <- n1 + n2

  if (usePropPnts && useNumShifts) {
    stop("Must provide either propPnts or numShifts, but not both.")
  }

  actualShifts <- shiftThrshld
  if (usePropPnts) {
    actualShifts <- round(propPnts * n_pooled)
  } else if (useNumShifts) {
    if (numShifts >= n_pooled) {
      stop("number of shifts larger than the combined sample sizes!")
    }
    actualShifts <- numShifts
  } else {
    actualShifts <- ifelse(n_pooled < shiftThrshld, n_pooled, shiftThrshld)
  }

  # Handle testType
  finalNumRot <- if (testType == "toroidal") 1 else numRot
  finalNumShifts <- if (testType == "rotational") 1 else actualShifts

  # Map psiStat to engine type
  psiMap <- c(
    "CWA" = 0, "DWA" = 1, "UWA" = 2,
    "CWS" = 3, "DWS" = 4, "UWS" = 5
  )
  statType <- psiMap[psiStat]
  if (is.na(statType)) statType <- 3 # Default to CWS

  # Consistency hash for data order
  hash1 <- hashMat(data1)
  hash2 <- hashMat(data2)
  if (hash1 >= hash2) {
    data <- rbind(data1, data2)
    subjects <- rep(1:2, times = c(n1, n2))
  } else {
    data <- rbind(data2, data1)
    subjects <- rep(2:1, times = c(n2, n1))
  }

  # Center data around bivariate median
  medians <- apply(data, 2, median)
  data <- sweep(data, 2, medians)

  if (is.null(seedNum)) seedNum <- 42

  res <- distdiffR_engine(data, subjects, finalNumRot, finalNumShifts, statType, numPerms, seedNum)

  list(
    psiStat = res$psiStat,
    permPsi = res$permPsi,
    pval = mean(c((res$permPsi >= res$psiStat), 1))
  )
}
