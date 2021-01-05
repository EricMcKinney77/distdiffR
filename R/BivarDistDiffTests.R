library(Rcpp)
# library(microbenchmark)

cppFunction('
  NumericVector bcdf(NumericMatrix data, NumericMatrix eval) {
    NumericVector x(eval.nrow());
    int i, j;
    for(i = 0; i < eval.nrow(); ++i) {
      x(i) = 0;
      for(j = 0; j < data.nrow(); ++j) {
        if( data(j, 0) <= eval(i, 0) && data(j, 1) <= eval(i, 1) ) x(i) = x(i) + 1.0;
      }
      x(i) = x(i) / data.nrow();
    }
    return x;
  }
')

CalcPsiDWS <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n1 * sum((bcdf(data1, data1) - bcdf(data2, data1))^2) +
      n2 * sum((bcdf(data2, data2) - bcdf(data1, data2))^2)) / (n1 + n2)
}

CalcPsiUWS <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  sum((bcdf(data1, data1) - bcdf(data2, data1))^2) +
    sum((bcdf(data2, data2) - bcdf(data1, data2))^2)
}

CalcPsiRWS <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n2 * sum((bcdf(data1, data1) - bcdf(data2, data1))^2) +
      n1 * sum((bcdf(data2, data2) - bcdf(data1, data2))^2)) / (n1 + n2)
}

CalcPsiDWA <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n1 * sum(abs(bcdf(data1, data1) - bcdf(data2, data1))) +
      n2 * sum(abs(bcdf(data2, data2) - bcdf(data1, data2)))) / (n1 + n2)
}

CalcPsiUWA <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  sum(abs(bcdf(data1, data1) - bcdf(data2, data1))) +
    sum(abs(bcdf(data2, data2) - bcdf(data1, data2)))
}

CalcPsiRWA <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n2 * sum(abs(bcdf(data1, data1) - bcdf(data2, data1))) +
      n1 * sum(abs(bcdf(data2, data2) - bcdf(data1, data2)))) / (n1 + n2)
}

# CalcPsi <- function(data, subjects, psiType = c("DWS", "UWS", "RWS", "DWA", "UWA", "RWA")) {
#   data1 <- data[subjects == 1, ]
#   data2 <- data[subjects == 2, ]
#   n1 <- nrow(data1)
#   n2 <- nrow(data2)
#   bcdf1 <- bcdf(data1, data1) - bcdf(data2, data1)
#   bcdf2 <- bcdf(data2, data2) - bcdf(data1, data2)
#   n1plusn2 <- n1 + n2
#   psis <- vector(mode = "numeric")
#   if (any(psiType %in% c("DWS", "UWS", "RWS"))) {
#     bcdf1sq <- sum(bcdf1^2)
#     bcdf2sq <- sum(bcdf2^2)
#   }
#   if (any(psiType %in% c("DWA", "UWA", "RWA"))) {
#     bcdf1ab <- sum(abs(bcdf1))
#     bcdf2ab <- sum(abs(bcdf2))
#   }
#   if ("DWS" %in% psiType) psis <- c(psis, (n1 * bcdf1sq + n2 * bcdf2sq) / n1plusn2)
#   if ("UWS" %in% psiType) psis <- c(psis, bcdf1sq + bcdf2sq)
#   if ("RWS" %in% psiType) psis <- c(psis, (n2 * bcdf1sq + n1 * bcdf2sq) / n1plusn2)
#   if ("DWA" %in% psiType) psis <- c(psis, (n1 * bcdf1ab + n2 * bcdf2ab) / n1plusn2)
#   if ("UWA" %in% psiType) psis <- c(psis, bcdf1ab + bcdf2ab)
#   if ("RWA" %in% psiType) psis <- c(psis, (n2 * bcdf1ab + n1 * bcdf2ab) / n1plusn2)
#   names(psis) <- psiType
#   psis
# }

RotateData <- function(data, numRotations) {
  rotDataList <- vector("list", length = numRotations)  
  rotDataList[[1]] <- data
  radRot <- 2 * pi * 0:(numRotations - 1) / numRotations
  for (i in 2:numRotations) {
    radRoti <- radRot[i] # Rotation angle in radians
    sinr <- sin(radRoti)
    cosr <- cos(radRoti)
    rotMat <- matrix(c(cosr, -sinr, sinr, cosr), 2, 2)
    rotDataList[[i]] <- data %*% rotMat
  }
  rotDataList
}

BivarDistDiffTest <- function(data1, data2, numRot = 36, numPerms = 999, 
                              psiFun = CalcPsiDWS, seedNum = NULL) {
  
  # NOTE: Data cleaning must be done before applying this function, e.g., filter(X != 0 & Y != 0) 
  
  ## Combines the data from two subjects into one long matrix.
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  subjects <- rep(1:2, times = c(n1, n2))
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
  set.seed(seedNum)
  permPsi <- rep(0, numPerms)
  for (i in 1:numPerms) {
    permSubj <- sample(subjects, length(subjects), replace = FALSE)
    permPsi[i] <- mean(sapply(rotDataList, psiFun, permSubj))
  }
  
  list(psiStat = truePsi,
       permPsi = permPsi,
       pval = mean(c((permPsi >= truePsi), 1)))
}

cppFunction('
  List ToroShiftData(NumericMatrix data, int n1, int n2, float propPnts = 1) {
    NumericMatrix::Column X = data( _, 0);
    NumericMatrix::Column Y = data( _, 1);
    float dataWidth = max(X) - min(X);
    float dataHeight = max(Y) - min(Y);
    
    int n1plusn2 = n1 + n2;
    int numShifts = round(propPnts * n1plusn2);
    List shiftDataList = List::create();
    IntegerVector n1n2vec = seq(0, n1plusn2 - 1);
    IntegerVector sampledRows = sample(n1n2vec, numShifts);
    for (int x = 0; x < numShifts; ++x) {
      int sampledRow = sampledRows[x];
      NumericVector toroCntrPnt = data(sampledRow, _);
      NumericVector shiftedX = ifelse(X < toroCntrPnt[0], X + dataWidth, X);
      NumericVector shiftedY = ifelse(Y < toroCntrPnt[1], Y + dataHeight, Y);
      shiftDataList.push_back(cbind(shiftedX, shiftedY));
    }
    return shiftDataList;
  }
')

BivarToroDiffTest <- function(data1, data2, propPnts = 1, numPerms = 999,
                              psiFun = CalcPsiDWS, seedNum = NULL) {

  # NOTE: Data cleaning must be done before applying this function, e.g., filter(X != 0 & Y != 0)

  ## Combines the data from two subjects into one long matrix.
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  subjects <- rep(1:2, times = c(n1, n2))
  data <- rbind(data1, data2)
  srtDatOrdr <- order(data[, 1], data[, 2])
  data <- data[srtDatOrdr, ]
  subjects <- subjects[srtDatOrdr]

  ## Center the data around the bivariate median of the combined data sets.
  medians <- apply(data, 2, median)
  data <- sweep(data, 2, medians)

  ## Applies toroidal shifts to the data and stores the shifted data frames in a list.
  set.seed(seedNum)
  shiftDataList <- ToroShiftData(data, n1, n2, propPnts)

  ## Calculate psi for the real data
  truePsi <- mean(sapply(shiftDataList, psiFun, subjects))

  # Calculate psi for all permutations of the rotated data.
  permPsi <- rep(0, numPerms)
  for (i in 1:numPerms) {
    permSubj <- sample(subjects, length(subjects), replace = FALSE)
    permPsi[i] <- mean(sapply(shiftDataList, psiFun, permSubj))
  }

  list(psiStat = truePsi,
       permPsi = permPsi,
       pval = mean(c((permPsi >= truePsi), 1)))
}

# library(dplyr)
# n1 = 15
# X <- rnorm(n1, mean = 0, sd = 2) %>% round(digits = 3)
# Y <- rnorm(n1, mean = 0, sd = 2) %>% round(digits = 3)
# # samp1 <- runifpoint(n1, win = owin(c(-5, 5), c(-5, 5)))
# # samp1 <- data.frame(X = samp1$x, Y = samp1$y)
# samp1 <- data.frame(X, Y)
# 
# n2 = 16
# X <- rnorm(n2, mean = 0, sd = 2) %>% round(digits = 3)
# Y <- rnorm(n2, mean = 0, sd = 2) %>% round(digits = 3)
# samp2 <- data.frame(X, Y)
# 
# samps <- rbind(samp1, samp2)
# samps <- cbind(Subject = c(rep(1, n1), rep(2, n2)), samps)
# rm(X, Y, n1, n2)
# 
# samp1 <- as.matrix(samp1)
# samp2 <- as.matrix(samp2)

RotateToroDiffTest <- function(data1, data2, numRot = 36, propPnts = 0.9,
                               numPerms = 999, psiFun = CalcPsiDWS, seedNum = NULL) {
  
  # NOTE: Data cleaning must be done before applying this function, e.g., filter(X != 0 & Y != 0) 
  
  ## Combines the data from two subjects into one long matrix.
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  subjects <- rep(1:2, times = c(n1, n2))
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
  set.seed(seedNum)
  lstOfRotShiftDataLists <- lapply(rotDataList, ToroShiftData, n1, n2, propPnts)
  
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

# ## Sample data for development and testing purposes
# set.seed(84321)
# 
# n1 = 500
# X <- round(rnorm(n1, mean = 3, sd = 2), digits = 3)
# Y <- round(rnorm(n1, mean = 3, sd = 2), digits = 3)
# data1 <- as.matrix(data.frame(X, Y))
# 
# n2 = 500
# X <- round(rnorm(n2, mean = 0, sd = 2), digits = 3)
# Y <- round(rnorm(n2, mean = 0, sd = 2), digits = 3)
# data2 <- as.matrix(data.frame(X, Y))
# 
# 
# junk2 <- BivarDistDiffTest(data1, data2, numRot = 36, numPerms = 99,
#                              psiFun = CalcPsiDWS, seedNum = 123)
# junk2
# 
# microbenchmark(
#   junk2 <- BivarDistDiffTest(data1, data2, numRot = 36, numPerms = 99,
#                              psiFun = CalcPsiDWS, seedNum = 123),
#   times = 5) # median 27.6 sec on CHPC
# 
# junk3 <- BivarToroDiffTest(data1, data2, propPnts = 0.1, numPerms = 99,
#                            psiFun = CalcPsiDWS, seedNum = 123)
# junk3
# 
# microbenchmark(
#   junk3 <- BivarToroDiffTest(data1, data2, propPnts = 0.1, numPerms = 99,
#                              psiFun = CalcPsiDWS, seedNum = 123),
#   times = 5) # median 66.9 sec on CHPC
# 

# diff1 <- function(X) {
#   max(X) - min(X)
# }
# 
# diff2 <- function(X) {
#   diff(range(X))
# }
# 
# X <- runif(1000)
# 
# diff1(X)
# diff2(X)
# 
# microbenchmark(
#   junk1 <- diff1(X),
#   junk2 <- diff2(X),
#   times = 100)
