#' Create rotated versions of the data
#'
#' This function produces a list of rotated versions of the two-column input
#' matrix. Specifically, the number of rotations (i.e., an integer passed to
#' `numRotations`) divides a complete circle into `numRotations` equal angles,
#' and `numRotations` rotated versions of the input data are output in a list.
#'
#' @param data A two column matrix of the bivariate combined samples
#' @param numRotations A non-negative integer specifying the number of rotations to be applied to the data within 360 degrees.
#'
#' @return A list of matrices containing the coordinates for each version of the rotated data (including the original data, which is the first matrix of the list)
#' @export
#' @examples
#' data(iris)
#' sample1 <- as.matrix(iris[iris$Species == "setosa", 1:2])
#'
#' # Generate five rotated versions of sample1 (every 72 degrees) within 360 degrees.
#' RotateData(sample1, 5)
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
