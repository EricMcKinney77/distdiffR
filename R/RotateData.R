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
