CalcPsiRWS <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n2 * sum((bcdf(data1, data1) - bcdf(data2, data1))^2) +
      n1 * sum((bcdf(data2, data2) - bcdf(data1, data2))^2)) / (n1 + n2)
}
