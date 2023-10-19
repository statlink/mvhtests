sarabai <- function(x1, x2) {
  ## x and y are high dimensional datasets
  n1 <- dim(x1)[1]    ;    n2 <- dim(x2)[1]  ## sample sizes
  m1 <- Rfast::colmeans(x1)   ;   m2 <- Rfast::colmeans(x2)  ## sample means
  n <- n1 + n2 - 2
  z1 <- t(x1) - m1
  z2 <- t(x2) - m2
  Sn <- tcrossprod( z1 ) + tcrossprod( z2 )
  ## Sn is the pooled covariance matrix
  trSn <- sum( z1^2 )/n + sum( z2^2 ) /n
  trSn2 <- sum(Sn^2)/n^2
  Bn <- sqrt( n^2/( (n + 2) * (n - 1) ) * (trSn2 - trSn^2/n) )
  up <- n1 * n2/(n1 + n2) * sum( (m1 - m2)^2 ) - trSn
  down <- sqrt( 2 * (n + 1)/n ) * Bn
  Z <- up / down  ## test statistic
  pvalue <- pnorm(Z, lower.tail = FALSE)
  res <- c(Z, pvalue)
  names(res) <- c('Z', 'p-value')
  res
}
