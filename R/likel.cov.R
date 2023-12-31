likel.cov <- function(x, ina, a = 0.05) {
  ## x is the data set
  ## ina is a numeric vector indicating the groups of the data set
  ## a is the level of significance, set to 0.05 by default
  ina <- as.numeric(ina)
  p <- dim(x)[2]  ## dimension of the data set
  n <- dim(x)[1]  ## total sample size
  k <- max(ina)  ## number of groups
  nu <- tabulate(ina) ## the sample size of each group
  t1 <- rep( (nu - 1)/nu, each = p^2 )
  t2 <- rep(nu - 1, each = p^2 )
  s <- array( dim = c(p, p, k) )
  ## the next 3 lines create the pooled covariance matrix
  ## and calculate the covariance matrix of each group
  for (i in 1:k)  s[, , i] <- Rfast::cova( x[ina == i, ] )
  mat <- t1 * s
  mat1 <- t2 * s
  Sp <- colSums( aperm(mat1) ) / n
  deta <- apply(mat, 3, det)
  pame <- det(Sp) / deta
  test <- sum(nu * log(pame))  ## test statistic
  dof <- 0.5 * p * (p + 1) * (k - 1) ## degrees of freedom of chi-square
  pvalue <- pchisq(test, dof, lower.tail = FALSE)  ## p-value of test statistic
  crit <- qchisq(1 - a, dof)  ## critical value of chi-square distribution
  res <- c(test, pvalue, dof, crit)
  names(res) <- c('test', 'p-value', 'df', 'critical')
  res
}
