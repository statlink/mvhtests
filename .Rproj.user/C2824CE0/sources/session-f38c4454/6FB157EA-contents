equal.cov <- function(x, Sigma, a = 0.05) {
  ## x is the data set
  ## Sigma is the assumed covariance matrix
  ## a is the level of significance set by default to 0.05
  Sigma <- as.matrix(Sigma)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size
  S <- cova(x)  ## sample covariance matrix
  mesa <- solve(Sigma, S)
  test <- n * sum( diag(mesa) ) - n * log( det(mesa) ) - n * p ## test statistic
  dof <- 0.5 * p * (p + 1)  ## the degrees of freedom of  chi-square distribution
  pvalue <- pchisq(test, dof, lower.tail = FALSE)  ## p-value
  crit <- qchisq(1 - a, dof)  ## critical value of chi-square distribution
  res <- c(test, dof, pvalue, crit)
  names(res) <- c("test", "df", "p-value", "critical")
  res
}
