################################
#### Relationship between James and one sample Hotelling's T^2
####  8/2015
#### mtsagris@yahoo.gr
#### References: Sarah Emerson (2009)
#### Small sample performance and calibration of the Empirical Likelihood method
#### PhD thesis, pages: 76-81
################################
james.hotel <- function(x1, x2) {
  ## x and y are the multivariate samples
  n1 <- dim(x1)[1]  ;   n2 <- dim(x2)[1]  ## sample sizes
  m1 <- Rfast::colmeans(x1)  ## sample mean vector of the first sample
  m2 <- Rfast::colmeans(x2)  ## sample mean vector of the second sample
  dbar <- m2 - m1  ## difference of the two mean vectors
  s1 <- Rfast::cova(x1)      ;    s2 <- Rfast::cova(x2)
  A1 <- s1/n1   ;    A2 <- s2/n2
  V <- A1 + A2  ## covariance matrix of the difference
  test <- as.numeric(   dbar %*% solve(V, dbar) )
  a1inv <- chol2inv( chol(A1) )
  a2inv <- chol2inv( chol(A2) )
  mc <- solve( a1inv + a2inv ) %*% ( a1inv %*% m1 + a2inv %*% m2 )
  a1 <- a1inv / n1   ;   a2 <- a2inv / n2
  funa <- function(m) {
    n1 * (m1 - m) %*% a1 %*% ( m1 - m ) +  ## Hotelling's test statistic
    n2 * (m2 - m) %*% a2 %*% ( m2 - m )
  }
  bar <- optim( (m1 + m2)/2, funa )
  bar <- optim( bar$par, funa )
  tests <- c(test, bar$value )
  names(tests) <- c("James test", "T1(mu) + T2(mu)")
  list(tests = tests, mathematics.mean = t(mc), optimised.mean = bar$par )
}
