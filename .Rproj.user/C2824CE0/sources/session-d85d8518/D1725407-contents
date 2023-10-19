###############################
#### Relationship between James's MANOVA and one sample Hotelling's T^2
####  8/2015
#### mtsagris@yahoo.gr
#### References: Sarah Emerson (2009)
#### Small sample performance and calibration of the Empirical Likelihood method
#### PhD thesis, pages: 76-81
################################
maovjames.hotel <- function(x, ina) {
  ## contains the data
  ## ina is a grouping  variable indicating the groups
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  nu <- tabulate(ina)
  mi <- Rfast2::colGroup(x, ina) / nu  ## group mean vectors
  p <- dim(x)[2]
  si <- array( dim = c(p, p, g) )
  for (i in 1:g)  si[, , i] <- solve( cov( x[ina == i, ] ) )
  f <- numeric(g)
  funa <- function(m) {
    for (i in 1:g)  f[i] <- nu[i] * ( mi[i, ] - m ) %*% si[, , i] %*% (mi[i, ] - m)
    sum(f)
  }
  bar <- optim(Rfast::colmeans(x), funa, control = list(maxit=2000) )
  bar <- optim(bar$par, funa, control = list(maxit=2000) )
  bar <- optim(bar$par, funa, control = list(maxit=2000) )
  list(test = bar$value, mc = bar$par)
}
