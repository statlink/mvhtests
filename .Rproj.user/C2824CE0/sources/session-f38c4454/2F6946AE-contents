dcor.mmpc <- function(y, x, max_k = 3, alpha = 0.05, B = 999, backward = TRUE) {

  runtime <- proc.time()

  p <- dim(x)[2]
  ini.stat <- ini.pvalue <- numeric(p)

  for ( vim in 1:p )  {
    a <- dcov::dcor.test(y, x[, vim],R = B, type = "U")
    ini.stat[vim] <- a$statistic
    ini.pvalue[vim] <- a$p.values
  }

  pval <- ini.pvalue
  vars <- which(pval < alpha)
  if ( length(vars) > 0 ) {
    sela <- which.min(pval)
  } else  sela <- vars
  vars <- setdiff(vars, sela)

  while ( length(vars) > 0 ) {
    pval2 <- numeric(p)
    for ( i in 1:min( max_k, length(sela) ) ) {
      if ( length(sela) == 1 ) {
        cand <- matrix(sela, nrow = 1)
      } else  cand <- Rfast::comb_n(sort(sela), i)
      j <- 1
      while ( length(vars) > 0  &  j <= dim(cand)[2] ) {
        for (vim in vars)   pval2[vim] <- dcov::pdcor.test(y, x[, vim], x[, cand[, j]], R = B, type = "U" )$p.values
        pval[vars] <- pmax(pval[vars], pval2[vars])
        ide <- which(pval[vars] < alpha)
        vars <- vars[ide]
        j <- j + 1
      }  ## end  while ( length(vars) > 0  &  j <= dim(cand)[2] ) {
    }  ## end  for ( i in 1:min( max_k, length(sela) ) ) {
    sel <- which.min(pval[vars])
    sela <- c(sela, vars[sel] )
    vars <- setdiff(vars, vars[sel])
  } ## end  while ( length(vars) > 0 ) {

  pvalue <- pval

  runtime <- proc.time() - runtime

  sela <- which(pvalue < alpha)
  len <- sum(sela > 0)
  if (len > 0) {
    res <- cbind(sela, pvalue[sela])
  } else  res <- matrix(c(0, 0), ncol = 2)
  colnames(res) <- c("Vars", "p-value") 
  
  if ( backward ) {
    b <- dcor.bsmmpc(y, x[, sela])
    res <- cbind(res, b$res[, 2])
    runtime <- runtime + b$runtime
    colnames(res)[3] <- "BS p-value"
  }
  
  list(runtime = runtime, pvalue = pvalue, res = res)
}


