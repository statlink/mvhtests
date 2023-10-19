dcor.fbed <- function(y, x, alpha = 0.05, K = 0, backward = TRUE) {

  runtime <- proc.time()
  dm <- dim(x)
  n <- dm[1]  ;  p <- dm[2]
  la <- log(alpha)  

  r <- as.vector( dcov::mdcor(y, x, type = "U") )
  ini.stat <- n * r + 1

  ini.pvalue <- pchisq( ini.stat, 1, lower.tail = FALSE, log.p = TRUE)
  pval <- ini.pvalue
  vars <- which(pval < la)
  if ( length(vars) > 0 ) {
    sela <- which.min(pval)
  } else  sela <- vars
  vars <- setdiff(vars, sela)

  while ( length(vars) > 0 ) {
    for (vim in vars)  {
      stat <- dcov::pdcor( y, x[, vim], x[, sela], type = "U" )
      pval[vim] <- pchisq(n * stat + 1, 1, lower.tail = FALSE, log.p = TRUE)
    }
    ide <- which(pval[vars] < la)
    vars <- vars[ide]
    sel <- which.min(pval[vars])
    sela <- c(sela, vars[sel] )
    vars <- setdiff(vars, vars[sel])
  } ## end  while ( length(vars) > 0 ) {
  pvalue <- pval[sela]
  
  if ( K > 0 ) {
    for ( j in 1: K ) {
      pval <- numeric(p)
      vars <- setdiff(1:p, sela)
      
      while ( length(vars) > 0 ) {
        pval <- numeric(p)
        for (vim in vars)  {
          stat <- dcov::pdcor( y, x[, vim], x[, sela], type = "U" )
          pval[vim] <- pchisq(n * stat + 1, 1, lower.tail = FALSE, log.p = TRUE)
        }
        ide <- which(pval < la)
        vars <- vars[ide]
        if ( length(ide) > 0 ) {
          sel <- which.min(pval)
          sela <- c(sela, sel)
          pvalue <- c(pvalue, pval[sel])
          vars <- setdiff(vars, vars[sel])
        }  
      } ##  end  while ( length(vars) > 0 ) {
    } ##  end  for (j in 1:K) {
  } ##  end  if ( K > 0 ) { 
       
  runtime <- proc.time() - runtime

  len <- sum(sela > 0)
  if (len > 0) {
    res <- cbind(sela, pvalue)
  } else  res <- matrix(c(0, 0), ncol = 2)
  colnames(res) <- c("Vars", "log p-value")

  if ( backward ) {
    b <- dcor.bs(y, x[, sela])
    res <- cbind(res, b$res[, 2])
    runtime <- runtime + b$runtime
    colnames(res)[3] <- "BS log p-value"
  }

  list(runtime = runtime, res = res)
}


