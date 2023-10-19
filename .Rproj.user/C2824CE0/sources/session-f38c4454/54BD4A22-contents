dcor.bsmmpc <- function(y, x, max_k = 3, alpha = 0.05, B = 999) {

  runtime <- proc.time()

  p <- dim(x)[2]  
  pval <- numeric(p)
  vars <- 1:p
  
  if ( p > 1 ) {
    for ( i in 1:min( max_k, p) ) {
      for (vim in 1:p) {
        pval2 <- NULL
        cand <- Rfast::comb_n(vars[-vim], i)
        for ( j in 1:dim(cand)[2] ) {  
          pval2 <- c(pval2, dcov::pdcor.test(y, x[, vim], x[, cand[, j]], R = B, type = "U" )$p.values )
        }
        pval[vim] <- max(pval[vim], pval2)
      } ##  end  for (vim in 1:p) {   
    } ##  end  for ( i in 1:min( max_k, p) ) {
  } ##  end  if ( p > 1 ) {      

  runtime <- proc.time() - runtime
  res <- cbind(vars, pval)
  colnames(res) <- c("Vars", "p-value") 
  list(runtime = runtime, res = res)
}
