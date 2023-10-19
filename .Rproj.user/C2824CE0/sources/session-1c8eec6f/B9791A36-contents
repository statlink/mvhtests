dcor.bs <- function(y, x, alpha = 0.05) {
  
  runtime <- proc.time()
  dm <- dim(x)
  n <- dm[1]  ;  p <- dm[2]
  la <- log(alpha)  
  poies <- 1:p  
  pval <- rep(-Inf, p)

  for ( vim in 1:p ) {
    stat <- dcov::pdcor( y, x[, vim], x[, -vim], type = "U" )
    pval[vim] <- pchisq(n * stat + 1, 1, lower.tail = FALSE, log.p = TRUE)
  }
  pvalue <- pval

  sel <- which.max(pval)
  ide <- pval[sel] > la

  if ( ide ) {

    sela <- sel
    pval <- rep(-Inf, p)
  
    while ( ide  &  length( poies[-sela] ) > 1 ) {
      for ( vim in poies[-sela] ) {
        stat <- dcov::pdcor( y, x[, vim], x[, -c(vim, sela)], type = "U" )
        pval[vim] <- pchisq(n * stat + 1, 1, lower.tail = FALSE, log.p = TRUE)
      }
      sel <- which.max(pval)
      ide <- pval[sel] > la
      if (ide) {
        sela <- c(sela, sel)
        pval[sela] <-  -Inf
      }
    } ##  end  while ( ide  &  length( poies[-sela] ) > 1 ) {

    if ( ide  &  length( poies[-sela] ) == 1 ) {
      vim <- poies[-sela]
      stat <- dcov::dcor( y, x[, vim], type = "U" )
      pval[vim] <- pchisq(n * stat + 1, 1, lower.tail = FALSE, log.p = TRUE)
      ide <- pval[vim] > la
      if (ide) {
        sela <- c(sela, vim)
        pval[vim] <-  -Inf
      }
    } ##  end  if ( ide  &  length( poies[-sela] ) == 1 ) {

    if ( length(sela) > 0 )  pvalue[sela] <- pval[sela]

  } ##  end  if ( ide ) {

  runtime <- proc.time() - runtime
   
  len <- which( is.infinite(pvalue) )
  if ( length(len) > 0 )  pvalue[len] <- 0
  res <- cbind(1:p, pvalue)
  colnames(res) <- c("Vars", "log p-value")

  list(runtime = runtime, res = res)
}

  

   
    





