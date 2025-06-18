
# https://internal-notes.hakyimlab.org/post/2022-10-25-how-to-calculate-proportion-of-variance-explained-by-variables/

calc_loglik_fast = function(h2, lambdavec, tyyCC2)
{
  # ## Need to pre-calculate the following variables that don't depend on h2
  # ## leverage the fact that KK has the same eigenvectors as XRM
  # tempo <- eigen(XRM)
  # lambdavec <- tempo$values
  # CC <- tempo$vectors
  # tyyCC <- t(yy) %*% CC
  # tyyCC2 <- tyyCC * tyyCC  

  ## this is what I want:   res = -0.5 * log(det(KK)) -0.5 * t(yy) %*% iKK %*% yy 
  lamplusvec = h2 * lambdavec + 1 - h2
  logdetKK = sum( log( lamplusvec ) )
    res = -0.5 * logdetKK - 0.5 * sum( tyyCC2 /lamplusvec)
    res
}

## approximate second derivative
d2h2_loglik = function(h2,lambdavec,tyyCC2,dh2=0.01)
{

  myfun = function(hh) calc_loglik_fast(hh, lambdavec, tyyCC2)
  (myfun(h2 + dh2) - 2*myfun(h2) + myfun(h2 - dh2) ) / (4*dh2^2)
}

sd_h2_mle <- function(h2,lambdavec,tyyCC2,...)
{
  sqrt( -1/d2h2_loglik(h2,lambdavec,tyyCC2,...) )
}

##The Fisher information is the expected value of the observed information given a single observation X. X distributed according to the hypothetical model with parameter \theta : https://en.wikipedia.org/wiki/Observed_information

## calc mle using optimize function
calc_mle = function(yy,XX,plotit=FALSE,rango=c(0.01,0.95))
{
  ## yy should be nsam x 1 vector (nsam is number of samples)
  ## XX should be nsam x mp matrix (mp is number of features)
  ## scale both yy and XX
  yy = scale(yy)
  XX = scale(XX)

  mp <- ncol(XX)
  nsam <- length(c(yy))
  
  XRM <- XX %*% t(XX) / mp
  tempo <- eigen(XRM)
  lambdavec <- tempo$values
  CC <- tempo$vectors
  tyyCC <- t(yy) %*% CC
  tyyCC2 <- tyyCC * tyyCC  
  
  myfun = function(hh) calc_loglik_fast(hh, lambdavec, tyyCC2)
  res <- optimize(myfun, c(0.01,0.99),  maximum=TRUE )
  res$se = sd_h2_mle(res$maximum,lambdavec,tyyCC2)
  res$pval = pnorm(-abs(res$maximum/res$se))*2
  if(plotit)
  {
    mm=20
    lvec = rep(NA,mm)
    hvec = seq(from=rango[1],to=rango[2],length.out=mm)
    for(cc in 1:mm)
      lvec[cc] = calc_loglik_fast(hvec[cc], lambdavec, tyyCC2)
    plot(hvec,lvec,xlab="h2 (=PVE)")
    title(main=glue("h2 mle = {signif(res$maximum,2)}; se = {signif(res$se,2)}; p = {signif(res$pval,2)}"))
    abline(v=c(res$maximum,res$maximum - 2*res$se, res$maximu + 2*res$se), lty='dotted')
  }
  res
}


# I am modifiying this because I don;'t want to start directly from the SNPs


calc_mle_from_grm = function(yy,grmXX,mp,plotit=FALSE,rango=c(0.01,0.95))
{
  ## yy should be nsam x 1 vector (nsam is number of samples)
  ## XX should be nsam x mp matrix (mp is number of features)
  ## scale both yy and XX
#   yy = scale(yy)
#   XX = scale(XX)

  #mp <- ncol(XX) # 
  nsam <- length(c(yy))
  XRM <- grmXX
  tempo <- eigen(XRM)
  lambdavec <- tempo$values
  CC <- tempo$vectors
  tyyCC <- t(yy) %*% CC
  tyyCC2 <- tyyCC * tyyCC  
  
  myfun = function(hh) calc_loglik_fast(hh, lambdavec, tyyCC2)
  res <- optimize(myfun, c(0.01,0.99),  maximum=TRUE )
  res$se = sd_h2_mle(res$maximum,lambdavec,tyyCC2)
  res$pval = pnorm(-abs(res$maximum/res$se))*2
  if(plotit)
  {
    mm=20
    lvec = rep(NA,mm)
    hvec = seq(from=rango[1],to=rango[2],length.out=mm)
    for(cc in 1:mm)
      lvec[cc] = calc_loglik_fast(hvec[cc], lambdavec, tyyCC2)
    # {
    #   plot(hvec,lvec,xlab="h2 (=PVE)")
    #   title(main=glue("h2 mle = {signif(res$maximum,2)}; se = {signif(res$se,2)}; p = {signif(res$pval,2)}"))
    #   abline(v=c(res$maximum,res$maximum - 2*res$se, res$maximu + 2*res$se), lty='dotted')
    # }

    lp <- ggplot2::ggplot() +
      geom_line(aes(x = hvec, y = lvec), color = 'grey') +
      geom_point(aes(x = hvec, y = lvec), color = 'red', pch = 21) +
      labs(title = glue("h2 mle = {signif(res$maximum,2)}; se = {signif(res$se,2)}; p = {signif(res$pval,2)}")) +
      xlab(expression(h[est]^2)) + ylab('Log Likelihood') +
      geom_vline(aes(xintercept = c(res$maximum)), linetype = 'dashed', color = 'grey') +
      theme_minimal()
      # geom_vline(aes(xintercept = c(res$maximum, res$maximum - 2*res$se, res$maximu + 2*res$se)), linetype = 'dotted')

    likelihood_plots <<- append(likelihood_plots, list(lp))
  }
  res
}