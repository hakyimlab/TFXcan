




qq_generic <- function(data_points, distribution=c('normal', 'chisq', 'uniform', 't', 'poisson', 'neg_binom', 'gamma'), data_type='p_values', BH=F, CI=F, mlog10_p_thres=30,
                       neg_log10_values=F, mean=0, sd=1, df=1, min=0, max=1, lambda=0.5, shape=20, rate=1, params_to_plot=NULL, params_to_abline=NULL, params_to_legend=NULL){
    
    n <- length(data_points)
    data_points <- na.omit(data_points)
    
    if(length(distribution) > 1){
        stop('Input one distribution.')
    }
    
    p_thres = 10^{-mlog10_p_thres}
    if( sum( data_points < p_thres) ){
      warning(paste("thresholding p to ", p_thres) )
      data_points = pmax(data_points, p_thres)
    }

    theoretical_quantiles <- switch(as.character(distribution), 
                                    'normal' = qnorm(ppoints(n), mean = mean, sd = sd)[order(order(data_points))],
                                    'chisq' = qchisq(ppoints(n), df=df)[order(order(data_points))],
                                    'uniform' = qunif(ppoints(n), min=min, max=max)[order(order(data_points))],
                                    't' = qt(ppoints(n), df=df)[order(order(data_points))],
                                    'poisson' = qpois(ppoints(n), lambda = lambda)[order(order(data_points))],
                                    'gamma' = qgamma(ppoints(n), shape = shape)[order(order(data_points))],
                                    #'neg_binom' = qnbinom(ppoints(n), )
                                    stop('The distribution is invalid. Input a valid distribution.')
    )
    
    xlab_use <- 'Theoretical Quantiles'
    ylab_use <- 'Sample Quantiles'
    
    if(neg_log10_values == T){
        theoretical_quantiles <- -log10(theoretical_quantiles)
        data_points <- -log10(data_points) 

        if(data_type == 'p_values'){
            xlab_use <- bquote('Expected -'~log[10](italic(p)))
            ylab_use <- bquote('Observed -'~log[10](italic(p)))
        } else {
            xlab_use <- '-log10(Theoretical quantiles)'
            ylab_use <- '-log10(Sample quantiles)'
        }
    }
    
    do.call(base::plot, c(list(x=theoretical_quantiles, y=data_points, xlab = xlab_use, ylab = ylab_use), params_to_plot))
    #plot(theoretical_quantiles, data_points, xlab = xlab_use, ylab = ylab_use, ...)
    #do.call(graphics::abline, c(list(a = 0, b = 1, col = "red"), params_to_abline))
    if(BH == TRUE){
        do.call(graphics::abline, c(list(a=-log10(0.05),b=1, col='orange',lty=1)))
        do.call(graphics::abline, c(list(a=-log10(0.10),b=1, col='green',lty=2)))
        do.call(graphics::abline, c(list(a=-log10(0.25),b=1, col='blue',lty=3)))
    }
    abline(h=-log10(0.05/n), col='red')
    do.call(graphics::legend, c(list(legend=c("FDR = 0.05","FDR = 0.10","FDR = 0.25"), col=c('orange','green','blue'), lty=1:3, cex=1), params_to_legend))
}


qqR2 <- function(corvec,nn,pad_neg_with_0 = FALSE,...)
{
## nn is the sample size, number of individuals used to compute correlation.
## needs correlation vector as input.
## nullcorvec generates a random sample from correlation distributions, under the null hypothesis of 0 correlation using Fisher's approximation.
  if(pad_neg_with_0) corvec[corvec < 0 | is.na(corvec) ]=0
  mm <- length(corvec)
  nullcorvec = tanh(rnorm(mm)/sqrt(nn-3)) ## null correlation vector
  qqplot(nullcorvec^2,corvec^2,...); abline(0,1, col='red')#; grid()
}

qqR <- function(corvec,nn,...)
{
## nn is the sample size, number of individuals used to compute correlation.
## needs correlation vector as input.
## nullcorvec generates a random sample from correlation distributions, under the null hypothesis of 0 correlation using Fisher's approximation.
  mm <- length(corvec)
  nullcorvec = tanh(rnorm(mm)/sqrt(nn-3)) ## null correlation vector
  qqplot(nullcorvec,corvec,...); abline(0,1)#; grid()
}


qqunif = 
  function(p,BH=T,CI=T,mlog10_p_thres=30,...)
  {
    ## thresholded by default at 1e-30
    p=na.omit(p)
    nn = length(p)
    xx =  -log10((1:nn)/(nn+1))
    
    p_thres = 10^{-mlog10_p_thres}
    if( sum( p < p_thres) )
    {
      warning(paste("thresholding p to ",p_thres) )
      p = pmax(p, p_thres)
    }
    plot( xx,  -sort(log10(p)),
          xlab=expression(Expected~~-log[10](italic(p))),
          ylab=expression(Observed~~-log[10](italic(p))),
          cex.lab=1,mgp=c(2,1,0),
          ... )
    abline(0,1,col='gray')
    if(BH)
    {
      abline(-log10(0.05),1, col='red',lty=1)
      abline(-log10(0.10),1, col='green',lty=2)
      abline(-log10(0.25),1, col='blue',lty=3)
      legend('right', c("FDR = 0.05","FDR = 0.10","FDR = 0.25"),
             col=c('red','green','blue'),lty=1:3, cex=1)
      abline(h=-log10(0.05/nn)) ## bonferroni
    }
    if(CI)
    {
      ## create the confidence intervals
      c95 <- rep(0,nn)
      c05 <- rep(0,nn)
      ## the jth order statistic from a
      ## uniform(0,1) sample
      ## has a beta(j,n-j+1) distribution
      ## (Casella & Berger, 2002,
      ## 2nd edition, pg 230, Duxbury)
      ## this portion was posted by anonymous on
      ## http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
      
      for(i in 1:nn)
      {
        c95[i] <- qbeta(0.95,i,nn-i+1)
        c05[i] <- qbeta(0.05,i,nn-i+1)
      }
      
      lines(xx,-log10(c95),col='gray')
      lines(xx,-log10(c05),col='gray')
    }
  }
