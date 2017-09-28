# conf_pval.fns.R stripped of all functions that aren't directly involved in computing the conf pval

#-----------------------------------------------------

getConfPvalLat <- function(x0, sample, conf=alpha, confT=qchisq(conf, 3, lower.tail=FALSE),
                        mins, maxs, minScale, maxScale) {
  
  # Returns the MLE-based point estimate of the p-value of x0 assuming a 3-Gamma null and tries
  # to find a lower confidence bound on it.
  
  targetF <- function(tht) {
    # Returns the 3-Gamma(tht) likelihood of x0 if tht lies within the confidence set
    # but 1+2llr otherwise
    llr2 <- 2 * (fit$loglike - sum( dgamma(sample-tht[1], scale=tht[2], shape=tht[3], log=TRUE) ) )
    if (llr2 > confT)
      -llr2
    else
      pgamma(x0-tht[1], scale=tht[2], shape=tht[3], lower.tail=FALSE)
  }
  fit <- fsgam3sb(sample, mins=mins, maxs=maxs)
  tht0 <- fit$Estimates
  names(tht0) <- NULL
  if (missing(minScale) || missing(maxScale)) {
    parLb <- c(scale=0, shape=mins*(1-1e-6))
    parUb <- c(shape=maxs*(1+1e-6), loc=min(sample)*(1-.Machine$double.eps))
  } else {
    parLb <- c(scale=min(tht0[2],minScale)*(1-1e-6), shape=mins*(1-1e-6))
    parUb <- c(scale=max(tht0[2],maxScale)*(1+1e-6), shape=maxs*(1+1e-6), loc=min(sample)*(1-.Machine$double.eps))
  }
  t <- consOptimPrep(c("loc", "scale", "shape"), parLb, parUb)
  maxP <- constrOptim(tht0, targetF, ui = t$M, ci = t$v, method="Nelder-Mead",
                      control=list(reltop=1e-8, maxit=1000, fnscale=-1))$value
  shapeLPs = c(mins*(1+1e-3), sqrt(mins*maxs))
  locLPs = c(tht0[1], min(sample)-tht0[2], min(sample)-tht0[2]*100)
  if (missing(minScale) || missing(maxScale))
    scaleLPs = c(0.1, 1, 10)
  else
    scaleLPs = c(minScale, maxScale)
#  keep <- "max at MLE" # debug
  for (loc in locLPs)
    for (shape in shapeLPs)
      for (scale in scaleLPs) {
        pv <- constrOptim(c(loc,scale,shape), targetF, ui = t$M, ci = t$v, method="Nelder-Mead",
                    control=list(reltop=1e-8, maxit=1000, fnscale=-1))$value
        if (pv > maxP) #{
          maxP <- pv
          #for debug
          #keep <- c(loc, scale, shape) }
      }
#  cat("max p-value", maxP, "at", keep, "\n")
  ePval = pgamma(x0-tht0[1], scale=tht0[2], shape=tht0[3], lower.tail=FALSE)
  cat(sprintf("\nThe MLE of the p-value is %.2g and its confidence interval is (0, %.2g)\n\n", ePval, maxP))
  c("maxPval" = maxP, "ePval"=ePval)
}


#-----------------------------------------------------
getMLEpVal <- function(x0, sample, mins=1, maxs=500) {
  
  # Returns the MLE-based point estimate of the p-value of x0 assuming a 3-Gamma null
  # mins,maxs are bounds on the shape parameter
  
  fit <- fsgam3sb(sample, mins=mins, maxs=maxs)
  tht0 <- fit$Estimates
  names(tht0) <- NULL
  c("ePval" = pgamma(x0-tht0[1], scale=tht0[2], shape=tht0[3], lower.tail=FALSE))
}


#	Fit a shifted gamma

fsgam3b <- function(data, minb=2*min(data)-max(data), maxb=min(data)*(1-.Machine$double.eps)) {
#
# Fit a shifted gamma data using the fact that (a,b) can be expressed in terms of nu
# at a critical point.
#
	cons_ll <- function(nu) {	# WARNING: if nu is ridiculously big the data washes out
		a <- - 1 / ( ( n / (mean(data) - nu) ) / sum(1/(data-nu)) - 1 )
		b <- - n / (1-a) * 1  / sum(1/(data-nu));
		sdata <- data - nu;
		sum(dgamma(sdata, shape=a, scale=b, log=TRUE))
	}
	maxb <- min(maxb, min(data)*(1-.Machine$double.eps)) # overide the user
	n <- length(data)
	maxstruc <- optimize(cons_ll, lower=minb, upper=maxb, maximum=TRUE, tol=1e-08)
	nu <- maxstruc$maximum
	ll <- maxstruc$objective
	a <- - 1 / ( ( n / (mean(data) - nu) ) / sum(1/(data-nu)) - 1 )
	b <- - n / (1-a) * 1  / sum(1/(data-nu));
	list("Estimates"=c("loc"=nu, "scale"=b, "shape"=a), "loglike"=ll)
}


fsgam3sb <- function(data, mins=1, maxs=1e6) {
#
# Translates shape bounds to location bounds and calls fsgam3b
#
	n <- length(data)
	a_crit <- function(nu) - 1 / ( ( n / (mean(data) - nu) ) / sum(1/(data-nu)) - 1 )
# Translate the shape bounds to ones on the location parameter
	if (mins == 1)
		maxl <- min(data)*(1-.Machine$double.eps)
	else
		maxl <- uniroot(function(nu) a_crit(nu) - mins, 
				c(-1e6,min(data)*(1-10*.Machine$double.eps)), tol=1e-8)$root
        lb <- -1e6
        while (a_crit(lb) < maxs)
          lb <- lb*10
	minl <- uniroot(function(nu) a_crit(nu) - maxs, 
				c(lb,min(data)*(1-10*.Machine$double.eps)), tol=1e-8)$root
	fsgam3b(data, minb=minl, maxb=maxl)
}


consOptimPrep <- function(pns, par_lb, par_ub) {
	M <- NULL; bs <- NULL		# <- vector(mode="numeric")
	for (i in 1:length(pns)) {
		t <- par_bnd(par_lb, pns[i], i, TRUE)
		M <- rbind(M, t$Mrow)
		bs <- c(bs, t$entry)
	}
	for (i in 1:length(pns)) {
		t <- par_bnd(par_ub, pns[i], i, FALSE)
		M <- rbind(M, -t$Mrow)
		bs <- c(bs, -t$entry)
	}
	list("M"=M, "v"=bs)
}


par_bnd <- function(par_b, pname, pnum, lowerb) {
	if (is.na(par_b[pname])) {
		b_row <- c(0,0,0)
		if (lowerb)
			b <- -1
		else
			b <- 1
	} else {
		b_row <- c(0,0,0)
		b_row[pnum] <- 1
		b <- par_b[pname]
	}
	list("Mrow"=b_row, "entry"=b)
}

