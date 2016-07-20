phylolm <- function(formula, data=list(), phy, 
	model=c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"),
	lower.bound=NULL, upper.bound=NULL, starting.value=NULL, measurement_error = FALSE,
	boot=0,full.matrix = TRUE, ...)
{

  ## initialize	
  if (!inherits(phy, "phylo")) stop("object \"phy\" is not of class \"phylo\".")
  model = match.arg(model)	
  if ((model=="trend")&(is.ultrametric(phy)))
   stop("the trend is unidentifiable for ultrametric trees.")
  if ((model=="lambda") && measurement_error)
    stop("the lambda transformation and measurement error cannot be used together: they are not distinguishable")
  if (is.null(phy$edge.length)) stop("the tree has no branch lengths.")
  if (is.null(phy$tip.label)) stop("the tree has no tip labels.")	
  tol = 1e-10	
  phy = reorder(phy,"pruningwise")
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  externalEdge <- des<=n
  mf = model.frame(formula=formula,data=data)
  if (nrow(mf)!=length(phy$tip.label))
    stop("number of rows in the data does not match the number of tips in the tree.")
  if (is.null(rownames(mf))) {
   warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
   data.names = phy$tip.label 
  }
  else data.names = rownames(mf)
  order = match(data.names, phy$tip.label)
  if (sum(is.na(order))>0) {
   warning("data names do not match with the tip labels.\n")
   rownames(mf) = data.names
  } else {
   tmp = mf
   rownames(mf) = phy$tip.label
   mf[order,] = tmp[1:nrow(tmp),]
  }
  X = model.matrix(attr(mf, "terms"), data=mf)
  y = model.response(mf)
  d = ncol(X)
  OU = c("OUrandomRoot","OUfixedRoot")
  flag = 0 # flag and D are used for OU model if tree is not ultrametric:
  D = NULL #            for the generalized 3-point structure

  ## preparing for OU model
  if (model %in% OU) {
    D = numeric(n)
    if (!is.ultrametric(phy)) {
      flag = 1
      dis = pruningwise.distFromRoot(phy)
      Tmax = max(dis[1:n])
      D = Tmax - dis[1:n]
      D = D - mean(D)
      phy$edge.length[externalEdge] <- phy$edge.length[externalEdge] + D[des[externalEdge]]
      ## phy is now ultrametric, with height Tmax:
      Tmax = Tmax + min(D)
    }
  }
	
  ## preparing for trend model
  if (model == "trend") {
    trend = pruningwise.distFromRoot(phy)[1:n]
    X = cbind(X,trend)
    d = d+1
  }

  ## calculate Tmax = average distance from root to tips,
  ## to choose appropriate starting values later # fixit
  dis = pruningwise.distFromRoot(phy)[1:n]
  Tmax = mean(dis)

  ## Default bounds
  bounds.default = matrix(c(1e-7/Tmax,50/Tmax,1e-7,1,1e-6,1,1e-5,3,-3/Tmax,0,1e-16,1e16), ncol=2, byrow=TRUE)
  rownames(bounds.default) = c("alpha","lambda","kappa","delta","rate","sigma2_error")
  colnames(bounds.default) = c("min","max")

  ## Default starting values
  starting.values.default = c(0.5/Tmax,0.5,0.5,0.5,-1/Tmax,1) 
  names(starting.values.default) = c("alpha","lambda","kappa","delta","rate","sigma2_error")

  ## User defined bounds and starting values
  if (is.null(lower.bound)) {
    if (model %in% OU) 
      lower.bound = bounds.default[1,1]
    if (model=="lambda") lower.bound = bounds.default[2,1]
    if (model=="kappa") lower.bound = bounds.default[3,1]
    if (model=="delta") lower.bound = bounds.default[4,1]
    if (model=="EB") lower.bound = bounds.default[5,1]
  }
  if (is.null(upper.bound)) {
    if (model %in% OU) 
      upper.bound = bounds.default[1,2]
    if (model=="lambda") upper.bound = bounds.default[2,2]
    if (model=="kappa") upper.bound = bounds.default[3,2]
    if (model=="delta") upper.bound = bounds.default[4,2]
    if (model=="EB") upper.bound = bounds.default[5,2]
  }
  if (is.null(starting.value)) {
    if (model %in% OU) 
      starting.value = starting.values.default[1]
    if (model=="lambda") starting.value = starting.values.default[2]
    if (model=="kappa") starting.value = starting.values.default[3]
    if (model=="delta") starting.value = starting.values.default[4]
    if (model=="EB") starting.value = starting.values.default[5]
  }

  if (measurement_error) {
    lower.bound = c(lower.bound, bounds.default[6,1])
    upper.bound = c(upper.bound, bounds.default[6,2])
    starting.value = c(starting.value, starting.values.default[6])
  }


  ## preparing for general use of "parameter" for branch length transformation
  prm = list(myname = starting.value[1])
  names(prm) = model # good for lambda, kappa, delta
  if (model %in% OU) names(prm) = "alpha"
  if (model == "EB") names(prm) = "rate"
  if (measurement_error) {
    if (model %in% c("BM","trend")) names(prm) = "sigma2_error"
    else prm[["sigma2_error"]] = starting.value[2]
  }

  ## log-likelihood, computation using the three-point structure
  ole= 4 + 2*d + d*d # output length
  loglik <- function(parameters,y,X) {
    tree = transf.branch.lengths(phy,model,parameters=parameters,
      check.pruningwise=F,check.ultrametric=F,D=D,check.names=F)$tree
    if (flag) { # need diagonal terms in D
      y = exp(-parameters$alpha*D) * y # variables local to this function
      X = exp(-parameters$alpha*D) * X
    }
    tmp=.C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
      as.integer(1),as.integer(d),as.integer(ROOT),as.double(tree$root.edge),as.double(tree$edge.length),
      as.integer(des), as.integer(anc), as.double(as.vector(y)), as.double(as.vector(X)),
      result=double(ole))$result # tmp has, in this order:
    ## logdetV, 1'V^{-1}1, y'V^{-1}1, y'V^{-1}y, X'V^{-1}1, X'V^{-1}X, X'V^{-1}y
    comp = list(vec11=tmp[2], y1=tmp[3], yy=tmp[4], X1=tmp[5:(4+d)],
                XX=matrix(tmp[(5+d):(ole-d)], d,d),Xy=tmp[(ole-d+1):ole],logd=tmp[1])
    invXX = solve(comp$XX)
    betahat = invXX%*%comp$Xy
    sigma2hat = as.numeric((comp$yy - 2*t(betahat)%*%comp$Xy + t(betahat)%*%comp$XX%*%betahat)/n)
    if (sigma2hat<0) {
      resdl = X%*%betahat - y
      tmpyy=.C("threepoint", as.integer(N),as.integer(n),as.integer(phy$Nnode),
        as.integer(1),as.integer(d),as.integer(ROOT),as.double(tree$root.edge),as.double(tree$edge.length),
        as.integer(des), as.integer(anc), as.double(as.vector(resdl)), as.double(as.vector(X)),
        result=double(ole))$result[4]
      sigma2hat = tmpyy/n
    }
    n2llh = as.numeric( n*log(2*pi) + n + n*log(sigma2hat) + comp$logd) # -2 log-likelihood
    if (flag)
      n2llh = n2llh + parameters$alpha * 2*sum(D)
    ## because diag matrix used for generalized 3-point structure is exp(alpha diag(D))
    vcov = sigma2hat*invXX*n/(n-d)
    return(list(n2llh=n2llh, betahat = as.vector(betahat), sigma2hat=sigma2hat,vcov=vcov))
  }

  ## Fitting
  lower = lower.bound
  upper = upper.bound
  start = starting.value

  if ((model %in% c("BM","trend"))&&(!measurement_error)) {
    BMest = loglik(prm, y, X) # root edge taken to be 0
    results <- list(coefficients=BMest$betahat, sigma2=BMest$sigma2hat, optpar=NULL, sigma2_error = 0,
                    logLik=-BMest$n2llh/2, p=1+d, aic=2*(1+d)+BMest$n2llh, vcov = BMest$vcov)
  } else {
    ##------- Optimization of phylogenetic correlation parameter is needed ---------#
    # first: checks of bounds
    if (sum(lower>start)+sum(upper<start)>0 )
      stop("The starting value is not within the bounds of the parameter.")

    # def of function to be optimized
    # minus2llh_sinvar: return the -loglik given a single variable: measure-err
    #    logvalue=log of measure-err variance
    #    if not BM: prm[[1]] needs up-to-date before calling this fcn:
    #               prm[[1]] = alpha or lamda or kappa or rate
    minus2llh_sinvar=function(logvalue) {
      if(model %in% c("BM","trend")){
        prm[[1]] = exp(logvalue)
      } else
        prm[[2]] = exp(logvalue)
      loglik(prm, y, X)$n2llh
    }
    # minus2llh: returns -loglik given 'alpha' and possibly measurement error
    minus2llh <- function(logvalue) {
      if (model == "EB") prm[[1]]=logvalue[1] # which is 'rate', not log(rate)
      else prm[[1]]=exp(logvalue[1]) # first element is the parameter of the model
      if ((measurement_error)&&(!(model %in% c("BM","trend")))){
        prm[[2]]=exp(logvalue[2])
      }
      loglik(prm, y, X)$n2llh
    }

    # get objects to start the search, and to bound the search
    if (lower[1]==upper[1] && !(measurement_error)) { # no optimization in fact
      prm[[1]] = lower[1]
      BMest = loglik(prm, y, X)
    }  else {
      if(model !="EB"){   # logstart = (log(alpha/lambda...), log(m.e. variance)) or just log(alpha,...)
        logstart = log(start)
        loglower = log(lower)
        logupper = log(upper)
      } else{
        # logstart = (rate, log(m.e. variance)) or just rate
        logstart = start  # do *not* log transform 'rate' because it is <= 0
        loglower = lower
        logupper = upper
        if(measurement_error){
          logstart[2] = log(start[2])
          loglower[2] = log(lower[2])
          logupper[2] = log(upper[2])
        }
      }
      if (!(model %in% c("BM","trend")) && lower[1] != upper[1]){
        opt <- optim(logstart, fn = minus2llh, method = "L-BFGS-B",lower=loglower, upper = logupper, ...)
        # next: get (co)variance parameters alpha/lambda... and m.e. variance into "prm"
        #       to be used later in loglik to get the estimated beta and sigma2.
        if (model == "EB") MLEvalue = as.numeric(opt$par[1]) else MLEvalue = as.numeric(exp(opt$par[1]))
        prm[[1]] = MLEvalue
        if ((isTRUE(all.equal(MLEvalue,lower[1], tol=tol)))||(isTRUE(all.equal(MLEvalue,upper[1],tol=tol)))) {
          matchbound = 1
          if ((model %in% c("lambda","kappa"))&&(MLEvalue == 1)) matchbound=0
          if ((model == "EB")&&(MLEvalue == 0)) matchbound=0
          if (matchbound)
            warning(paste("the estimation of", names(prm)[1],
                          'matches the upper/lower bound for this parameter.
                          You may change the bounds using options "upper.bound" and "lower.bound".\n'))
        }
        if (measurement_error){
            MLEsigma2_error = as.numeric(exp(opt$par[2]))
            prm[[2]] = MLEsigma2_error
        }
      } else { # then there must be measurement error
          # if BM or trend, logstart = log(m.e. variance), already set above.
          if (!(model %in% c("BM","trend"))) { # then we must have (lower[1]==upper[1])
            prm[[1]]=lower[1]
            logstart=logstart[2]
            loglower=loglower[2]
            logupper=logupper[2]
          }
          opt <- optim(logstart, fn = minus2llh_sinvar,method = "L-BFGS-B",lower=loglower, upper = logupper, ...)
          MLEsigma2_error = as.numeric(exp(opt$par[1]))
          if (model %in% c("BM","trend")){
            prm[[1]] = MLEsigma2_error
          } else {
            prm[[2]] = MLEsigma2_error
          }
      }
      # estimate beta and sigma2:
      BMest = loglik(prm, y, X)
    }
    # rescaling measurement error by sigma2
    sigma2_errorhat = 0
    if (measurement_error)
     sigma2_errorhat = MLEsigma2_error * BMest$sigma2hat
    # rescale sigma2 if OU, because it was "gamma" originally: sigma2 = 2 alpha gamma:
    if (model %in% OU)
      BMest$sigma2hat = 2*prm[[1]] * BMest$sigma2hat
      results <- list(coefficients=BMest$betahat, sigma2=BMest$sigma2hat, optpar=prm[[1]], sigma2_error = sigma2_errorhat,
                    logLik=-BMest$n2llh/2, p=2+d, aic=2*(2+d)+BMest$n2llh, vcov = BMest$vcov)
    if (model %in% c("BM","trend")) { # adjust for BM and trend models
      results$optpar = NULL
      results$p = results$p - 1
      results$aic = results$aic - 2
    }
    if (measurement_error) {
      results$p = results$p + 1 # adjust the number of parameters
      results$aic = results$aic + 2 # adjust AIC value
    }
  } # end of cases that are not BM only

  names(results$coefficients) = colnames(X)
  colnames(results$vcov) = colnames(X)
  rownames(results$vcov) = colnames(X)
  results$fitted.values = drop(X %*% results$coefficients)
  # drop: simplifies to vector, instead of matrix with 1 column
  results$residuals = y - results$fitted.values
  results$mean.tip.height = Tmax
  results$y = y
  results$X = X
  results$n = n
  results$d = d
  results$formula = formula
  results$call = match.call()
  results$model = model
  results$boot = boot

  ## starting the bootstrap
  if (boot>0) {
  # Turn off warnings
    options(warn=-1)
    # simulate all bootstrap data sets
    simmodel=model
    OU = c("OUrandomRoot","OUfixedRoot")
    if (model %in% OU){
      simmodel="OU"
    }
    # prm has "alpha" and "sigma2_error". Here we want "alpha" and "sigma2"
    prmsimul = list(sigma2 = results$sigma2)
    if (!(simmodel %in% c("BM","trend"))){
      Noname=prm[[1]]
      names(Noname)=NULL
      prmsimul[[2]] = Noname
      names(prmsimul)[2] = names(prm)[1]
    }

    # booty = bootstrap y response. parameter list depends on the model
    booty <- results$fitted.values +
     rTrait(n = boot, phy = phy,model = simmodel, parameters=prmsimul)
     # by default: ancestral.state=0 and optimal.value=0 --used for OU model only
    if (measurement_error){
      booty <- booty + rnorm(boot*n, mean=0, sd=sqrt(results$sigma2_error))
    }

    # analyze these bootstrapped data
    ncoeff = length(results$coefficients)
    colnumberalpha=ncoeff+2
    # nparam_var: number of parameters that affect the variances and covariances
    nparam_var = 1 # for sigma2
    if (measurement_error){
        nparam_var = nparam_var+1
        colnumberalpha=colnumberalpha+1
    }
    if (!(model %in% c("BM","trend"))) { nparam_var=nparam_var+1}
    # initialize bootmatrix: matrix with bootstrap estimated parameters
    bootmatrix <- matrix(NA, boot, ncoeff + nparam_var)
    colnames(bootmatrix) <- paste0("v",1:(ncoeff + nparam_var)) # temporary names, to initialize
    colnames(bootmatrix)[1:ncoeff] <- names(results$coefficients)
    colnames(bootmatrix)[ncoeff+1] <- "sigma2"
    if (measurement_error){
      colnames(bootmatrix)[ncoeff+2] <- "sigma2_error"
    }
    if (!(model %in% c("BM","trend"))) {
      colnames(bootmatrix)[colnumberalpha] <- names(prm)[1]
    }

    for (i in 1:boot){
      y = booty[,i]
      # all other cases: first get 'prm' up-to-date.
      if (!(model %in% c("BM","trend")) && lower[1]==upper[1])
        bootmatrix[i,colnumberalpha] = lower[1] # will be used later to update 'prm'

      # otherwise: something needs to be optimized: m.e., alpha, or both.
      # below: storing 'standardized' estimated of m.e. variance in MLEsigma2_error. Rescaled later.
      if (!(model %in% c("BM","trend")) && lower[1] != upper[1]){
        try(opt <- optim(logstart, fn = minus2llh, method = "L-BFGS-B",lower=loglower, upper = logupper, ...),silent=TRUE)
        if (!inherits(opt, 'try-error')){
          if (model == "EB") {
            bootmatrix[i,colnumberalpha] = as.numeric(opt$par[1]); # will be used later to update 'prm'
          } else {
            bootmatrix[i,colnumberalpha] = as.numeric(exp(opt$par[1]));
          }
          if(measurement_error)  MLEsigma2_error = as.numeric(exp(opt$par[2]))
        } else {
          if(measurement_error){
            try(opt <- optim(logstart, fn = minus2llh_sinvar,method = "L-BFGS-B",lower=loglower, upper = logupper, ...),silent = TRUE)
            if (!inherits(opt, 'try-error')){
              MLEsigma2_error = as.numeric(exp(opt$par[1]))
            }
          }
        }
      }

      if (!(model %in% c("BM","trend")))
        prm[[1]] = bootmatrix[i,colnumberalpha] # update of 'prm'
      if (measurement_error){
        if (!(model %in% c("BM","trend"))) {
          prm[[2]] = MLEsigma2_error
        } else {
          prm[[1]] = MLEsigma2_error
        }
      }
      BMest = loglik(prm, y, X)
      if (model %in% OU)
             BMest$sigma2hat = 2*prm[[1]] * BMest$sigma2hat # was "gamma" originally: sigma2 = 2 alpha gamma
      if (measurement_error)
             bootmatrix[i,ncoeff+2] <- MLEsigma2_error * BMest$sigma2hat
            bootmatrix[i,1:ncoeff] <- BMest$betahat
            bootmatrix[i,ncoeff+1] <- BMest$sigma2hat
    } # end of loop over bootstrap reps
    # summarize bootstrap estimates
    ind.na <- which(is.na(bootmatrix[,1]))
    # indices of replicates that failed: phylolm had an error
    if (length(ind.na)>0) {
      bootmatrix <- bootmatrix[-ind.na,]
      numOnes <- range(apply(booty[,ind.na],2,sum))
    }
    bootmean <- apply(bootmatrix, 2, mean)
    bootsd <- apply(bootmatrix, 2, sd)
    bootconfint95 <- apply(bootmatrix, 2, quantile, probs = c(.025, .975))
    bootmeansdLog = matrix(NA,2,nparam_var) # will give NaN for rate, but won't be printed.
    colnames(bootmeansdLog) = colnames(bootmatrix)[(ncoeff+1):(ncoeff+nparam_var)]
    for (i in 1:nparam_var){
      bootmeansdLog[1,i] <- mean(log(bootmatrix[, ncoeff + i]))
      bootmeansdLog[2,i] <- sd(log(bootmatrix[, ncoeff + i]))
    }

    results$bootmean = bootmean
    results$bootsd = bootsd
    results$bootconfint95 = bootconfint95
    results$bootmeansdLog = bootmeansdLog
    results$bootnumFailed = length(ind.na)
    if (full.matrix) results$bootstrap = bootmatrix

    ### Turn on warnings
    options(warn=0)
  }
  class(results) = "phylolm"
  return(results)
}

################################################
################################################

print.phylolm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  aiclogLik = c(x$aic,x$logLik)
  names(aiclogLik) = c("AIC","logLik")
  print(aiclogLik, digits = digits)
  cat("\nParameter estimate(s) using ML:\n")
  if (!is.null(x$optpar)) {
    if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha:",x$optpar)
    if (x$model %in% c("lambda","kappa","delta")) cat(x$model,":",x$optpar)
    if (x$model=="EB") cat("rate:",x$optpar)
    cat("\n")
  }
  cat("sigma2:",x$sigma2,"\n")
  if (x$sigma2_error > 0) cat("sigma2_error:",x$sigma2_error,"\n")
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
################################################
summary.phylolm <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  tval <- coef(object) / se

  if (object$boot == 0)
    TAB <- cbind(Estimate = coef(object), StdErr = se, t.value = tval,
                 p.value = 2*pt(-abs(tval), df=object$n - object$d))
  else
    TAB <- cbind(Estimate = coef(object), StdErr = se, t.value = tval,
                 lowerbootCI = object$bootconfint95[1,1:object$d],
                 upperbootCI = object$bootconfint95[2,1:object$d],
                 p.value = 2*pt(-abs(tval), df=object$n - object$d)) # need p-value last for printCoefmat

  res <- list(call=object$call, coefficients=TAB,
              residuals = object$residuals, sigma2 = object$sigma2,
              optpar=object$optpar, sigma2_error = object$sigma2_error, logLik=object$logLik,
              df=object$p, aic=object$aic, model=object$model,
              mean.tip.height=object$mean.tip.height,
              bootNrep = ifelse(object$boot>0, object$boot - object$bootnumFailed, 0))
  if (res$bootNrep>0) {
    res$bootmean = object$bootmean
    res$bootsd = object$bootsd
    res$bootconfint95 = object$bootconfint95
    res$bootmeansdLog <- object$bootmeansdLog
  }
  class(res) = "summary.phylolm"
  res
}
################################################
print.summary.phylolm <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  aiclogLik = c(x$aic,x$logLik)
  names(aiclogLik) = c("AIC","logLik")
  print(aiclogLik, digits = digits)
  r <- zapsmall(quantile(x$residuals), digits + 1)
  names(r) <- c("Min", "1Q", "Median", "3Q", "Max")
  cat("\nRaw residuals:\n")
  print(r, digits = digits)

  cat("\nMean tip height:",x$mean.tip.height)
  cat("\nParameter estimate(s) using ML:\n")
  if (!is.null(x$optpar)) {
    if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha:",x$optpar)
    if (x$model %in% c("lambda","kappa","delta")) cat(x$model,":",x$optpar)
    if (x$model=="EB") cat("rate:",x$optpar)
    cat("\n")
  }

  cat("sigma2:",x$sigma2,"\n")
  if (x$sigma2_error > 0) cat("sigma2_error:",x$sigma2_error,"\n")

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
  if (!is.null(x$optpar)) {
    cat("\nNote: p-values are conditional on ")
    if (x$model %in% c("OUrandomRoot","OUfixedRoot")) cat("alpha=",x$optpar,".",sep="")
    if (x$model %in% c("lambda","kappa","delta")) cat(x$model,"=",x$optpar,".",sep="")
    if (x$model=="EB") cat("rate=",x$optpar,".",sep="")
    cat("\n")
  }

  if (x$bootNrep > 0) {
    cat("\n")
    ncoef = nrow(x$coefficients)
    nparam_var=length(x$bootmean)-ncoef # number of variance/covariance parameters

    for (i in 1:nparam_var){
      tmp = colnames(x$bootmeansdLog)[i]# sigma2, sigma2_error or lambda/alpha/...
      tmp = ifelse(tmp %in% c("sigma2", "sigma2_error"), tmp, "optpar")
      cat(colnames(x$bootmeansdLog)[i],": ",x[[tmp]],"\n", sep="")
      cat("      bootstrap mean: ",    x$bootmean[ncoef+i]   ," (on raw scale)","\n",sep="")
      if (!(x$model %in% c("EB","lambda")) || tmp != "optpar")
        cat("                      ",exp(x$bootmeansdLog[1,i])," (on log scale, then back transformed)","\n",sep="")
      cat("      bootstrap 95% CI: (",x$bootconfint95[1,ncoef+i],",",x$bootconfint95[2,ncoef+i],")\n\n", sep="")
    }
    cat("Parametric bootstrap results based on",x$bootNrep,"fitted replicates\n")
  }
}
################################################
residuals.phylolm <-function(object,type=c("response"), ...){
  type <- match.arg(type)
  object$residuals	 
}
################################################
vcov.phylolm <- function(object, ...){
  object$vcov
}
################################################
logLik.phylolm <- function(object, ...){
  res = list(logLik = object$logLik, df = object$p)
  class(res) = "logLik.phylolm"
  res
}
print.logLik.phylolm <- function (x, ...) {
  cat("'log Lik.' ",x$logLik," (df=",x$df,")\n", sep = "")
}
AIC.logLik.phylolm <- function(object, k=2, ...) {
  return(k*object$df - 2*object$logLik)
}
AIC.phylolm <- function(object, k=2, ...) {
  return(AIC(logLik(object),k))
}
extractAIC.phylolm <- function(fit, scale, k=2, ...) {
    c(fit$p, - 2*fit$logLik + k * fit$p)
}
nobs.phylolm <- function(object, ...){
  return(object$n)
}
################################################
predict.phylolm <- function(object, newdata=NULL, ...){
  if (object$model=="trend")
    stop("Predicting for trend model has not been implemented.")
  if(is.null(newdata)) y <- fitted(object)
  else{			
    X = model.matrix(delete.response(terms(formula(object))),data = newdata)
    y <- X %*% coef(object)
  }
  y
}
################################################
plot.phylolm <-function(x, ...){
  plot(x$y, fitted(x), xlab = "Observed value", ylab = "Fitted value", ...)
}
################################################
