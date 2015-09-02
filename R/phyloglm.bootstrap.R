bootstrap <- function(formula, data, phy,
                      method = c("logistic_MPLE", "logistic_IG10"),
                      btol = 10, log.alpha.bound = 4, start.beta = NULL,
                      start.alpha = NULL, B = 100, full.matrix = FALSE)
{
  fit <- phyloglm(formula = formula, data = data, phy = phy, method = method, btol = btol, log.alpha.bound = log.alpha.bound, start.beta = start.beta, start.alpha = start.alpha)

  # simulate all bootstrap data sets
  mf <- model.frame(formula = formula, data = data)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  bootobject <- rbinTrait(n = B, phy = phy, beta = fit$coefficients,
                          alpha = fit$alpha, X = X, model = "LogReg")
  responsecolumn <- attr(attr(mf, "terms"), "response")

  # analyze these bootstrapped data
  ncoeff = length(fit$coefficients)
  bootmatrix <- matrix(NA, B, ncoeff + 1)
  colnames(bootmatrix) <- c(names(fit$coefficients), "alpha")
  for (i in 1:B){
    mf[, responsecolumn] = bootobject[,i]
    bootfit <- try(phyloglm(formula = formula, data = mf, phy = phy,
                            method = method, btol = btol,
                            log.alpha.bound = log.alpha.bound,
                            start.beta = start.beta, start.alpha = start.alpha),
                   silent=TRUE)
    if (!inherits(bootfit, 'try-error')){
      bootmatrix[i, 1:ncoeff] <- bootfit$coefficients
      bootmatrix[i, ncoeff + 1] <- bootfit$alpha
    }
  }
	
  # summarize bootstrap estimates
  ind.na <- which(is.na(bootmatrix[,1]))
  # indices of replicates that failed: phyloglm had an error
  if (length(ind.na)>0) {
    bootmatrix <- bootmatrix[-ind.na,]
    numOnes <- range(apply(bootobject[,ind.na],2,sum))
  }
  means <- apply(bootmatrix, 2, mean)
  sd <- apply(bootmatrix, 2, sd)
  confint95 <- apply(bootmatrix, 2, quantile, probs = c(.025, .975))
	
  meanAlog <- mean(log(bootmatrix[, ncoeff + 1]))
  sdAlog <- sd(log(bootmatrix[, ncoeff + 1]))	

  results <- list(mean=means, sd=sd, confint95=confint95,
                  log.alpha.mean = meanAlog, log.alpha.sd = sdAlog,
                  numFailed = length(ind.na))
  if(full.matrix)
    results$bootstrap = bootmatrix
  return(results)
}			
