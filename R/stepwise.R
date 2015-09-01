################################################
### Stepwise selection for phylolm using AIC
################################################

phylostep <- function(formula, starting.formula = NULL, data=list(), phy, 
                    model=c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"),
                    direction = c("both", "backward", "forward"), trace = TRUE,
                    lower.bound=NULL, upper.bound=NULL, starting.value=NULL, ...) 
{
  ## initialize  
  model = match.arg(model)	
  direction = match.arg(direction)  
  plm.full = phylolm(formula, data, phy, model, lower.bound, upper.bound, starting.value, ...)
  response = plm.full$formula[[2]] # name of the response
  covariates = names(plm.full$coefficients) # name of the covariates
  if (length(data) == 0) {
    data = as.data.frame(plm.full$X)
  } # check if there is no data structure
  
  if (is.null(covariates)) stop("Covariates has no name.")
  p = length(covariates)
  if (p==1) stop("Your full model only has intercept. Model selection is not needed.")

  ## create a formula from current model
  create.formula <- function(plm) {
    on = which(plm==1)
    str = paste(response," ~ 1",sep="")
    for (i in 1:length(on))
      if (i>1) str = paste(str," + ",covariates[i],sep="")
    return(str)
  }
  
  ## fit phylolm
  fit <- function(plm) {
    return(phylolm(create.formula(plm), data, phy, model, 
                   lower.bound, upper.bound, starting.value, ...))
  }
  
  ## plm.current is a binary vector of length p
  ## where 0 at i-th position means excluding i-th covariate
  
  if (direction == "forward") {
    plm.current = c(1,rep(0,p-1)) # only intercept
    fit.current = fit(plm.current)
  } else {
    plm.current = rep(1,p) # all covariates
    fit.current = plm.full
  }
  
  if (!is.null(starting.formula)) {
    fit.current = phylolm(starting.formula, data, phy, model, lower.bound, upper.bound, starting.value, ...)
    covariates.current = names(fit.current$coefficients)
    plm.current = c(1,rep(0,p-1))
    if (length(covariates.current)>1)
      for (i in 2:length(covariates.current)) 
        plm.current[which(covariates==covariates.current[i])] = 1  
  }

if (trace) {
  cat("----------\n")
  cat(paste("Starting model: ",create.formula(plm.current),"\n",sep=""))
  cat(paste("Direction: ",direction,"\n",sep=""))
  cat(paste("AIC: ",AIC(fit.current),"\n",sep=""))
  cat("----------\n")
}
    
  flag = 0 # flag of termination
  count = 1
  while (flag == 0) {
    flag = 1
    plm.best = plm.current
    fit.best = fit.current
    for (i in 2:p) {
      
      plm.propose = plm.current
      if ((plm.current[i]==1)&&(direction %in% c("both", "backward"))) 
        plm.propose[i] = 0  # remove i-th covariate
      if ((plm.current[i]==0)&&(direction %in% c("both", "forward"))) 
        plm.propose[i] = 1 # add i-th covariate
        
      ## check if proposed model is better
      fit.propose = fit(plm.propose)
      if (AIC(fit.propose) < AIC(fit.best)) {
        plm.best = plm.propose
        fit.best = fit.propose
        flag = 0
        }
      }
    
    ## Set current model as the best model
    plm.current = plm.best
    fit.current = fit.best
    
    if (trace) {
      cat("Step ",count,"\n",sep="")
      cat(paste("Current model: ",create.formula(plm.current),"\n",sep=""))
      cat(paste("AIC: ",AIC(fit.current),"\n",sep=""))
      cat("---\n")
    }
    count = count + 1
    }
  if (trace) cat("END\n")
  return(fit.current)
}  