################################################
### Stepwise selection for phylolm using AIC
################################################

phyloglmstep <- function(formula, starting.formula = NULL, data=list(), phy, 
                      method = c("logistic_MPLE","logistic_IG10"),
                      direction = c("both", "backward", "forward"), trace = 2,
                      btol = 10, log.alpha.bound = 4, start.beta=NULL, 
                      start.alpha=NULL, boot = 0, full.matrix = TRUE,
                      k=2, ...) 
{
    ## initialize  
    method = match.arg(method)	
    direction = match.arg(direction)  
    fit.full = phyloglm(formula, data, phy, method, btol, log.alpha.bound, start.beta, start.alpha, boot, full.matrix, ...)
    response = fit.full$formula[[2]] # name of the response
    covariates = attr(terms(formula), "term.labels") # name of the covariates
    
    p = length(covariates)
    if (p==0) stop("Your full model only has intercept. Model selection is not needed.")
    plm.full = rep(1,p)
    
    ## create a formula from current model
    create.formula <- function(plm) {
        on = which(plm==1)
        str = paste(response," ~ 1",sep="")
        if (length(on) > 0) 
            for (i in 1:length(on)) str = paste(str," + ",covariates[on[i]],sep="")
        return(str)
    }
    
    ## fit phylolm
    fit <- function(plm) {
        return(phyloglm(create.formula(plm), data, phy, method, btol, log.alpha.bound, start.beta, start.alpha, boot, full.matrix, ...))
    }
    
    ## plm.current is a binary vector of length p
    ## where 0 at i-th position means excluding i-th covariate
    
    if (direction == "forward") {
        plm.current = rep(0,p) # only intercept
        fit.current = fit(plm.current)
    } else {
        plm.current = plm.full # all covariates
        fit.current = fit.full
    }
    
    if (!is.null(starting.formula)) {
        fit.current = phylolm(starting.formula, data, phy, method, btol, log.alpha.bound, start.beta, start.alpha, boot, full.matrix, ...)
        covariates.current = attr(terms(starting.formula), "term.labels")
        plm.current = rep(0,p)
        position = match(covariates.current,covariates)
        if (any(is.na(position))) stop("The starting model is not a submodel of the full model.")
        plm.current[position] = 1
    }
    
    if (trace>0) {
        cat("----------\n")
        cat(paste("Starting model: ",create.formula(plm.current),"\n",sep=""))
        cat(paste("Direction: ",direction,"\n",sep=""))
        cat(paste("AIC(k=",k,"): ",AIC(fit.current,k),"\n",sep=""))
        #     cat("----------\n")
    }
    
    flag = 0 # flag of termination
    count = 1
    while (flag == 0) {
        flag = 1
        plm.best = plm.current
        fit.best = fit.current
        
        ### to preserve hierarchy priciple
        terms.add = add.scope(formula(create.formula(plm.current)), 
                              formula(create.formula(plm.full)))
        terms.drop = drop.scope(formula(create.formula(plm.current)))
        
        for (i in 1:p) {
            
            plm.propose = plm.current
            do.update = FALSE
            if ((plm.current[i]==1)&&(direction %in% c("both", "backward"))
                &&(covariates[i] %in% terms.drop)) {
                plm.propose[i] = 0  # remove i-th covariate
                do.update = TRUE
            }
            if ((plm.current[i]==0)&&(direction %in% c("both", "forward"))
                &&(covariates[i] %in% terms.add)) {
                plm.propose[i] = 1 # add i-th covariate
                do.update = TRUE
            }
            
            ## check if proposed model is better
            if (do.update){
                fit.propose = fit(plm.propose)
                if (trace>1) {
                    cat(paste0("\tProposed: ",create.formula(plm.propose),"\n"))
                    cat(paste0("\tAIC(k=",k,"): ",AIC(fit.propose,k),"\n"))
                }
                if (AIC(fit.propose,k) < AIC(fit.best,k)) {
                    plm.best = plm.propose
                    fit.best = fit.propose
                    flag = 0
                }
            }
        }
        
        ## Set current model as the best model
        plm.current = plm.best
        fit.current = fit.best
        
        if (trace>0) {
            cat("----------\nStep ",count,"\n",sep="")
            cat(paste("Current model: ",create.formula(plm.current),"\n",sep=""))
            cat(paste("AIC(k=",k,"): ",AIC(fit.current,k),"\n",sep=""))
        }
        count = count + 1
    }
    if (trace>0) cat("---END\n")
    return(fit.current)
}  
