mace <- function(trait, phy, lambda = NULL) {
  
  trait = as.matrix(trait)
  m = ncol(trait) # each column is a trait
  n = nrow(trait) # number of species
  C = cor(trait, use = "pairwise") # trait correlation matrix
  csum = sum(C)
  cmax = max(C-diag(1,m))
  
  standardize <- function() {
    
    Z = matrix(NA, nrow = n, ncol = m)
    Ys2 = rep(NA, m)
    Ymu = rep(NA, m)
    
    for (i in 1:m) {
      fit = phylolm(trait[,i] ~ 1, phy = phy)
      Ys2[i] = fit$sigma2
      Ymu[i] = fit$coefficients
      Z[,i] = trait[,i]/sqrt(Ys2[i])
    }
    
    return(list(Z = Z, Ymu = Ymu, Ys2 = Ys2))
    
  }
  
  stand = standardize() # standardize
  Zmu = stand$Ymu/sqrt(stand$Ys2) # \mu^{mle,\mb{Z}}
  
  t = transf.branch.lengths(phy = phy, model="BM")
  tmp = three.point.compute(t$tree,trait[,1],trait[,1],t$diagWeight)
  v = tmp$vec11 # compute 1'V^{-1}1
  
  # upper bound for lambda
  lambda_upper = (m-1)*(1-cmax)/((m-1)^2*(max(Zmu) - min(Zmu))^2 + (m^2 - csum)/v)
  
  if (is.null(lambda)) lambda = lambda_upper/2
  tmp = 2*lambda*(sum(Zmu))/(v + 2*lambda*m)
  res = v*stand$Ymu/(v + 2*lambda*m) + tmp*sqrt(stand$Ys2) # our estimates
    
  return(res)
}