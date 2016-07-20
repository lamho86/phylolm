library(testthat)

test_that("Log likelihood of 1d OU model",
{
  
  n = 100
  tr = rtree(n) # tree
  vcv.mat = vcv.phylo(tr)
  pat.dist.mat = cophenetic(tr)
  root.to.tip = diag(vcv.mat)
  alpha = 1
  sigma2 = 1
  sigma2_error = 0.5
  ancestral.state = 0
  optimal.value = 1
  
  trait = rTrait(n = 1, tr, model = "OU", 
                           parameters = list(ancestral.state=ancestral.state, alpha=alpha,
                                             sigma2=sigma2,sigma2_error=sigma2_error,
                                             optimal.value=optimal.value))
  root.to.tip <- diag(vcv.mat)
  X = (1-exp(-alpha * root.to.tip))*optimal.value + exp(-alpha * root.to.tip)*ancestral.state
  names(X) = names(trait)
  
  V = sigma2/2/alpha *(1 - exp(-2*alpha*vcv.mat)) * exp(-alpha * pat.dist.mat) + sigma2_error * diag(n) 
  
  loglik = as.numeric(- (n*log(2*pi) + log(det(V)) + t(trait-X)%*%solve(V)%*%(trait-X))/2)
  
  OUloglik = OU1d.loglik(trait=trait, phy=tr, model="OUfixedRoot", 
                         parameters=list(ancestral.state=ancestral.state, alpha=alpha,sigma2=sigma2,
                                         sigma2_error=sigma2_error,optimal.value=optimal.value))
  expect_equal(0.0,abs(loglik-OUloglik), 1E-10)
  
}
          )
