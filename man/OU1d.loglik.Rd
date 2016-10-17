\name{OU1d.loglik}
\alias{OU1d.loglik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log likelihood of an one-dimensional Ornstein-Uhlenbeck model
}
\description{
computes log likelihood of an one-dimensional Ornstein-Uhlenbeck model with an algorithm that is linear in the number of tips in the tree. 
}
\usage{
OU1d.loglik(trait, phy, model = c("OUrandomRoot", "OUfixedRoot"), parameters = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trait}{
a vector of trait values.
}
  \item{phy}{
a phylogenetic tree of type phylo with branch lengths.
}
  \item{model}{
an Ornstein-Uhlenbeck model.
}
  \item{parameters}{
List of parameters for the model
}
}
\author{Lam Si Tung Ho}
\examples{
tr = rtree(100)
alpha = 1
sigma2 = 1
sigma2_error = 0.5
ancestral.state = 0
optimal.value = 1
  
trait = rTrait(n = 1, tr, model = "OU", 
              parameters = list(ancestral.state=ancestral.state, alpha=alpha,
                                sigma2=sigma2,sigma2_error=sigma2_error,
                                optimal.value=optimal.value))
OU1d.loglik(trait=trait, phy=tr, model="OUfixedRoot", 
            parameters=list(ancestral.state=ancestral.state, alpha=alpha,sigma2=sigma2,
                            sigma2_error=sigma2_error,optimal.value=optimal.value))
}
