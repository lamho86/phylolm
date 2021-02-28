\encoding{latin1}
\name{phylostep}
\alias{phylostep}

\title{Stepwise model selection for Phylogenetic Linear Model}
\description{Performs stepwise model selection for phylogenetic linear
  models, using the criterion -2*log-likelihood + k*npar, where npar is
  the number of estimated parameters and k=2 for the usual AIC.}
\usage{
phylostep(formula, starting.formula = NULL, keeping.formula = NULL, data = list(), 
       phy, model = c("BM", "OUrandomRoot","OUfixedRoot", 
       "lambda", "kappa", "delta", "EB", "trend"),
       direction = c("both", "backward", "forward"), trace = 2,
       lower.bound = NULL, upper.bound = NULL, 
       starting.value = NULL, k=2, ...)
}

\arguments{
  \item{formula}{formula of the full model.}
  \item{starting.formula}{optional formula of the starting model.}
  \item{keeping.formula}{optional formula of the keeping model. Covariates of the keeping model are always included in the model.}
  \item{data}{a data frame containing variables in the model. If not
  found in \code{data}, the variables are taken from current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{a model for the phylogenetic covariance of residuals.}
  \item{direction}{direction for stepwise search, can be \code{both}, \code{forward}, and \code{backward}.}
  \item{trace}{if positive, information on each searching step is
    printed. Larger values may give more detailed information.}
  \item{lower.bound}{optional lower bound for the optimization of the phylogenetic model parameter.}
  \item{upper.bound}{optional upper bound for the optimization of the phylogenetic model parameter.}
  \item{starting.value}{optional starting value for the optimization of
    the phylogenetic model parameter.}
  \item{k}{optional weight for the penalty.}
  \item{\dots}{further arguments to be passed to the function \code{optim}.}
}
\details{The default \eqn{k=2} corresponds to the usual AIC penalty.
  Use \eqn{k=\log(n)}{k=log(n)} for the usual BIC, although it is
  unclear how BIC should be defined for phylogenetic regression.
  
  See \code{\link{phylolm}} for details on the possible
  phylogenetic models for the error term, for default bounds on the
  phylogenetic signal parameters, or for matching tip labels between the
  tree and the data.

}
\value{
   A phylolm object correponding to the best model is returned.
}

\author{Lam Si Tung Ho and Cecile Ane}
\seealso{
\code{\link{phylolm}}.
}

\examples{
set.seed(123456)
tre = rcoal(60)
taxa = sort(tre$tip.label)
b0=0; b1=1;
x1 = rTrait(phy=tre,model="BM",
           parameters=list(ancestral.state=0,sigma2=10))
x2 = rTrait(phy=tre,model="BM",
            parameters=list(ancestral.state=0,sigma2=10))
x3 = rTrait(phy=tre,model="BM",
            parameters=list(ancestral.state=0,sigma2=10))
y <- b0 + b1*x1 + 
     rTrait(n=1,phy=tre,model="BM",parameters=list(
              ancestral.state=0,sigma2=1))
dat = data.frame(trait=y[taxa],pred1=x1[taxa],pred2=x2[taxa],pred3=x3[taxa])
fit = phylostep(trait~pred1+pred2+pred3,data=dat,phy=tre,model="BM",direction="both")
summary(fit)
}


