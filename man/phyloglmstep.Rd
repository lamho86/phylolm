\encoding{latin1}
\name{phyloglmstep}
\alias{phyloglmstep}

\title{Stepwise model selection for Phylogenetic Generalized Linear Model}
\description{Performs stepwise model selection for phylogenetic generalized
  linear models, using the criterion -2*log-likelihood + k*npar, where npar 
  is the number of estimated parameters and k=2 for the usual AIC.}
\usage{
phyloglmstep(formula, starting.formula = NULL, data=list(), phy, 
       method = c("logistic_MPLE","logistic_IG10"),
       direction = c("both", "backward", "forward"), trace = 2,
       btol = 10, log.alpha.bound = 4, start.beta=NULL, 
       start.alpha=NULL, boot = 0, full.matrix = TRUE,
       k=2, ...)
}

\arguments{
  \item{formula}{formula of the full model.}
  \item{starting.formula}{optional formula of the starting model.}
  \item{data}{a data frame containing variables in the model. If not
  found in \code{data}, the variables are taken from current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{method}{The "logistic_IG10" method optimizes a GEE approximation to the penalized likelihood of the logistic regression. "logistic_MPLE" maximizes the penalized likelihood of the logistic regression. In both cases, the penalty is Firth's correction.}
  \item{direction}{direction for stepwise search, can be \code{both}, \code{forward}, and \code{backward}.}
  \item{trace}{if positive, information on each searching step is
    printed. Larger values may give more detailed information.}
  \item{btol}{bound on the linear predictor to bound the searching space.}
  \item{log.alpha.bound}{bound for the log of the parameter alpha.}
  \item{start.beta}{starting values for beta coefficients.}
  \item{start.alpha}{starting values for alpha (phylogenetic correlation).}
  \item{boot}{number of independent bootstrap replicates, \code{0} means no bootstrap.}
  \item{full.matrix}{if \code{TRUE}, the full matrix of bootstrap estimates (coefficients and alpha) will be returned.}
  \item{k}{optional weight for the penalty.}
  \item{\dots}{further arguments to be passed to the function \code{optim}.}
}
\details{The default \eqn{k=2} corresponds to the usual AIC penalty.
  Use \eqn{k=\log(n)}{k=log(n)} for the usual BIC, although it is
  unclear how BIC should be defined for phylogenetic regression.
  
  See \code{\link{phyloglm}} for details on the possible
  phylogenetic methods for the error term, for default bounds on the
  phylogenetic signal parameters, or for matching tip labels between the
  tree and the data.

}
\value{
   A phyloglm object correponding to the best model is returned.
}

\author{Rutger Vos}
\seealso{
\code{\link{phyloglm}}.
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
X = cbind(rep(1,60), x1)
y = rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
dat = data.frame(trait=y[taxa],pred1=x1[taxa],pred2=x2[taxa],pred3=x3[taxa])
fit = phyloglmstep(trait~pred1+pred2+pred3,data=dat,phy=tre,method="logistic_MPLE",direction="both")
summary(fit)
}

