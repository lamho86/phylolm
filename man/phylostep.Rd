\encoding{latin1}
\name{phylostep}
\alias{phylostep}

\title{Stepwise model selection for Phylogenetic Linear Model}
\description{Performs stepwise model selection for Phylogenetic Linear Model.}
\usage{
phylostep(formula, starting.formula = NULL, data = list(), 
       phy, model = c("BM", "OUrandomRoot","OUfixedRoot", 
       "lambda", "kappa", "delta", "EB", "trend"),
       direction = c("both", "backward", "forward"), trace = TRUE,
       lower.bound = NULL, upper.bound = NULL, 
       starting.value = NULL, ...)
}

\arguments{
  \item{formula}{formula of the full model.}
  \item{starting.formula}{formula of the starting model.}
  \item{data}{a data frame containing variables in the model. If not
  found in \code{data}, the variables are taken from current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{a model for the covariance (see Details).}
  \item{direction}{direction for stepwise search, can be \code{both}, \code{forward}, and \code{backward}.}
  \item{trace}{if \code{TRUE}, the information of each searching step is printed}
  \item{lower.bound}{optional lower bound for the optimization of the phylogenetic model parameter.}
  \item{upper.bound}{optional upper bound for the optimization of the phylogenetic model parameter.}
  \item{starting.value}{optional starting value for the optimization of the phylogenetic model parameter.}
  \item{\dots}{further arguments to be passed to the function \code{optim}.}
}
\details{This function uses an algorithm that is linear in the number of
  tips in the tree to calculate the likelihood. Possible phylogenetic
  models for the error term are the Brownian motion model (\code{BM}), the
  Ornstein-Uhlenbeck model with an ancestral state to be estimated at
  the root (\code{OUfixedRoot}), the Ornstein-Uhlenbeck model with the
  ancestral state at the root having the stationary distribution
  (\code{OUrandomRoot}), Pagel's \eqn{\lambda}{lambda} model
  (\code{lambda}), Pagel's \eqn{\kappa}{kappa} model (\code{kappa}),
  Pagel's \eqn{\delta}{delta} model (\code{delta}), the early burst
  model (\code{EB}), and the Brownian motion model with a trend
  (\code{trend}).

  By default, the bounds on the phylogenetic parameters are
  \eqn{[10^{-7}/T,50/T]} for \eqn{\alpha}{alpha}, 
  \eqn{[10^{-7},1]} for \eqn{\lambda}{lambda},
  \eqn{[10^{-6},1]} for \eqn{\kappa}{kappa},
  \eqn{[10^{-5},3]} for \eqn{\delta}{delta} and
  \eqn{[-3/T,0]} for \code{rate}, where \eqn{T} is the mean root-to-tip distance.
}
\value{
   A phylolm object correponding to the best model is returned.
}

\author{Lam Si Tung Ho}
\note{
The tip labels in the tree are matched to the data names (row names in
the data frame). If no data frame is provided through the argument data,
taxon labels in the tree are matched to taxon labels in the response
variable based on the row names of the response vector, and the taxa are
assumed to come in the same order for all variables in the model.

For the delta model, the tree is rescaled back to its original height
after each node's distance from the root is raised to the power
delta. This is to provide a stable estimate of the variance parameter
\eqn{\sigma^2}{sigma^2}. For non-ultrametric trees, the tree is rescaled
to maintain the longest distance from the root to its original value.

The trend model can only be used with non-ultrametric trees. For this
model, one predictor variable is added to the model whose values are the
distances from the root to every tip of the tree. The estimate of the
coefficent for this variable forms the trend value.
}
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


