\encoding{latin1}
\name{rTrait}
\alias{rTrait}

\title{Continuous trait simulation}
\description{Simulates a continuous trait along a tree from various phylogenetic models.}
\usage{
rTrait(n=1, phy, model=c("BM","OU","lambda","kappa","delta","EB","trend"),
       parameters = NULL, plot.tree=FALSE)
}
\arguments{
  \item{n}{number of independent replicates}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{a phylogenetic model. Default is "BM", for Brownian motion.
    Alternatives are "OU", "lambda", "kappa", "delta", "EB" and "trend".}
  \item{parameters}{List of parameters for the model (see Note).}
  \item{plot.tree}{If TRUE, the tree with transformed branch lengths
    will be shown, except for the OU model.}
}
\details{
  Possible phylogenetic models are the Brownian motion
  model (BM), the Ornstein-Uhlenbeck model (OU),
  Pagel's \eqn{\lambda}{lambda} model (lambda), Pagel's \eqn{\kappa}{kappa} model (kappa),
  Pagel's \eqn{\delta}{delta} model (delta), the early burst model (EB), and the
  Brownian motion model with a trend (trend).
}
\value{If \code{n=1}, a numeric vector with names from the tip labels in
  the tree. For more than 1 replicate, a matrix with the tip labels as
  row names, and one column per replicate.}
\note{
  The default choice for the parameters are as follows:
  \code{ancestral.state=0},
  \code{sigma2=1}, \code{optimal.value=0} for the OU model,
  \code{alpha=0} for the selection strength in the OU model,
  \code{lambda=1}, \code{kappa=1}, \code{delta=1}, \code{rate=0} for the
  EB model, \code{trend=0}. These default choices correspond to the BM model.
}

\author{Lam Si Tung Ho and C. Ane}

\seealso{
\code{\link[ape]{rTraitCont}}.
}

\examples{
tre = rtree(50)
y = rTrait(n=1, phy=tre, model="OU",
           parameters=list(optimal.value=2,sigma2=1,alpha=0.1))
}


