\encoding{latin1}
\name{rbinTrait}
\alias{rbinTrait}
\title{Binary trait simulation}
\description{Simulates a binary trait along a phylogeny, according to the model in Ives and Garland (2010).}
\usage{
rbinTrait(n=1, phy, beta, alpha, X = NULL, model = c("LogReg"))
}
\arguments{
  \item{n}{number of independent replicates}
  \item{phy}{a phylogenetic tree of type phylo with brach lengths.}
  \item{beta}{a vector of coefficients for the logistic regression model.}
  \item{alpha}{the phylogenetic correlation parameter.}
  \item{X}{a design matrix with one row for each tip in the phylogenetic tree}
  \item{model}{Currently, only phylogenetic logistic regression is implemented.}
}

\value{If \code{n=1}, a numeric vector of 0-1 values with names from the
  tip labels in the tree. For more than 1 replicate, a matrix with the
  tip labels as row names, and one column per replicate.}

\note{In the absence of a design matrix \code{X}, a single intercept is
  used. In this case \code{beta} should be a vector of length one and
  the model reduces to a 2-state Markov process on the tree with
  stationary mean \eqn{e^\beta/(1+e^\beta)}{exp(beta)/(1+exp(beta))}.
  If a design matrix \code{X} is provided, the length of \code{beta} should be
  equal to the number of columns in \code{X}.}

\author{Lam Si Tung Ho and C. Ané }

\references{
Ives, A. R. and T. Garland, Jr. 2010. "Phylogenetic logistic regression
for binary dependent variables". Systematic Biology \bold{59}:9-26.
}
\seealso{
\code{\link[phylolm]{rTrait}}.
}
\examples{
tre = rtree(50)
x = rTrait(n=1,phy=tre)
X = cbind(rep(1,50),x)
y = rbinTrait(n=1, phy=tre, beta=c(-1,0.5), alpha=1, X=X)
}
