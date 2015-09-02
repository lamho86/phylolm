\encoding{latin1}
\name{bootstrap}
\alias{bootstrap}

\title{Bootstrap Method for Binary Traits}
\description{Runs phylogenetic logistic regression using \code{phyloglm}
  and a parametric bootstrap analysis. Binary trait simulation of
  bootstrap replicates is done with \code{rBinTrait} under the model
  estimated from the original data.}
\usage{
bootstrap(formula, data, phy, method = c("logistic_MPLE",
          "logistic_IG10"), btol = 10, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL, B = 100, full.matrix = FALSE)
}

\arguments{
  \item{formula}{a model formula.}
  \item{data}{a data frame containing variables in the model. If not
  found in \code{data}, the variables are taken from the current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{method}{The "logistic_IG10" method optimizes a GEE approximation
    to the penalized likelihood of the logistic regression.
    "logistic_MPLE" maximizes the penalized likelihood of the logistic
    regression. In both cases, the penalty is Firth's correction.} 
  \item{btol}{bound on the linear predictor to bound the searching space.}
  \item{log.alpha.bound}{bound for the log of the parameter alpha}
  \item{start.beta}{starting values for beta coefficients}
  \item{start.alpha}{starting values for alpha (phylogenetic correlation)}
  \item{B}{number of independent bootstrap replicates}
  \item{full.matrix}{if \code{TRUE}, the full matrix of bootstrap
  estimates (coefficients and alpha) will be returned.}
}
\value{
  \item{mean}{means of the parameters estimated from bootstrap replicates.}
  \item{sd}{standard deviations of the estimated parameters.}
  \item{confint95}{95\% confidence interval.}
  \item{log.alpha.mean}{mean of the logs of the estimated alphas.}
  \item{log.alpha.sd}{standard deviation of the logs of the estimated alphas.}
  \item{numFailed}{number of independent bootstrap replicates for which
  \code{phyloglm} failed. These failures may be due to the bootstrap
  data having too few 0's or too few 1's.}
  \item{bootstrap}{( \code{full.matrix} = \code{TRUE} only) matrix of
  all bootstrap estimates.}
}
\references{
Ho, L. S. T. and Ané, C. 2014. "A linear-time algorithm for Gaussian and
  non-Gaussian trait evolution models".  Systematic Biology \bold{63}(3):397-408.

Ives, A. R. and T. Garland, Jr. 2010. "Phylogenetic logistic regression
for binary dependent variables".  Systematic Biology \bold{59}:9-26.

Paradis E. and Claude J. 2002. "Analysis of Comparative Data Using Generalized Estimating Equations". Journal of Theoretical Biology \bold{218}:175-185.
}
\author{Rachel Feldman and Cécile Ané}
\seealso{
\code{\link{phyloglm}}.
}

\examples{
set.seed(123456)
tre <- rtree(50)
x <- rTrait(n=1,phy=tre)
X <- cbind(rep(1,50),x)
y <- rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
dat <- data.frame(trait01=y, predictor=x)
fit <- bootstrap(formula = trait01~predictor, data=dat, phy=tre,
                 method = "logistic_MPLE", full.matrix=TRUE)
fit$confint95
fit$bootstrap
}
