\encoding{latin1}
\name{phyloglm}
\alias{phyloglm}

\title{Phylogenetic Generalized Linear Model}
\description{Fits the phylogenetic logistic regression described in
  Ives and Garland (2010) and the Poisson regression described in Paradis and Claude (2002). The computation uses an algorithm that is
  linear in the number of tips in the tree.}
\usage{
phyloglm(formula, data, phy, method = c("logistic_MPLE","logistic_IG10",
         "poisson_GEE"), btol = 10, log.alpha.bound = 4,
         start.beta=NULL, start.alpha=NULL,
         boot = 0, full.matrix = TRUE)
}

\arguments{
  \item{formula}{a model formula.}
  \item{data}{a data frame containing variables in the model. If not
  found in \code{data}, the variables are taken from the current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{method}{The "logistic_IG10" method optimizes a GEE approximation to the penalized likelihood of the logistic regression. "logistic_MPLE" maximizes the penalized likelihood of the logistic regression. In both cases, the penalty is Firth's correction. The "poisson_GEE" method solves the generalized estimating equations (GEE) for Poisson regression.}
  \item{btol}{(logistic regression only) bound on the linear predictor to bound the searching space.}
  \item{log.alpha.bound}{(logistic regression only) bound for the log of the parameter alpha.}
  \item{start.beta}{starting values for beta coefficients.}
  \item{start.alpha}{(logistic regression only) starting values for alpha (phylogenetic correlation).}
  \item{boot}{number of independent bootstrap replicates, \code{0} means no bootstrap.}
  \item{full.matrix}{if \code{TRUE}, the full matrix of bootstrap estimates (coefficients and alpha) will be returned.}
}
\details{
  This function uses an algorithm that is linear in the number of tips in the tree.

  Bootstrapping can be parallelized using the \code{future} package on any future 
  compatible back-end. For example, run \code{library(future); plan(multiprocess))}, 
  after which bootstrapping will automatically occur in parallel. See 
  \code{\link[future]{plan}} for options.
}
\value{
  \item{coefficients}{the named vector of coefficients.}
  \item{alpha}{(logistic regression only) the phylogenetic correlation parameter.}
  \item{scale}{(Poisson regression only) the scale parameter which estimates the overdispersion.}
  \item{sd}{standard deviation for the regression coefficients.}
  \item{vcov}{covariance matrix for the regression coefficients.}
  \item{logLik}{(logistic regression only) log likelihood.}
  \item{aic}{(logistic regression only) AIC.}
  \item{penlogLik}{(logistic regression only) penalized log likelihood, using Firth's penalty for coefficients.}
  \item{y}{response.}
  \item{n}{number of observations (tips in the tree).}
  \item{d}{number of dependent variables.}
  \item{formula}{the model formula.}
  \item{call}{the original call to the function.}
  \item{method}{the estimation method.}
  \item{convergence}{An integer code with '0' for successful optimization. With logistic_MPLE, this is the convergence code from the \code{optim} routine.}
  \item{X}{a design matrix with one row for each tip in the phylogenetic tree.}
  \item{bootmean}{(\code{boot > 0} only) bootstrap means of the parameters estimated.}
  \item{bootsd}{(\code{boot > 0} only) bootstrap standard deviations of the estimated parameters.}
  \item{bootconfint95}{(\code{boot > 0} only) bootstrap 95\% confidence interval.}
  \item{bootmeanAlog}{(\code{boot > 0} only) bootstrap mean of the logs of the estimated alphas.}
  \item{bootsdAlog}{(\code{boot > 0} only) bootstrap standard deviation of the logs of the estimated alphas.}
  \item{bootnumFailed}{(\code{boot > 0} only) number of independent bootstrap replicates for which
  \code{phyloglm} failed. These failures may be due to the bootstrap data having too few 0's or too few 1's.}
  \item{bootstrap}{(\code{boot > 0} and \code{full.matrix} = \code{TRUE} only) matrix of all bootstrap estimates.}
}
\references{
Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for Gaussian and
  non-Gaussian trait evolution models".  Systematic Biology \bold{63}(3):397-408.

Ives, A. R. and T. Garland, Jr. 2010. "Phylogenetic logistic regression
for binary dependent variables".  Systematic Biology \bold{59}:9-26.

Paradis E. and Claude J. 2002. "Analysis of Comparative Data Using Generalized Estimating Equations". Journal of Theoretical Biology \bold{218}:175-185.
}
\author{Lam Si Tung Ho, Robert Lachlan, Rachel Feldman and C?cile An?}
\note{
The tip labels in the tree are matched to the data names (row names in
the data frame). If no data frame is provided through the argument \code{data},
taxon labels in the tree are matched to taxon labels in the response
variable based on the row names of the response vector, and the taxa are
assumed to come in the same order for all variables in the model.

The logistic regression method of Ives and Garland (2010) uses alpha to estimate the level of phylogenetic correlation. The GEE method for Poisson regression does not estimate the level of phylogenetic correlation but takes it from the existing branch lengths in the tree.

The standard deviation and the covariance matrix for the coefficients of logistic regression are conditional on the estimated value of the phylogenetic correlation parameter \eqn{\alpha}{alpha}.

The default choice \code{btol=10} constrains the fitted values, i.e. the probability of "1" predicted by the model, to lie within
\eqn{1/(1+e^{ 10})=0.000045}{1/(1+exp(+10)) = 0.000045} and
\eqn{1/(1+e^{-10})=0.999955}{1/(1+exp(-10)) = 0.999955}.

The log of \eqn{\alpha}{alpha} is bounded in the interval
\eqn{-\log(T) \pm \mathrm{log.alpha.bound}}{(-log(T) - log.alpha.bound,
  -log(T) + log.alpha.bound)} 
where \eqn{T} is the mean of the distances from the root to tips. In
other words, \eqn{\alpha T}{alpha*T} is constrained to lie within
\eqn{\exp(\pm\mathrm{log.alpha.bound})}{exp(+/- log.alpha.bound)}.
}
\seealso{
\code{\link[ape]{compar.gee}}.
}

\examples{
set.seed(123456)
tre = rtree(50)
x = rTrait(n=1,phy=tre)
X = cbind(rep(1,50),x)
y = rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
dat = data.frame(trait01 = y, predictor = x)
fit = phyloglm(trait01~predictor,phy=tre,data=dat,boot=100)
summary(fit)
coef(fit)
vcov(fit)
}
