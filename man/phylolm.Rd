\encoding{latin1}
\name{phylolm}
\alias{phylolm}

\title{Phylogenetic Linear Model}
\description{Fits a phylogenetic linear regression model. The likelihood is calculated with an algorithm that is linear in the number of tips in the tree.}
\usage{
phylolm(formula, data = list(), phy, model = c("BM", "OUrandomRoot",
       "OUfixedRoot", "lambda", "kappa", "delta", "EB", "trend"),
       lower.bound = NULL, upper.bound = NULL,
       starting.value = NULL, measurement_error = FALSE,
       boot=0,full.matrix = TRUE, ...)
}

\arguments{
  \item{formula}{a model formula.}
  \item{data}{a data frame containing variables in the model. If not
  found in \code{data}, the variables are taken from current environment.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{a model for the covariance (see Details).}
  \item{lower.bound}{optional lower bound for the optimization of the phylogenetic model parameter.}
  \item{upper.bound}{optional upper bound for the optimization of the phylogenetic model parameter.}
  \item{starting.value}{optional starting value for the optimization of the phylogenetic model parameter.}
  \item{measurement_error}{a logical value indicating whether there is measurement error \code{sigma2_error} (see Details).}
  \item{boot}{number of independent bootstrap replicates, 0 means no bootstrap.}
  \item{full.matrix}{if \code{TRUE}, the full matrix of bootstrap estimates (coefficients and covariance parameters) will be returned.}
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

  Using measurement error means that the covariance matrix is taken
  to be \eqn{\sigma^2*V + \sigma^2_{error}*I}{sigma^2 * V + sigma^2_error * I}
  where \eqn{V} is the phylogenetic covariance matrix from the chosen model,
  \eqn{I} is the identity matrix, and \eqn{\sigma^2_{error}}{sigma^2_error} is the
  variance of the measurement error (which could include environmental variability,
  sampling error on the species mean, etc.).

  By default, the bounds on the phylogenetic parameters are
  \eqn{[10^{-7}/T,50/T]} for \eqn{\alpha}{alpha}, 
  \eqn{[10^{-7},1]} for \eqn{\lambda}{lambda},
  \eqn{[10^{-6},1]} for \eqn{\kappa}{kappa},
  \eqn{[10^{-5},3]} for \eqn{\delta}{delta} and
  \eqn{[-3/T,0]} for \code{rate}, where \eqn{T} is the mean root-to-tip distance.
  \eqn{[10^{-16}, 10^{16}]} for the ratio \code{sigma2_error}/\code{sigma2} (if measurement errors is used).
  
  Bootstrapping can be parallelized using the \code{future} package on any future 
  compatible back-end. For example, run \code{library(future); plan(multiprocess))}, 
  after which bootstrapping will automatically occur in parallel. See 
  \code{\link[future]{plan}} for options.
}
\value{
   \item{coefficients}{the named vector of coefficients.}
   \item{sigma2}{the maximum likelihood estimate of the variance rate \eqn{\sigma^2}{sigma^2}.}
   \item{sigma2_error}{the maximum likelihood estimate of the variance of the measurement errors.}
   \item{optpar}{the optimized value of the phylogenetic correlation parameter (alpha, lambda, kappa, delta or rate).}
   \item{logLik}{the maximum of the log likelihood.}
   \item{p}{the number of all parameters of the model.}
   \item{aic}{AIC value of the model.}
   \item{vcov}{covariance matrix for the regression coefficients, given
  the phylogenetic correlation parameter (if any).}
   \item{fitted.values}{fitted values}
   \item{residuals}{raw residuals}
   \item{y}{response}
   \item{X}{design matrix}
   \item{n}{number of observations (tips in the tree)}
   \item{d}{number of dependent variables}
   \item{mean.tip.height}{mean root-to-tip distance, which can help
  choose good starting values for the correlation parameter.}
   \item{formula}{the model formula}
   \item{call}{the original call to the function}
   \item{model}{the phylogenetic model for the covariance}
   \item{bootmean}{(\code{boot > 0} only) bootstrap means of the parameters estimated.}
  \item{bootsd}{(\code{boot > 0} only) bootstrap standard deviations of the estimated parameters.}
  \item{bootconfint95}{(\code{boot > 0} only) bootstrap 95\% confidence interval.}
  \item{bootmeansdLog}{(\code{boot > 0} only) bootstrap mean and standard deviation of the logs of the estimated covariance parameters.}
  \item{bootnumFailed}{(\code{boot > 0} only) number of independent bootstrap replicates for which
  \code{phylolm} failed.}
  \item{bootstrap}{(\code{boot > 0} and \code{full.matrix} = \code{TRUE} only) matrix of all bootstrap estimates.}
  \item{r.squared}{The r^2 for the model.}
  \item{adj.r.squared}{The adjusted r^2 for the model.}
}
\references{
Ho, L. S. T. and Ane, C. 2014.  "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology \bold{63}(3):397-408.

Butler, M. A. and King, A. A. 2004. "Phylogenetic comparative analysis:
A modeling approach for adaptive evolution". The American Naturalist \bold{164}:683-695.

Hansen, T. F. 1997. "Stabilizing selection and the comparative analysis
of adaptation". Evolution \bold{51}:1341-1351.

Harmon, L. J. et al. 2010. "Early bursts of body size and shape
evolution are rare in comparative data". Evolution \bold{64}:2385-2396.

Ho, L. S. T. and Ane, C. 2013. "Asymptotic theory with hierarchical
autocorrelation: Ornstein-Uhlenbeck tree models". The Annals of Statistics \bold{41}(2):957-981.

Pagel, M. 1997. "Inferring evolutionary processes from
phylogenies". Zoologica Scripta \bold{26}:331-348.

Pagel, M. 1999. "Inferring the historical patterns of biological
evolution". Nature \bold{401}:877-884.

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

Pagel's \eqn{\lambda}{lambda} model and measurement error cannot be used together:
the parameters \eqn{\lambda}{lambda}, \eqn{\sigma^2}{sigma^2} and
\eqn{\sigma^2_{error}}{sigma^2_error} are not distinguishable (identifiable) from each other.
}
\seealso{
\code{\link[ape]{corBrownian}}, \code{\link[ape]{corMartins}},
\code{\link[ape]{corPagel}}, \code{\link[geiger]{fitContinuous}}, \code{\link[caper]{pgls}}.
}

\examples{
set.seed(123456)
tre = rcoal(60)
taxa = sort(tre$tip.label)
b0=0; b1=1;
x <- rTrait(n=1, phy=tre,model="BM",
            parameters=list(ancestral.state=0,sigma2=10))
y <- b0 + b1*x + 
     rTrait(n=1,phy=tre,model="lambda",parameters=list(
              ancestral.state=0,sigma2=1,lambda=0.5))
dat = data.frame(trait=y[taxa],pred=x[taxa])
fit = phylolm(trait~pred,data=dat,phy=tre,model="lambda")
summary(fit)

# adding measurement errors and bootstrap
z <- y + rnorm(60,0,1)
dat = data.frame(trait=z[taxa],pred=x[taxa])
fit = phylolm(trait~pred,data=dat,phy=tre,model="BM",measurement_error=TRUE,boot=100)
summary(fit)

}


