\encoding{latin1}
\name{transf.branch.lengths}
\alias{transf.branch.lengths}
\title{Creates a tree with branch lengths to represent the 3-point structure of a covariance matrix}
\description{Creates a phylogenetic tree with branch lengths and a diagonal matrix to represent a (generalized) 3-point structured covariance matrix from a trait evolution model on a known tree.}
\usage{
transf.branch.lengths(phy, model = c("BM", "OUrandomRoot",
       "OUfixedRoot", "lambda", "kappa", "delta", "EB", "trend"),
       parameters = NULL, check.pruningwise = TRUE,
       check.ultrametric=TRUE, D=NULL, check.names = TRUE)
}
\arguments{
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{model}{a phylogenetic model. Default is "BM", for Brownian motion. Alternatives are "OU", "lambda", "kappa", "delta", "EB" and "trend".}
  \item{parameters}{List of parameters for the model (see Note).}
  \item{check.pruningwise}{if FALSE, the tree is assumed to be in pruningwise order.}
  \item{check.ultrametric}{if FALSE, the tree is assumed to be
    ultrametric and \code{D} needs to be provided. This is used for the
    OU transformations only.}
  \item{D}{vector of ajustments for the external edge lengths, to make
    the tree ultrametric. Used for the OU transformations only.}
  \item{check.names}{\code{D} needs to have names that match tip labels
    unless check.names=FALSE, in which case the elements in \code{D} are
    assumed to come in the same order as tip labels in the tree.}
}
\details{
Possible phylogenetic models are the Brownian motion model (BM), the
Ornstein-Uhlenbeck model (OU), Pagel's lambda model (lambda), Pagel's
kappa model (kappa), Pagel's delta model (delta), the early burst model
(EB), and the Brownian motion with a trend (trend). Edge lengths are
unchanged under BM and the trend model.
Under the kappa model, each branch length \eqn{\ell}{l} is transformed
to \eqn{\ell^\kappa}{l^kappa}. 
If the time from the root to a node is \eqn{t} in \code{phy},
it is transformed to
\eqn{ T * (t/T)^\delta}{T * (t/T)^delta} 
under the delta model, where \eqn{T} is the maximum root-to-tip
distance. The transformed tree has the same \eqn{T}. 
Under EB, \eqn{t} is transformed to 
\eqn{(e^{\mathrm{rate}*t}-1)/\mathrm{rate}}{(exp(rate*t)-1)/rate}, 
which is very close to \eqn{t} (i.e. to the BM model)
when \code{rate} is close to 0. 
Under the lambda model, the time \eqn{t} from the
root to a node is transformed to 
\eqn{\lambda t}{lambda * t} for an internal node and
is unchanged for a tip. 
Under "OUrandomRoot", \eqn{t} is transformed to 
\eqn{\exp(-2\alpha (T-t))}{exp(-2 alpha (T-t))}, 
where \eqn{T} is now the mean root-to-tip distance. 
Under "OUfixedRroot", \eqn{t} is transformed to 
\eqn{\exp(-2\alpha (T-t))(1-\exp(-2 \alpha t))}{exp(-2 alpha (T-t)) * (1-exp(-2 alpha t))}. 
Under the OU models, \code{diagWeight} contains \eqn{\exp(\alpha
  D_i)}{exp(alpha D_i)} for tip \eqn{i}, where \eqn{D_i} is the extra
length added to tip \eqn{i} to make the tree ultrametric.
}
\value{
  \item{tree}{a rooted tree with a root edge and transformed branch lengths.}
  \item{diagWeight}{a vector containing the diagonal elements of the
    diagonal matrix  for the generalized 3-point structure.}
}
\references{
  Ho, L. S. T. and Ane, C. \emph{A linear-time algorithm for Gaussian
    and non-Gaussian trait evolution models}. Systematic Biology \bold{63}(3):397-408.
}
\author{Lam Si Tung Ho}
\note{
The default choice for the parameters are as follows: \code{alpha=0} for
the selection strength in the OU model, \code{lambda=1}, \code{kappa=1},
\code{delta=1}, \code{rate=0} for the EB model, \code{sigma2_error=0} for the variance of measurement errors. These default choices
correspond to the BM model.

Edges in the output tree are in pruningwise order.

If \code{model="BM"} or \code{model="trend"}, the output tree is the same as the
input tree except that the output tree is in pruningwise order.
}

\seealso{
\code{\link[phylolm]{three.point.compute}}.
}

\examples{
set.seed(123456)
tre1 = rcoal(10)
tre2 = transf.branch.lengths(phy=tre1, model="OUrandomRoot",
                             parameters = list(alpha=1))
par(mfrow = c(2,1))
plot(tre1)
plot(tre2$tree,root.edge=TRUE)
}
