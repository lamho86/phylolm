\encoding{latin1}
\name{three.point.compute}
\alias{three.point.compute}
\title{Computations with a (generalized) three-point structured tree}
\description{Computes \eqn{P'V^{-1}Q} and the \eqn{\log(\det V)}{log(det
    V)} of a (generalized) three-point structured matrix.}
\usage{
three.point.compute(phy, P, Q = NULL, diagWeight = NULL, 
            check.pruningwise = TRUE, check.names = TRUE, check.precision = TRUE)
}
\arguments{
  \item{phy}{a rooted phylogenetic tree of type phylo with branch
    lengths, to represent the 3-point structured matrix \eqn{V_0}. Note
    that the matrix of interest is \eqn{V = D V_0 D}.}
  \item{P,Q}{two matrices.}
  \item{diagWeight}{a vector containing the diagonal elements of the
    diagonal matrix \eqn{D} if \eqn{V} has a generalized 3-point structure:
    \eqn{V = D V_0 D}}
  \item{check.pruningwise}{If FALSE, the tree is assumed to be in pruningwise order.}
  \item{check.names}{if FALSE, the row names of \code{P}, \code{Q}, and
    the names of \code{diagWeight} are assumed to be the same as the
    labels of the tips in the tree.} 
  \item{check.precision}{if FALSE, diagWeight will be allowed to be below Machine epsilon. Recommended to remain TRUE.} 
}
\value{
	\item{vec11}{\eqn{1'V^{-1}1}.}
	\item{P1}{\eqn{P'V^{-1}1}.}
	\item{PP}{\eqn{P'V^{-1}P}.}
	\item{Q1}{\eqn{Q'V^{-1}1}.}
	\item{QQ}{\eqn{Q'V^{-1}Q}.}
	\item{PQ}{\eqn{P'V^{-1}Q}.}
	\item{logd}{\eqn{\log(\det V)}{log(det V)}.}
}
\references{
Ho, L. S. T. and An?, C. (2014). "A linear-time algorithm for Gaussian and
non-Gaussian trait evolution models". Systematic Biology \bold{63}(3):397-408.
}
\author{Lam Si Tung Ho, Robert Lachlan}
\note{
The matrix \eqn{V} is assumed to be \eqn{V = D V_0 D} where \eqn{D} is the diagonal
matrix with non-zero diagonal elements in \code{diagWeight}, and where
\eqn{V_0} is the 3-point structured covariance matrix
determined by \code{phy} and its branch lengths. Note that \eqn{D} do
not correspond to measurement error terms. 

The number of rows in \code{P} and \code{Q} and the length of \code{diagWeight} need
to be the same as the number of tips in the tree. When Q = NULL, the
function only returns \eqn{1'V^{-1}1}, \eqn{P'V^{-1}1} and \eqn{P'V^{-1}P}.
}

\seealso{
\code{\link[phylolm]{transf.branch.lengths}}.
}

\examples{
tre1 = rtree(500)
tre2 = transf.branch.lengths(phy=tre1, model="OUrandomRoot",
                             parameters = list(alpha = 0.5))
Q = rTrait(n=2,tre1)
y = rTrait(n=1,tre1)
P = cbind(1,y)
three.point.compute(tre2$tree,P,Q,tre2$diagWeight)
}

