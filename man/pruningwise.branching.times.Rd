\name{pruningwise.branching.times}
\alias{pruningwise.branching.times}
\title{Calculates internal node ages in
  an ultrametric "pruningwise" tree}
\description{Calculates the branching times, or ages, of all internal nodes in
  an ultrametric tree whose internal representation is in "pruningwise" order.}
\usage{
pruningwise.branching.times(phy)
}
\arguments{
  \item{phy}{an ultrametric phylogenetic tree of type phylo with branch
    lengths, already in "pruningwise" order.}
}
\value{
  {a vector of node ages, with the original internal node names if
    those were available in \code{phy}, or otherwise named by the node numbers
    in \code{phy}.}
}
\author{Lam Si Tung Ho}

\seealso{
\code{\link{pruningwise.distFromRoot}}, \code{\link[ape]{branching.times}}.
}

\examples{
tre = reorder(rcoal(50),"pruningwise")
pruningwise.branching.times(tre)
}
