\name{pruningwise.distFromRoot}
\alias{pruningwise.distFromRoot}
\title{Calculates node distance from the root in
  an "pruningwise" tree}
\description{Calculates the distance from the root to all nodes, in
  a tree whose internal representation is in "pruningwise" order.}
\usage{
pruningwise.distFromRoot(phy)
}
\arguments{
  \item{phy}{a phylogenetic tree of type phylo with branch
    lengths, already in "pruningwise" order.}
}
\value{
  {a vector of distances, with the original tip labels and internal node names if
    internal node names were available, or otherwise named by the node numbers
    in \code{phy}.}
}
\author{Lam Si Tung Ho}

\seealso{
\code{\link{pruningwise.branching.times}}, \code{\link[stats]{cophenetic}}.
}

\examples{
tre = reorder(rtree(50),"pruningwise")
pruningwise.distFromRoot(tre)
}
