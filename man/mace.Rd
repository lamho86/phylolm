\name{mace}
\alias{mace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multi-task learning for ancestral state estimation.
}
\description{
Estimate the ancestral states of multiple traits simultaneously using a regularized maximum likelihood objective.
}
\usage{
mace(trait, phy, lambda = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trait}{a matrix of trait values. Each row is one species and each column is a trait.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{lambda}{regularizer parameter.}
}
\details{
  Traits are assumed to evolve according to the Brownian motion model.
}
\value{
a numeric vector of estimated ancestral states.
}
\references{
Ho, Lam Si Tung, Vu Dinh, and Cuong V. Nguyen. "Multi-task learning improves ancestral state reconstruction." Theoretical Population Biology 126 (2019): 33-39.
}
\author{
Lam Si Tung Ho
}
\note{
The default choice for \code{lambda} was proposed by Ho et al. (2019).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
m = 3
anc = c(0, 8, 16)
sig2 = c(1, 1, 2)
tree = rtree(50)

trait = rTrait(n = 1, phy = tree, model = "BM",
               parameters=list(ancestral.state = anc[1], sigma2 = sig2[1]))
for (i in 2:m) {
  trait = cbind(trait,rTrait(n = 1, phy = tree, model = "BM",
                             parameters=list(ancestral.state = anc[i], sigma2 = sig2[i])))
}
res = mace(trait, tree)
print(res)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.