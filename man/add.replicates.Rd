\encoding{latin1}
\name{add.replicates}
\alias{add.replicates}
\title{Add replicates as tips of the tree}
\description{
Create a new tree such that each replicate is a tip steming from its parent node
that is a spcies of the original tree, with a branch lenght of zero.
This aims to represent multiple measurement from a given species.}
\usage{
add.replicates(phy, data,
    species_id = "species", sample_id = "sample",
    eps = .Machine$double.eps^2)
}
\arguments{
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{data}{either a vector containing replicated species names, or a data frame containing at least two columns, one with sample ids named \code{sample_id}, and one with species names for each samples, named \code{species_id}.}
  \item{species_id}{optional name of the column in \code{data} containing species names. Default to \code{"species"}.}
  \item{sample_id}{optional name of the column in \code{data} containing sample ids. Default to \code{"sample"}.}
  \item{eps}{length of terminal branches, that should be close to zero. Default to \code{.Machine$double.eps^2}.}
}
\details{
Note that this tree can only be used with models that include within-species variation
(e.g., "lambda" or "GC") or with \code{measurement_error=TRUE}.

This function requires function \code{\link[phytools]{bind.tip}} from package \code{phytools}.
}
\value{
  a rooted tree with one tip for each replicate in \code{data}.
}
\author{Paul Bastide}
\note{
The new tip branch lenghts \code{eps} should be close to zero,
but not strictly equal to zero to avoid numerical instabilities.
}

\seealso{
\code{\link[phylolm]{phylolm}}, \code{\link[phytools]{bind.tip}}.
}

\examples{
## GC model example on the wild tomato dataset
data(tomato)
# tree without the replicates
plot(tomato$phy)
# create a tree with individual measurements as tips
tree_rep <- add.replicates(tomato$phy,tomato$data,species_id="Sp_AccessionID",sample_id="SampleID")
plot(tree_rep)
# fit GC model with measurement error on corolla diameter
fit <- phylolm(CD~1,data=tomato$data,phy=tree_rep,model="GC",measurement_error=TRUE)
summary(fit)
}
