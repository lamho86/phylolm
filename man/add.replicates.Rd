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
  \item{data}{a data frame containing at least two columns, one with sample ids named \code{sample_id}, and one with species names named for each samples, named \code{species_id}.}
  \item{species_id}{name of the column in \code{data} containing species names. Default to \code{"species"}.}
  \item{sample_id}{name of the column in \code{data} containing sample ids. Default to \code{"sample"}.}
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
set.seed(1289)
## simulate species tree
ntips <- 10
tree <- ape::rphylo(ntips, 0.1, 0)
plot(tree)
## Add some replicates for each species
reps <- sample(1:3, ntips, replace = TRUE)
rep_ids <- make.unique(rep(tree$tip.label, times = reps), sep = "_")
## match species to samples
data <- data.frame(species = rep(tree$tip.label, times = reps),
                   sample = rep_ids)
## add replicates on the tree
tree_rep <- add.replicates(tree, data)
plot(tree_rep)
}
