\encoding{latin1}
\name{add.individuals}
\alias{add.individuals}
\title{Add individuals as tips of the tree}
\description{
Create a new tree such that each individual within a species (or population)
is a tip stemming from its parent node
that is a species of the original tree, with a branch of length (or close to) zero.
This can be used if there are several measurements from a single species.}
\usage{
add.individuals(phy, data,
    species_id = "species", sample_id = "sample",
    eps = .Machine$double.eps^2)
}
\arguments{
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{data}{either a vector containing replicated species names, or a data frame containing at least two columns, one with unique sample identifiers named \code{sample_id}, and one with species names for each sample, named \code{species_id}.}
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
  a rooted tree with one tip for each individual in \code{data}.
}
\author{Paul Bastide}
\note{
For some applications, the new external (tip) branch lengths \code{eps} should be
close to zero, but not strictly equal to zero to avoid numerical instabilities.
}

\seealso{
\code{\link[phylolm]{phylolm}}, \code{\link[phytools]{bind.tip}}.
}

\examples{
## Adding individuals on a tree from a vector
# simulate species tree
set.seed(1289)
ntips <- 10
tree <- ape::rphylo(ntips, 0.1, 0)
plot(tree)
# vector specifying how many individuals there are per species
reps <- sample(1:3, ntips, replace=TRUE)
species_vec <- rep(tree$tip.label, times=reps)
# add individuals on the tree
tree_rep <- add.individuals(tree, species_vec)
plot(tree_rep)

## GC model example on the wild tomato dataset
data(tomato)
# species tree
plot(tomato$phy)
# create a tree with individuals as tips using the data frame
tree_rep <- add.individuals(tomato$phy, tomato$data,
                            species_id="Sp_AccessionID", sample_id="SampleID")
plot(tree_rep)
# fit GC model with measurement error on corolla diameter
fit <- phylolm(CD~1, data=tomato$data, phy=tree_rep, model="GC", measurement_error=TRUE)
summary(fit)
}
