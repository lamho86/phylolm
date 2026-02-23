\name{tomato}
\docType{data}
\alias{tomato}
\title{Phylogenetic tree and trait data for 12 wild tomato populations.}
\format{
A list containing two objects:
\describe{
  \item{phy}{The phylogenetic tree for the 12 wild tomato populations, pruned to match
  the species in the trait dataset}
  \item{data}{A data frame containing individual measurements for floral trait, with columns:
  \describe{
  \item{Sp_AccessionID}{Species name and accession ID of the population.}
  \item{SampleID}{Individual ID for each sample.}
  \item{CD}{Corolla diameter (mm).}
  \item{AL}{Anther length (mm).}
  \item{SL}{Stigma length (mm).}
  }
}}
}
\usage{data(tomato)}
\description{
Dataset and phylogeny for floral traits of wild tomatoes.
}
\references{
Hibbins, Mark S., Breithaupt, Lara C. and Hahn, Matthew W. 2023.
"Phylogenomic comparative methods: Accurate evolutionary inferences in the presence of gene tree discordance."
Proceedings of the National Academy of Sciences, 120(22).
}
\keyword{datasets}
