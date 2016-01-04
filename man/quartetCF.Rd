\name{quartetCF}
\docType{data}
\alias{quartetCF}
\title{Quartet concordance factors across Arabidopsis thaliana}
\description{Concordance factors of quartets for 29 A. thaliana accessions and A. lyrata,
obtained from chromosome 4 using MDL to delimit loci
then % MrBayes on each locus and
BUCKy on each 4-taxon set.}

\usage{data(quartetCF)}

\format{Data frame with 7 variables and 27,405 rows.
Each row corresponds to one 4-taxon set
(choosing 4 taxa out of 30 makes 27,405 combinations).
The first four variables, named 'taxon1' through 'taxon4',
give the names of the 4 taxa in each 4-taxon set.
Variables 5 through 7 are named CF12.34, CF13.24 and CF14.23,
and give the estimated concordance factors of the 3 quartets on
each set of 4 taxa: taxon 1 + taxon 2 versus taxon3 + taxon 4, etc.
These concordance factors are the proportion of loci that have a
given quartet tree.
They were obtained from chromosome 4 using MDL to delimit loci
then BUCKy to estimate quartet concordance factors (CFs).
}

\references{
Stenz, Noah W. M., Bret Larget, David A. Baum and Cécile Ané (2015).
Exploring tree-like and non-tree-like patterns using genome sequences:
An example using the inbreeding plant species \emph{Arabidopsis thaliana} (L.) Heynh.
Systematic Biology, \bold{64}(5):809-823.

TICR pipeline: github.com/nstenz/TICR
}
\keyword{datasets}