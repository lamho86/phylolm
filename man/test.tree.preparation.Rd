\encoding{latin1}
\name{test.tree.preparation}
\alias{test.tree.preparation}

\title{data structure preparation for testing a population tree}
\description{Takes a guide tree and quartet concordance factor data,
 and makes preliminary calculations to speed up the test of adequacy
 of a population tree with \code{\link{test.one.species.tree}}.}
\usage{
test.tree.preparation(cf, guidetree)
}

\arguments{
  \item{cf}{data frame containing one row for each 4-taxon set, with taxon names
   in columns 1-4.}
  \item{guidetree}{tree of class phylo on the same taxon set as those in \code{cf}.}
}
\value{
  \item{quartet2edge}{matrix of booleans with as many rows as in \code{cf}
   (one row for each 4-taxon set) and with as many columns as edges in \code{guidetree}.
   For a given set of 4 taxa and for edge \code{i}, the corresponding entry is true
   if the guide tree pruned to the 4 taxa includes edge \code{i} as part of its internal edge.}
  \item{dominant}{Vector with as many entries as rows in \code{cf}, one for each 4-taxon set.
   Entries are 1, 2 or 3 depending on the quartet displayed in the guide tree: 1 for the quartet
   12|34 (i.e. taxon1 + taxon2 versus taxon3 + taxon4), 2 for the quartet 13|24 and 3 for the
   quartet 14|23.}
}
\seealso{\code{\link{test.one.species.tree}}.}