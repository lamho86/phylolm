\name{test.one.species.tree}
\alias{test.one.species.tree}

\title{Tests the fit of a population tree to quartet concordance factor data}
\description{From a set of quartet concordance factors obtained from genetic data
 (proportion of loci that truly have a given quartet), this function tests
 the adequacy of the coalescent process on a given population tree, where branch
 lengths indicate coalescent units.}
\usage{
test.one.species.tree(cf, guidetree, prep, edge.keep,
                      plot=TRUE, shape.correction = TRUE)
}

\arguments{
  \item{cf}{data frame containing one row for each 4-taxon set, with taxon names
   in columns 1-4 and concordance factors in columns 5-7.}
  \item{guidetree}{tree of class phylo on the same taxon set as those in \code{cf},
   with branch lengths in coalescent units.}
  \item{prep}{result of \code{test.tree.preparation(cf,guidetree)}.
   If the test is repeated multiple times on various resolutions of the guide tree
   (see \code{edge.keep}), it saves time to do this only once.}
  \item{edge.keep}{Indices of edges to keep in the guide tree.
   All other edges are collapsed to reflect ancestral panmixia.
   In the tested population tree, the collapsed edges have length set to 0.}
  \item{plot}{boolean. If TRUE, a number of plots are output.}
  \item{shape.correction}{boolean. If TRUE, the shapes of all Dirichlet distributions
   are corrected to be greater or equal to 1. This correction avoids Dirichlet densities
   going near 0 or 1. It applies when the \eqn{\alpha}{alpha} parameter is estimated and
   when the outlier p-values are calculated.}
}
\value{
  \item{alpha}{estimated \eqn{\alpha}{alpha} parameter.}
  \item{negPseudoLoglik}{Negative pseudo log-likelihood of the population tree.}
  \item{X2}{Chi-square statistic, from comparing the counts of outlier p-values
   (in \code{outlier.table}) to the expected counts.}
  \item{chisq.pval}{p-value from the chi-square test, obtained from the comparing the \code{X2}
   value to a chi-square distribution with 3 df.}
  \item{chisq.conclusion}{character string. If the chi-square test is significant, this
   statement says if there is an excess (or deficit) of outlier 4-taxon sets.}
  \item{outlier.table}{Table with 2 rows (observed and expected counts) and 4 columns:
   number of 4-taxon sets with p-values \eqn{p\leq 0.01}{p <= 0.01},
   \eqn{0.01<p\leq 0.05}{0.01<p<=0.05}, \eqn{0.05<p\leq 0.10}{0.05<p<=0.10}
   or \eqn{p>0.10}.}
  \item{outlier.pvalues}{Vector of outlier p-values, with as many entries as there
   are rows in \code{cf}, one for each set of 4 taxa.}
  \item{cf.exp}{Matrix of concordance factors expected from the estimated population tree,
   with as many rows as in \code{cf} (one row for each 4-taxon set) and 3 columns
   (one for each of the 3 possible quartet trees).}
}
\references{
Stenz, Noah W. M., Bret Larget, David A. Baum and Cécile Ané (2015).
Exploring tree-like and non-tree-like patterns using genome sequences:
An example using the inbreeding plant species \emph{Arabidopsis thaliana} (L.) Heynh.
Systematic Biology, \bold{64}(5):809-823.
}
\author{Cécile Ané}
\seealso{\code{\link{stepwise.test.tree}}, \code{\link{test.tree.preparation}}.}

\examples{
data(quartetCF)
data(guidetree)
\donttest{prelim <- test.tree.preparation(quartetCF,guidetree) # takes 5-10 seconds

# test of panmixia: all edges collapsed, none resolved.
panmixia <- test.one.species.tree(quartetCF,guidetree,prelim,edge.keep=NULL)
panmixia[1:6]

# test of full tree: all internal edges resolved, none collapsed.
Ntaxa = length(guidetree$tip.label)
# indices of internal edges:
internal.edges = which(guidetree$edge[,2] > Ntaxa)
fulltree <- test.one.species.tree(quartetCF,guidetree,prelim,edge.keep=internal.edges)
fulltree[1:6]

# test of a partial tree, some edges (but not all) collapsed
edges2keep <- c(1,2,4,6,7,8,11,14,20,21,23,24,31,34,35,36,38,39,44,47,53)
partialTree <- test.one.species.tree(quartetCF,guidetree,prelim,edge.keep=edges2keep)
partialTree[1:5]
partialTree$outlier.table
# identify taxa most responsible for the extra outlier quartets
outlier.4taxa <- which(partialTree$outlier.pvalues < 0.01)
length(outlier.4taxa) # 483 4-taxon sets with outlier p-value below 0.01
q01 = as.matrix(quartetCF[outlier.4taxa,1:4])
sort(table(as.vector(q01)),decreasing=TRUE)
# So: Cnt_1 and Vind_1 both appear in 239 of these 483 outlier 4-taxon sets.
sum(apply(q01,1,function(x){"Cnt_1" \%in\% x | "Vind_1" \%in\% x}))
# 266 outlier 4-taxon sets have either Cnt_1 or Vind_1
sum(apply(q01,1,function(x){"Cnt_1" \%in\% x & "Vind_1" \%in\% x}))
# 212 outlier 4-taxon sets have both Cnt_1 and Vind_1
}
}
