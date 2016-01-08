\name{stepwise.test.tree}
\alias{stepwise.test.tree}

\title{Fits a population tree to data from quartet concordance factors}
\description{From a set of quartet concordance factors obtained from genetic data
 (proportion of loci that truly have a given quartet) and from a guide tree,
 this functions uses a stepwise search to find the best resolution of that guide tree.
 Any unresolved edge corresponds to ancestral panmixia, on which the coalescent process
 is assumed.}
\usage{
stepwise.test.tree(cf, guidetree, search="both", method="PLL", kbest=5,
                   maxiter=100, startT="panmixia", shape.correction=TRUE)
}

\arguments{
  \item{cf}{data frame containing one row for each 4-taxon set and containing taxon names
  in columns 1-4, and concordance factors in columns 5-7.}
  \item{guidetree}{tree of class phylo on the same taxon set as those in \code{cf},
   with branch lengths in coalescent units.}
  \item{search}{one of "both" (stepwise search both forwards and backwards at each step),
   or "heuristic" (heuristic shallow search: not recommended). }
  \item{method}{Only "PLL" is implemented. The scoring criterion to rank population trees is the
   pseudo log-likelihood (ignored if search="heuristic").}
  \item{kbest}{Number of candidate population trees to consider at each step
   for the forward and for the backward phase (separately). Use a lower value for faster
   but less thorough search.}
  \item{maxiter}{Maximum number of iterations. One iteration consists of considering
   multiple candidate population trees, using both a forward step and a backward step.}
  \item{startT}{starting population tree. One of "panmixia", "fulltree", or a numeric
   vector of edge numbers to keep resolved. The other edges are collapsed for panmixia.}
  \item{shape.correction}{boolean. If true, the shapes of all Dirichlet distributions
   used to test the adequacy of a population tree
   are corrected to be greater or equal to 1. This correction avoids Dirichlet densities
   going near 0 or 1. It is applied both when the \eqn{\alpha}{alpha} parameter is estimated and
   when the outlier p-values are calculated.}
}
\value{
  \item{Nedge}{Number of edges kept resolved in the guide tree. Other edges
   are collapsed to model ancestral panmixia.}
  \item{edges}{Indices of edges kept resolved in the guide tree.}
  \item{notincluded}{Indices of edges collapsed in the guide tree, to model ancestral panmixia.}
  \item{alpha}{estimated \eqn{\alpha}{alpha} parameter.}
  \item{negPseudoLoglik}{Negative pseudo log-likelihood of the final estimated population tree.}
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
\seealso{\code{\link{test.one.species.tree}}.}

\examples{
data(quartetCF)
data(guidetree)
\donttest{resF <- stepwise.test.tree(quartetCF,guidetree,startT="fulltree") # takes ~ 1 min
resF[1:9]
}
}
