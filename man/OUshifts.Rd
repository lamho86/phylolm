\encoding{latin1}
\name{OUshifts}
\alias{OUshifts}

\title{Detections of shifts in the OU process along a phylogeny.}

\description{
  Trait data is fitted to a phylogeny using an Ornstein-Uhlenbeck (OU) process, such that the mean (or selection optimum) of the process may change in one or more edges in the tree. The number and location of changes, or shifts, is estimated using an information criterion.
}

\usage{
OUshifts(y, phy, method = c("mbic", "aic", "bic", "saic", "sbic"),
         nmax, check.pruningwise = TRUE)
}

\arguments{
  \item{y}{values for the trait data.}
  \item{phy}{a phylogenetic tree of type phylo with branch lengths.}
  \item{method}{a method for model selection (see details below).}
  \item{nmax}{maximum allowed number of shifts.}
  \item{check.pruningwise}{if TRUE, the algorithm checks if the ordering
  of the edges in phy are in pruningwise order.}
}

\details{This function does not accept multivariate data (yet): \code{y} should be a vector named with species labels.  The data \code{y} and the tree \code{phy} need to contain the same species. The user can choose among various information criteria.  Each criterion seeks to minimize the value of \eqn{-2 \log[}likelihood\eqn{(y, M)] +} penalty\eqn{(M)}, where \eqn{M} is an OU model with \eqn{m} shifts, placed on various edges along the phylogeny.  All models use \eqn{3+m} parameters: \eqn{\alpha}{alpha}, \eqn{\sigma^2}{sigma^2}, and \eqn{m+1} parameters to describe the expected trait values in each of the \eqn{ m+1} regimes.  The AIC penalty is \eqn{2*(3+m)}. The BIC penalty is  \eqn{(3+m) \log(n)}{(3+m) * \log(n)} where \eqn{n} is the numer of species.  If one considers the position of the \eqn{m} shifts in the phylogeny as parameters (even though they are discrete parameters), we get the sAIC penalty \eqn{2*(3+2m}) (used in SURFACE), and the sBIC penalty \eqn{(3+2*m)*\log(n)}.  The default penalty (model = 'mbic') is defined as \eqn{3*\log(n)+(2m-1)\log(n)+\sum_{i=0}^{m}(\log(n_i))}. A lower value of \code{nmax} will make the search faster, but if the estimated number of shifts is found equal to nmax, then the output model is probably not optimal.  Re-running with a larger \code{nmax} would take longer, but would likely return a more complex model with a better score.
}

\value{
  \item{y}{the input trait.}
  \item{phy}{the input tree.}
  \item{score}{the information criterion value of the optimal model.}
  \item{nmax}{maximum allowed number of shifts.}
  \item{nshift}{estimated number of shifts.}
  \item{alpha}{estimated selection strength of the optimal model.}
  \item{sigma2}{estimated variance of the optimal model.}
  \item{mean}{estimated the expected value of the trait in lineages without shift.}
  \item{pshift}{positions of shifts, i.e. indicies of edges where the estimated shifts occurred. The same ordering of edges is used as in phy.}
  \item{shift}{estimated shifts in the expected value of the trait.}
}

\note{
  The tip labels in the tree are matched to the data names. The default choice for the parameters are as follows:
  \code{method = "mbic"},
  \code{check.pruningwise = TRUE}
  
  Due to unidentifiability, the parameters are the expected value of the trait and their shifts instead of the ancestral trait, the optimal values and shifts in optimal values.
}

\references{
Ho, L. S. T. and Ane, C. 2014.  "Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models". Methods in Ecology and Evolution. \bold{5}(11):1133-1146.

Ingram, T. and Mahler, D.L. 2013. "SURFACE: detecting convergent
  evolution from comparative data by fitting Ornstein-Uhlenbeck models
  with step-wise Akaike information criterion". Methods in Ecology and
  Evolution \bold{4}:416-425.

Zhang, N.R. and Siegmund, D.O. 2007. "A modified Bayes information
  criterion with applications to the analysis of comparative genomic
  hybridization data". Biometrics \bold{63}:22-32.
}

\author{Lam Si Tung Ho}

\examples{
data(flowerSize)
data(flowerTree)
result <- OUshifts(flowerSize$log_transformed_size, flowerTree, 
                   method = "mbic", nmax = 1)
plot.OUshifts(result,show.tip.label=FALSE)
}


