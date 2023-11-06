\name{phylolm-methods}
\alias{print.phylolm}
\alias{summary.phylolm}
\alias{print.summary.phylolm}
\alias{residuals.phylolm}
\alias{predict.phylolm}
\alias{vcov.phylolm}
\alias{logLik.phylolm}
\alias{print.logLik.phylolm}
\alias{AIC.phylolm}
\alias{AIC.logLik.phylolm}
\alias{extractAIC.phylolm}
\alias{nobs.phylolm}
\alias{plot.phylolm}
\alias{confint.phylolm}
\alias{model.frame.phylolm}

\title{Methods for class 'phylolm'.}
\description{These are method functions for class 'phylolm'.}
\usage{
\method{print}{phylolm}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{phylolm}(object, ...)
\method{nobs}{phylolm}(object, ...)
\method{residuals}{phylolm}(object,type = c("response"), ...)
\method{predict}{phylolm}(object, newdata = NULL, se.fit = FALSE, ...)
\method{vcov}{phylolm}(object, ...)
\method{logLik}{phylolm}(object, ...)
\method{AIC}{phylolm}(object, k=2, ...)
\method{plot}{phylolm}(x, ...)
\method{confint}{phylolm}(object, parm, level=0.95, ...)
\method{model.frame}{phylolm}(formula, ...)
}

\arguments{    
  \item{x}{an object of class \code{"phylolm"}.}
  \item{object}{an object of class \code{"phylolm"}.}
  \item{formula}{an object of class \code{"phylolm"}.}
  \item{digits}{number of digits to show in summary method.}
  \item{type}{Currently, only the "response" type is implemented. It returns the raw residuals, that is, the differences between the observed responses and the predicted values. They are phylogenetically correlated.}
  \item{newdata}{an optional data frame to provide the predictor values
  at which predictions should be made. If omitted, the fitted values are used. Currently, predictions are made for new species whose placement in the tree is unknown. Only their covariate information is used. The prediction for the trend model is not currently implemented.}
  \item{se.fit}{A switch indicating if standard errors are required.}
  \item{k}{numeric, the penalty per parameter to be used; the default \code{k = 2} is the classical AIC.}
  \item{parm}{a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
  \item{\dots}{further arguments to methods.}
}

\author{Lam Si Tung Ho}
\seealso{
\code{\link{phylolm}}
}

\examples{
set.seed(321123)
tre = rcoal(50)
y = rTrait(n=1,phy=tre,model="BM")
fit = phylolm(y~1,phy=tre,model="BM")
summary(fit)
vcov(fit)
}


