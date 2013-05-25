\name{parplot}
\alias{parplot}
\alias{parplot.OLScurve}
\title{Plot distribution of parameters}
\usage{
  parplot(object, ...)

  \method{parplot}{OLScurve} (object, type = "hist",
    group = NULL, breaks = NULL, prompt = TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{OLScurve}}

  \item{type}{type of plot to display; can be
  \code{'hist'}, \code{'boxplot'}, or \code{'splom'} for a
  histogram, boxplot, or scatter plot matrix}

  \item{group}{a \code{factor} grouping variable used to
  partition the results}

  \item{breaks}{number of breaks to be used in plotting the
  histogram}

  \item{prompt}{a logical variable indicating whether
  \code{devAskNewPage(ask=TRUE)} should be called}

  \item{...}{additional arguments to be passed}
}
\description{
  A plotting function for displaying the distribution of
  the OLS parameter estimates.
}
\examples{
\dontrun{
data <- t(t(matrix(rnorm(1000),200)) + 1:5)
group <- rep(c('Male', 'Female'), each=nrow(data)/2)
mod <- OLScurve(~ time, data = data)
parplot(mod)
parplot(mod, type = 'boxplot')
parplot(mod, type = 'splom')

parplot(mod, group=group)
parplot(mod, type='boxplot', group=group)
parplot(mod, type='splom', group=group)
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{OLS,}
\keyword{growth}

