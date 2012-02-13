\name{parplot}
\alias{parplot}
\alias{parplot.OLScurve}
\title{Plot distribution of parameters}
\usage{
  parplot(object, ...)

  \method{parplot}{OLScurve} (object, type = "hist",
    group = NULL, breaks = NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{OLScurve}}

  \item{type}{type of plot to display; can be
  \code{'hist'}, \code{'boxplot'}, or \code{'splom'} for a
  histogram, boxplot, or scatter plot matrix}

  \item{group}{a \code{factor} grouping variable used to
  parition the results}

  \item{breaks}{number of breaks to be used in plotting the
  histogram}

  \item{...}{additional arguments to be passed}
}
\description{
  A plotting function for displaying the distribution of
  the OLS parameter estimates.
}
\examples{
\dontrun{
data <- t(t(matrix(rnorm(1000),200)) + 1:5)
mod <- OLScurve(~ time, data = data)
parplot(mod)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{growth}
\keyword{OLS,}

