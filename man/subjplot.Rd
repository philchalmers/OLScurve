\name{subjplot}
\alias{subjplot}
\alias{subjplot.OLScurve}
\title{Plot individually estimated parameters}
\usage{
  subjplot(object, ...)

  \method{subjplot}{OLScurve} (object, group = NULL,
    layout = NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{OLScurve}}

  \item{group}{a \code{factor} grouping variable used to
  parition the results}

  \item{layout}{a variable to be passed to \code{xyplot} to
  adjust the graphical layout}

  \item{...}{additional arguments to be passed}
}
\description{
  A plotting function for displaying the individuals
  trajectories and their modelled functional form. Useful
  for detecting aberrant individual trajectories.
}
\examples{
\dontrun{
data <- t(t(matrix(rnorm(1000),200)) + 1:5)
mod <- OLScurve(~ time, data = data)
subjplot(mod)


}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{growth}
\keyword{OLS,}

