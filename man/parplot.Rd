\name{parplot}
\alias{parplot}
\alias{parplot.OLScurve}
\title{Plot individually estimated parameters}
\usage{
  parplot(object, ...)

  \method{parplot}{OLScurve} (object, type = "hist",
    group = NULL, breaks = NULL, ...)
}
\arguments{
  \item{object}{DESCRIPTION}

  \item{type}{DESCRIPTION}

  \item{group}{DESCRIPTION}

  \item{breaks}{DESCRIPTION}

  \item{...}{DESCRIPTION}
}
\description{
  Description of \code{parplot}
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

