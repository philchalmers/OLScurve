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
  \item{object}{DESCRIPTION}

  \item{group}{DESCRIPTION}

  \item{layout}{DESCRIPTION}

  \item{...}{DESCRIPTION}
}
\description{
  Description of \code{parplot}
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

